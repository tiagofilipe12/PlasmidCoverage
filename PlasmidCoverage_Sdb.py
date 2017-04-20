#!/usr/bin/env python

## Last update: 19/4/2017
## Author: T.F. Jesus
## This version runs with bowtie2 build 2.2.9, with multithreading for bowtie2-build

import argparse
import os
import re
import operator
import plotly
import plotly.graph_objs as go
from subprocess import Popen, PIPE, call
from Bio import SeqIO
from time import time
from datetime import datetime
import json


# TODO:
# 1) Set a minnimun read threshold to avoid computation on less than x reads mapped (depth_file?)
# 4) delete unsorted BAM files
# 5) Create directory structure to store all comparisons.

def alignmaxnumber(max_align, count_entries):
    if max_align:
        k_value = max_align
    else:
        k_value = str(count_entries)
    return k_value


def FolderExist(directory):
    if not directory.endswith("/"):
        directory = directory + "/"
    if not os.path.exists(os.path.join(directory)):
        os.makedirs(os.path.join(directory))
        print
        os.path.join(directory) + " does not exist. One will be created..."
    else:
        print
        os.path.join(directory) + " exists!"


def fastadict(fasta_file):
    if_handle = open(fasta_file, 'r')
    x = 0
    fasta_dic = {}
    problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/", "[", "]", ":", "{", "}"]
    sequence_list = []
    for line in if_handle:
        if len(line) > 0:
            line = line.splitlines()[0]
        if x == 0 and line.startswith(">"):
            sequence_list = []
            PlasmidName = "_".join(line[1:].split("|")[0:2])  ## stores only the gi for the reference
            for char in problematic_characters:
                PlasmidName = PlasmidName.replace(char, '_')
            x += 1
        elif x == 0 and not line.startswith(">"):
            print
            "Is this a fasta file? " + fasta_file
            print
            fasta_file + " will be ignored"
            break
        elif x >= 1 and line.startswith(">"):
            fasta_dic[
                PlasmidName] = sequence_list  # appends last sequence to be parsed before new structure for sequence
            sequence_list = []
            PlasmidName = "_".join(line[1:].split("|")[0:2])
            for char in problematic_characters:
                PlasmidName = PlasmidName.replace(char, '_')  ## stores only the gi for the reference
            x += 1
        else:
            sequence_list.append(line)
    if sequence_list:
        fasta_dic[PlasmidName] = sequence_list  # appends last sequence on the fasta
    if_handle.close()
    return fasta_dic


def SequenceLengthFromFasta(fasta_file, plasmid_length, fasta_path):
    fasta_dic = fastadict(fasta_file)
    print
    len(fasta_dic.keys())
    out_handle = open(os.path.join(fasta_path + ".temp"), 'w')
    for key in fasta_dic:
        key = "_".join(key.split("_")[0:2])  ## stores only the gi for the reference
        plasmid_length[key] = sum(len(s) for s in fasta_dic[key])
        out_handle.write('>' + key + '\n' + ''.join(fasta_dic[key]) + '\n')
    out_handle.close()
    return plasmid_length, len(fasta_dic.keys())


def ExtractFastaPlasmids(gbkfile, fastafile, plasmid_length):
    if_handle = open(gbkfile, 'r')
    gbkdata = SeqIO.read(if_handle, "genbank")
    problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/", "[", "]", ":", "{", "}"]
    plasmid_name = gbkdata.description[:-1]
    for char in problematic_characters:
        plasmid_name = plasmid_name.replace(char, '_')
    sequence = str(gbkdata.seq)
    if not sequence.endswith("\n"):
        sequence = sequence + "\n"
    out_handle = open(fastafile, 'w')
    out_handle.write('>' + plasmid_name + '\n')
    out_handle.write(sequence)
    if_handle.close()
    out_handle.close()
    plasmid_length[plasmid_name] = len(str(gbkdata.seq))
    print
    plasmid_length.values(), " values each gb file"
    print(" Wrote fasta file: " + fastafile)
    return plasmid_length


def CreateBowtieIdx(filename, dirname, threads):
    check = ""
    if not os.path.exists(os.path.join(dirname + "bowtie2idx")):
        os.makedirs(os.path.join(dirname + "bowtie2idx"))
    else:
        pass
    idx_file = os.path.join(dirname, 'bowtie2idx', os.path.splitext(filename)[0]) + '.idx'
    fasta_file = os.path.join(dirname, 'fasta', os.path.splitext(filename)[0]) + '.fasta'
    list_idxfiles = [idx_file + ".1.bt2", idx_file + ".2.bt2", idx_file + ".3.bt2", idx_file + ".4.bt2",
                     idx_file + ".rev.2.bt2", idx_file + ".rev.1.bt2"]
    for idx in list_idxfiles:
        if not os.path.isfile(idx):
            check = "yes"
    if "yes" in check:
        print("Creating " + idx_file)
        # Create bowtie index
        call('bowtie2-build -q ' + fasta_file + ' --threads ' + threads + ' ' + idx_file,
             shell=True)  # convert to popen
    else:
        print
        idx_file + " already exists!"

    return idx_file


def FastaConcatenation(dblist, output_name, plasmid_dir):
    print(dblist)
    main_filename = output_name + ".fasta"
    dirname = os.path.join(plasmid_dir, "fasta", "")
    if os.path.isfile(dirname + main_filename):
        print
        output_name + ".fasta already exists. Overriding file..."
    print
    "Saving to: " + main_filename
    concat = Popen("cat " + ' '.join(dblist) + "> " + dirname + main_filename, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = concat.communicate()
    print
    stderr
    return main_filename


## function to delete temporary fasta files
def deltemp(directory):
    files = os.listdir(directory)
    print
    "Deleting temporary fasta files in: " + directory
    for f in files:
        if f.endswith('.temp'):
            os.remove(os.path.join(directory, f))


## Creates an alignment for each .sam file (per reads file) and outputs a dictionary with keys = reads and values =
## Still not implemented
# def SamDictMultipleHits(samfile):
#	input_sam = open(samfile,"r+")
#	for line in input_sam:
#		if not line.startswith("@"):
#			tab_split = line.split("\t")
#			Read = tab_split[0]
#			Flag = tab_split[1]
#			Ref = tab_split[2]
#			Leftmost_pos = tab_split[3]
#			Seq_lenght = tab_split[9]
#			if Ref != "*":
#				sam_dict[Read]= [Flag, Ref, Leftmost_pos, Seq_lenght]
#	return sam_dict

def DepthFileReader(depth_file, plasmid_length):
    depth_info = open(depth_file, "r")
    depth_dic_coverage = {}
    for line in depth_info:
        tab_split = line.split("\t")
        Reference = "_".join(tab_split[0].strip().split("_")[0:2])  ## stores only the gi for the reference
        Position = tab_split[1]
        NumReadsAlign = float(tab_split[2].rstrip("\n"))
        if Reference not in depth_dic_coverage:
            depth_dic_coverage[Reference] = {}
        depth_dic_coverage[Reference][Position] = NumReadsAlign
    Percentage_BasesCovered, Mean = {}, {}
    for ref in depth_dic_coverage:
        Percentage_BasesCovered[ref] = float(len(depth_dic_coverage[ref])) / float(plasmid_length[ref])
        Mean[ref] = sum(depth_dic_coverage[ref].values()) / len(depth_dic_coverage[ref])
    return Percentage_BasesCovered, Mean


def bar_plot(trace_list, cutoff, number_plasmid, plasmid_db_out):
    trace_line = go.Scatter(x=number_plasmid, y=[cutoff] * len(number_plasmid), mode='lines', name='cut-off',
                            marker=dict(color='rgb(255, 0, 0)'))
    trace_list.append(trace_line)
    # x bars and 1 line plot
    layout = go.Layout(barmode='group', yaxis=dict(title='Percentage of plasmid length covered'))
    fig = go.Figure(data=trace_list, layout=layout)
    plot_url = plotly.offline.plot(fig, filename=plasmid_db_out + '.html', auto_open=False)


############# PLASMIDS ##################
def PlasmidProcessing(dblist, plasmids_path, plasmid_length, output_name):
    count_entries = 0
    print
    "========================================================================="
    print
    "Processing Plasmids in " + plasmids_path
    pct = 0
    for dirname, dirnames, filenames in os.walk(plasmids_path, topdown=True):
        dirnames[:] = [d for d in dirnames if d not in ["bowtie2idx", "fasta"]]
        for filename in filenames:
            # if it is a genebank file
            if filename.endswith('.gb'):
                pct += 1
                print
                "Plasmid (.gb) file found: " + filename
                print
                "#:" + str(pct)
                gbfile = os.path.join(dirname, filename)
                FolderExist(os.path.join(dirname + "fasta"))
                fasta_file = os.path.join(dirname, 'fasta', filename[:-len('.gb')]) + '.fasta'
                if not (os.path.exists(fasta_file)):
                    print
                    "Fasta file not file. Converting..."
                    # If fasta not found Transform .gb to .fasta
                    print
                    "Creating fasta: " + fasta_file + "..."
                    plasmid_length = ExtractFastaPlasmids(gbfile, fasta_file, plasmid_length)
                else:
                    print
                    "Fasta Found! No conversion needed"
                    # if there was a previous gb->fasta conversion
                    fasta_path = os.path.join(dirname, 'fasta', os.path.splitext(filename)[0]) + '.fasta'
                    plasmid_length, fasta_entries = SequenceLengthFromFasta(fasta_file, plasmid_length, fasta_path)
                    fasta_file = fasta_path + ".temp"
                    dblist.append(fasta_file)
            elif any(x in filename for x in [".fas", ".fasta", ".fna", ".fsa", ".fa"]):
                # if it is a fasta file extension
                FolderExist(os.path.join(dirname + "fasta"))
                pct += 1
                print
                "Plasmid file (.fasta) found: " + filename
                print
                "#:" + str(pct)
                fasta_file = os.path.join(dirname, os.path.splitext(filename)[0]) + os.path.splitext(filename)[1]
                fasta_path = os.path.join(dirname, 'fasta/', os.path.splitext(filename)[0]) + '.fasta'
                print
                "Fasta Found! No conversion needed (.fasta)"
                # copyfile(fasta_file, fasta_path)
                # if there was a previous .gb->.fasta conversion
                plasmid_length, fasta_entries = SequenceLengthFromFasta(fasta_file, plasmid_length, fasta_path)
                count_entries += fasta_entries
                fasta_file = fasta_path + ".temp"
                dblist.append(fasta_file)
            else:
                print
                "File extension " + os.path.splitext(filename)[
                    1] + " not recognized! Are you sure this is a fasta or gb file? Extensions autorized are .gb, .fas, .fasta, .fna, .fsa or .fa"

    print
    "========================================================================="
    print

    ## CONCATENATES ALL PLASMID FASTA INTO A SINGLE DB ######################################
    ## This avoids the mapping of reads to a "worst matching" plasmid########################
    ## Instead reads will now be matched against the best scoring plasmid in the entire db###
    main_db = FastaConcatenation(dblist, output_name, plasmids_path)
    return main_db, count_entries


def mapper(pair, idx_file, reads_file, threads, max_k, sam_file, maindb_path):
    if pair == True and "_pair" in reads_file[0]:
        btc = 'bowtie2 -x ' + idx_file + ' -1 ' + reads_file[0] + ' -2 ' + reads_file[
            1] + ' -p ' + threads + ' -k ' + max_k + ' -5 15 -S ' + sam_file
    else:
        btc = 'bowtie2 -x ' + idx_file + ' -U ' + reads_file + ' -p ' + threads + ' -k ' + max_k + ' -5 15 -S ' + sam_file
    print
    "1) " + btc
    proc1 = Popen(btc, stdout=PIPE, stderr=PIPE, shell=True)
    out, err = proc1.communicate()
    print
    err
    regex_match = re.search('[\d]{1}[.]{1}[\d]{2}% overall alignment rate', err)
    alignment_rate = regex_match.group(0).split('%')[0]
    if alignment_rate > 0.00:
        print
        "2) " + 'samtools faidx ' + maindb_path
        call('samtools faidx ' + maindb_path, shell=True)
        bam_file = sam_file[:-3] + 'bam'
        print(
            "3) " + 'samtools view -b -S -t ' + maindb_path + '.fai' + ' -@ ' + threads + ' -o ' + bam_file + ' ' + sam_file)
        call('samtools view -b -S -t ' + maindb_path + '.fai' + ' -@ ' + threads + ' -o ' + bam_file + ' ' + sam_file,
             shell=True)
        sorted_bam_file = bam_file[:-3] + 'sorted.bam'
        print("4) " + 'samtools sort' + ' -@ ' + threads + ' -o ' + sorted_bam_file + ' ' + bam_file)
        call('samtools sort' + ' -@ ' + threads + ' -o ' + sorted_bam_file + ' ' + bam_file, shell=True)
        print("5)" + 'samtools index ' + sorted_bam_file)
        call('samtools index ' + sorted_bam_file, shell=True)
        print("6) " + 'samtools depth ' + sorted_bam_file)
        depth_file = sorted_bam_file + '_depth.txt'
        print("Creating coverage Depth File: " + depth_file)
        proc2 = Popen('samtools depth ' + sorted_bam_file + ' > ' + depth_file, stdout=PIPE, stderr=PIPE, shell=True)
        out2, err2 = proc2.communicate()
        return depth_file


def main():
    parser = argparse.ArgumentParser(
        description="Outputs a coverage percentage for each Plasmid gbk in PlasmidDir using the reads presented in the directory structure in ReadsDir")
    parser.add_argument('-p', '--plasmid', dest='plasmid_dir', required=True,
                        help='Provide the path to the directory containing plasmid fastas')
    parser.add_argument('-r', '--read', dest='read_dir', required=True,
                        help='Provide the path to the directory containing reads fastas')
    parser.add_argument('-t', '--threads', dest='threads', default="1",
                        help="Specify the number of threads to be used by bowtie2")
    parser.add_argument('-2', '--pair', dest='paired', action='store_true',
                        help='Use this option if you have paired end reads. Note that you should append to the end of the fastq file "_pair1" and "_pair2" for the script to be able to recognize both pairs')
    parser.add_argument('-k', "--max_align", dest="max_align",
                        help="Specify the maximum number of alignments possible for each read. This options changes -k parameter of Bowtie2. By default this script will set -k to the number of fastas in reference directory (e.g. if you have 3 reference sequences the number of max_align allowed will automatically be set to 3.")
    parser.add_argument('-o', '--output', dest='output_name', default="plasmid_db_out",
                        help='Specify the output name you wish. No need for file extension!')
    parser.add_argument('-c', '--cutoff', dest='cutoff_number', default="0.00",
                        help='Specify the cutoff for percentage of plasmid coverage that reads must have to be in the output. This should be a number between 0.00-1.00')

    #	parser.add_argument('--only-plasmid', dest='only_plasmid', action='store_true', help='If you just want to have the result of plasmid sequences rather than whole reads (chromosomal + plasmid reads).') ## Still not implemented
    args = parser.parse_args()

    ## Lists and dictionaries

    plasmid_length = {}
    strain_list = []
    pidx2name = {}
    dblist = []
    sam_dict = {}

    ## Process plasmids references into a single fasta

    maindb, count_entries = PlasmidProcessing(dblist, args.plasmid_dir, plasmid_length, args.output_name)
    print
    count_entries
    maindb_path = os.path.join(args.plasmid_dir + "fasta/" + maindb)

    ## Deletes temporary fastas created during PlasmidProcessing function
    deltemp(os.path.join(args.plasmid_dir + "fasta/"))

    ##Create Bowtie Idx files for plasmid references
    idx_file = CreateBowtieIdx(maindb, args.plasmid_dir, args.threads)

    ### READS#########################
    output_txt = open(args.output_name + ".txt", "w")

    master_keys = []
    trace_list = []
    counter = 0  # counter used to control output file
    file_reset = False
    for dirname, dirnames, filenames in os.walk(args.read_dir):
        for subdirname in dirnames:
            for dirname2, dirnames2, filenames2 in os.walk(os.path.join(dirname, subdirname)):
                for filename in filenames2:
                    if filename.find('fastq') != -1:
                        fn = filename.split('.')[0]
                        strain_list.append(fn)
                        print
                        print
                        "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                        print
                        print
                        "Filename: " + filename
                        print
                        datetime.fromtimestamp(time()).strftime('%Y-%m-%d %H:%M:%S')
                        # print "Mapping "+ filename+" vs "+ maindb_path
                        sam_file = dirname2 + '/' + args.output_name + '_' + subdirname + '.sam'
                        reads_file = os.path.join(dirname2, filename)
                        threads = args.threads
                        max_k = alignmaxnumber(args.max_align, count_entries)
                        if "_pair" in reads_file and args.paired and file_reset != True:
                            reads_list = [reads_file]
                            file_reset = True
                        elif "_pair" in reads_file and args.paired and file_reset == True:
                            reads_list.append(reads_file)
                            depth_file = mapper(args.paired, idx_file, reads_list, threads, max_k, sam_file,
                                                maindb_path)
                            reads_list = []
                            file_reset = False
                        else:
                            depth_file = mapper(args.paired, idx_file, reads_file, threads, max_k, sam_file,
                                                maindb_path)

                ## Compute descritptive statistics and prints to tabular txt file
                try:
                    Percentage_BasesCovered, Mean = DepthFileReader(depth_file, plasmid_length)
                    sorted_percCoverage_dic = sorted(Percentage_BasesCovered.items(), key=operator.itemgetter(1),
                                                     reverse=True)

                    tmp_list_k = []
                    tmp_list_v = []
                    list_all_k = []
                    list_all_v = []
                    if 0.00 <= float(args.cutoff_number) <= 1.00:
                        ## outputs a json file per input file
                        output_json = open(args.output_name + fn + ".json", "w")  ## new output json
                        json_dict = {}

                        if counter == 0:
                            output_txt.write(
                                "NOTE: outputed results for plasmids with more than " + args.cutoff_number + " mapping coverage.\n\n")
                            counter = 1
                        for k, v in sorted_percCoverage_dic:
                            if v >= float(args.cutoff_number):
                                tmp_list_k.append(k)
                                tmp_list_v.append(v)
                                json_dict[k] = v
                            if k not in master_keys:
                                master_keys.append(k)
                            list_all_v.append(v)
                            list_all_k.append(k)
                        ## COVERAGE PERCENTAGE ##
                        output_txt.write(
                            "Read name: " + fn + "\nReference Sequence\tCoverage Percentage\tMean mapping depth\n")
                        # count_x=0
                        for element in tmp_list_k:
                            output_txt.write(str(element) + "\t")  # outputs ref sequence
                            output_txt.write(
                                str(Percentage_BasesCovered[element]) + "\t")  # outputs coverage percentage
                            output_txt.write(str(Mean[element]) + "\n")  # outputs mean mapping depth
                        # count_x += 1

                        trace = go.Bar(x=list_all_k, y=list_all_v, name=fn)
                        trace_list.append(trace)
                    output_txt.write("\n")
                    print(json_dict)
                    output_json.write(json.dumps(json_dict))
                    output_json.close()
                except NameError:
                    depth_file = None

    ### Graphical outputs ###
    bar_plot(trace_list, float(args.cutoff_number), master_keys, args.output_name)
    output_txt.close()


if __name__ == "__main__":
    main()
