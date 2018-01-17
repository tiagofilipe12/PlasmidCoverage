#!/usr/bin/env python3

# Last update: 25/12/2017
# Author: T.F. Jesus
# This version runs with bowtie2 build 2.2.9, with multithreading for
# bowtie2-build

import argparse
import os
import re
import operator
import shutil
import plotly
import plotly.graph_objs as go
from subprocess import Popen, PIPE
from Bio import SeqIO
from time import time
from datetime import datetime
import json
from termcolor import cprint

class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def search_substing(string):
    '''
    A function to search for a substring containing the plasmid name
    Parameters
    ----------
    string: str
        The string to be searched. Usually a sequence header, a fasta entry.

    Returns
    -------
    plasmid_name: str
        The name of the plasmid if the regex returns anything related with
        "plasmid"

    '''
    plasmid_search = re.search('plasmid(.+?)__', string)
    if plasmid_search:
        plasmid_name = plasmid_search.group(1).replace("_", "")
        return plasmid_name

def alignmaxnumber(max_align, count_entries):
    '''
    A function that ascertains if max_align option was provided. If no
    max_align is provided by the user (with the -k option) then max_align
    will be set to the number of entries in fastas, i.e., the total number of
    sequences present in all fastas provided as input
    Parameters
    ----------
    max_align: str
        The args passed with -k, if any
    count_entries: float
        The number of entries in fasta

    Returns
    -------
    k_value: str
        The value to be set for -k option for bowtie2

    '''
    if max_align:
        k_value = max_align
    else:
        k_value = str(count_entries)
    return k_value


def folderexist(directory):
    '''
    A function that checks if a directory already exists
    Parameters
    ----------
    directory: str
        The path to the directory to be checked

    '''
    if not directory.endswith("/"):
        directory = directory + "/"
    if not os.path.exists(os.path.join(directory)):
        os.makedirs(os.path.join(directory))
        print(os.path.join(directory) + " does not exist. One will be "
                                        "created...")
    else:
        print(os.path.join(directory) + " exists!")

def fastadict(fasta_file):
    '''
    Creates an associative dictionary between the sequences name / header and
    the sequence itself
    Parameters
    ----------
    fasta_file: str
        The path to the fasta file to be parsed here

    Returns
    -------
    fasta_dic: dict
        A dictionary with keys as the sequence headers and the sequence as
        the corresponding value

    '''
    if_handle = open(fasta_file, 'r')
    x = 0
    fasta_dic = {}
    problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/",
                              "[", "]", ":", "{", "}"]
    sequence_list = []
    for line in if_handle:
        if len(line) > 0:
            line = line.splitlines()[0]
        if x == 0 and line.startswith(">"):
            sequence_list = []
            plasmidname = line.replace(">", "")  # stores only the
            #  acc for the reference
            for char in plasmidname:
                if char in problematic_characters:
                    plasmidname = plasmidname.replace(char, "_")
            x += 1
        elif x == 0 and not line.startswith(">"):
            print("Is this a fasta file? " + fasta_file)
            print(fasta_file + " will be ignored")
            break
        elif x >= 1 and line.startswith(">"):
            fasta_dic[plasmidname] = sequence_list  #
            # appends
            # last
            # sequence
            # to be
            # parsed before new structure for sequence
            sequence_list = []
            plasmidname = line.replace(">", "")  # stores only the
            for char in plasmidname:
                if char in problematic_characters:
                    plasmidname = plasmidname.replace(char, "_")  # stores only
                    #  the gi for the reference
            x += 1
        else:
            sequence_list.append(line)
    if sequence_list:
        fasta_dic[plasmidname] = sequence_list  # appends last sequence on the
        # fasta
    if_handle.close()
    return fasta_dic

def sequencelengthfromfasta(fasta_file, plasmid_length, fasta_path):
    '''
    A function to obtain sequence length when parsing fastas. Actually counts
    the number of bases in each accession
    Parameters
    ----------
    fasta_file: str
        The path to the fasta file to be parsed here
    plasmid_length: dict
        The dictionary that will store the key: value for accessions: lengths (empty)
    fasta_path: str
        The path to the fasta files

    Returns
    -------
    plasmid_length: dict
        The dictionary of accessions and lengths (populated)
    len(fasta_dic.keys())

    '''
    fasta_dic = fastadict(fasta_file)
    out_handle = open(os.path.join(fasta_path + ".temp"), "w")
    for key in fasta_dic:
        new_key = "_".join(key.split("_")[0:3])  # stores only the acc for the
        # reference
        plasmid_length[new_key] = sum(len(s) for s in fasta_dic[key])
        out_handle.write(">" + key + "\n" + "".join(fasta_dic[key]) + "\n")
    out_handle.close()
    return plasmid_length, len(fasta_dic.keys())

def extractfastaplasmids(gbkfile, fastafile, plasmid_length):
    '''
    A legacy function that was used to extract fastas from gb files
    Parameters
    ----------
    gbkfile: str
        The path to the genbank file
    fastafile: str
        A path to the fasta file
    plasmid_length: dict
        Stores the length of each sequence in fastas (empty)

    Returns
    -------
    plasmid_length: dict
        Stores the length of each sequence in fastas (populated)

    '''
    if_handle = open(gbkfile, 'r')
    gbkdata = SeqIO.read(if_handle, "genbank")
    problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/",
                              "[", "]", ":", "{", "}"]
    plasmid_name = gbkdata.description[:-1]
    for char in plasmid_name:
        if char in problematic_characters:
            plasmid_name = plasmid_name.replace(char, "_")
    sequence = str(gbkdata.seq)
    if not sequence.endswith("\n"):
        sequence = sequence + "\n"
    out_handle = open(fastafile, 'w')
    out_handle.write('>' + plasmid_name + '\n')
    out_handle.write(sequence)
    if_handle.close()
    out_handle.close()
    plasmid_length[plasmid_name] = len(str(gbkdata.seq))
    print(plasmid_length.values(), " values each gb file")
    print(" Wrote fasta file: " + fastafile)
    return plasmid_length

def createbowtieidx(filename, dirname, threads):
    '''
    A function to create the bowtie2 index in the case -idx option is not
    selected and -p option is provided
    Parameters
    ----------
    filename: str
        The path to the concatenated fasta file with all sequences
    dirname: str
        The path to the directory that contains the fasta files used to make
        the concatenated fasta and that here will be used to generate the
        bowtie2 index
    threads: str
        The number of threads to be used by bowtie2 to create the index files

    Returns
    -------
    idx_file: str
        A string with the common name for all 6 bowtie2 index files (given
        that bowtie2 automatically can detect all these file given the
        correct file name pattern)

    '''
    if not os.path.exists(os.path.join(dirname + "bowtie2idx")):
        os.makedirs(os.path.join(dirname + "bowtie2idx"))
    else:
        pass
    idx_file = os.path.join(dirname, "bowtie2idx",
                            os.path.splitext(filename)[0]) + ".idx"
    fasta_file = os.path.join(dirname, "fasta",
                              os.path.splitext(filename)[0]) + ".fasta"
    list_idxfiles = [idx_file + ".1.bt2", idx_file + ".2.bt2",
                     idx_file + ".3.bt2", idx_file + ".4.bt2",
                     idx_file + ".rev.2.bt2", idx_file + ".rev.1.bt2"]
    for idx in list_idxfiles:
        if not os.path.isfile(idx):
            print("Creating " + idx_file)
            # Create bowtie index
            bowtieidx_cmd = [
                "bowtie2-build",
                "-q",
                fasta_file,
                "--threads",
                threads,
                idx_file
            ]
            p = Popen(bowtieidx_cmd, stdout = PIPE, stderr = PIPE)
            p.wait()
            #call('bowtie2-build -q ' + fasta_file + ' --threads ' + threads +
            #     ' ' + idx_file, shell=True)  # convert to popen
        else:
            print(idx_file + " already exists!")

    return idx_file

def fastaconcatenation(dblist, output_name, plasmid_dir):
    '''
    A function to concatenate fasta files into a main fasta file
    Parameters
    ----------
    dblist: list
        A list that stores the paths to the files containing the plasmid
        sequences
    output_name: str
        The name of the outputs specified by the user
    plasmid_dir: str
        The path to fastas containing fastas

    Returns
    -------
    main_filename: str
        The name of the main concatenated fasta file

    '''
    #print(dblist)
    main_filename = output_name + ".fasta"
    #print(main_filename)
    dirname = os.path.join(plasmid_dir, "fasta", "")
    #print(dirname)
    if os.path.isfile(dirname + main_filename):
        print(output_name + ".fasta already exists. Overriding file...")
    print("Saving to: " + main_filename)

    python_cat(dblist, dirname + main_filename)

    #p = Popen("cat " + ' '.join(dblist) + " > " + dirname + main_filename,
      #        stdout=PIPE, stderr=PIPE, shell=True)

    #stdout, stderr = p.communicate()
    #print(stdout, stderr)
    #p.wait()

    return main_filename

def python_cat(dblist, output_name):
    '''
    A function to concatenate files
    Parameters
    ----------
    dblist: list
        A list that stores the paths to the files
        containing the plasmid sequences
    output_name: str
        The name of the outputs specified by the user

    '''
    destination = open(output_name, "wb") #write and binary
    for db in dblist:
        shutil.copyfileobj(open(db, "rb"), destination)
    destination.close()
    #doesn"t need to return... just create the file

# function to delete temporary fasta files
def deltemp(directory):
    '''
    A function to delete temporary directories
    Parameters
    ----------
    directory: str
        The directory to be removed

    '''
    files = os.listdir(directory)
    print("Deleting temporary fasta files in: " + directory)
    for f in files:
        if f.endswith(".temp"):
            os.remove(os.path.join(directory, f))

def depthfilereader(depth_file, plasmid_length):
    '''
    Function that parse samtools depth file and creates 3 dictionaries that
    will be useful to make the outputs of this script, both the tabular file
    and the json file that may be imported by pATLAS

    Parameters
    ----------
    depth_file: str
        the path to depth file for each sample
    plasmid_length: dict
        a dictionary that stores length of all plasmids in fasta given as input

    Returns
    -------
    percentage_basescovered: dict
            stores the percentage of the total sequence of a
            reference/accession (plasmid) in a dictionary
    mean: dict
        Stores the mean coverage for each reference/accession (plasmid) in a
        dictionary
    metadata: dict
        Stores information of species, plasmid_name and sequence length for
        each plasmid key in a List
    '''
    metadata = {}
    depth_info = open(depth_file, "r")
    depth_dic_coverage = {}
    for line in depth_info:
        tab_split = line.split("\t")
        reference = "_".join(tab_split[0].strip().split("_")[0:3])  # store
        species = "_".join(tab_split[0].strip().split("_")[3:5])
        plasmid_name = search_substing(line)
        # only the gi for the reference
        position = tab_split[1]
        numreadsalign = float(tab_split[2].rstrip("\n"))
        if reference not in depth_dic_coverage:
            depth_dic_coverage[reference] = {}
        depth_dic_coverage[reference][position] = numreadsalign
        metadata[reference] = [species, plasmid_name, plasmid_length[reference]]
    #print(depth_dic_coverage)
    percentage_basescovered, mean = {}, {}
    for ref in depth_dic_coverage:
        percentage_basescovered[ref] = float(len(depth_dic_coverage[ref])) / \
                                       float(plasmid_length[ref])
        mean[ref] = sum(depth_dic_coverage[ref].values()) / \
                    len(depth_dic_coverage[ref])
    return percentage_basescovered, mean, metadata

def bar_plot(trace_list, cutoff, number_plasmid, plasmid_db_out):
    '''
    A function to create bar plots for the results. This is now a legacy
    function and that may become deprecate in further implementations of this script
    Parameters
    ----------
    trace_list
    cutoff
    number_plasmid
    plasmid_db_out

    '''
    trace_line = go.Scatter(x=number_plasmid, y=[cutoff] * len(number_plasmid),
                            mode="lines", name="cut-off",
                            marker=dict(color="rgb(255, 0, 0)"))
    trace_list.append(trace_line)
    # x bars and 1 line plot
    layout = go.Layout(barmode="group", yaxis=dict(title="Percentage of "
                                                         "plasmid length "
                                                         "covered"))
    fig = go.Figure(data=trace_list, layout=layout)
    plotly.offline.plot(fig, filename=plasmid_db_out + ".html",
                                   auto_open=False)

# PLASMIDS #
def plasmidprocessing(dblist, plasmids_path, plasmid_length, output_name):
    '''
    Function to parse multiple fasta files and generate a main fasta to be
    used by bowtie2
    Parameters
    ----------
    dblist: list
        An empty list that will be used to store the paths to the files
        containing the plasmid sequences
    plasmids_path: str
        The path to fastas containing fastas
    plasmid_length: dict
        A dictionary that will store a reference between each accession
        number and its sequence length
    output_name: str
        The string naming the output

    Returns
    -------
    main_db: str
        The path to the concatenated fasta file with all sequences
    count_entries: float
        The number of sequences present in main_db

    '''
    count_entries = 0
    print("===================================================================")
    cprint("Processing Plasmids in " + plasmids_path, "green", attrs=["bold"])
    pct = 0
    for dirname, dirnames, filenames in os.walk(plasmids_path, topdown=True):
        dirnames[:] = [d for d in dirnames if d not in ["bowtie2idx", "fasta"]]
        for filename in filenames:
            # if it is a genebank file
            if filename.endswith(".gb"):
                pct += 1
                print("Plasmid (.gb) file found: " + filename)
                print("\n")
                print("#:" + str(pct))
                gbfile = os.path.join(dirname, filename)
                folderexist(os.path.join(dirname + "fasta"))
                fasta_file = os.path.join(dirname, "fasta",
                                          filename[:-len(".gb")]) + ".fasta"
                if not (os.path.exists(fasta_file)):
                    print("Fasta file not file. Converting...")
                    # If fasta not found Transform .gb to .fasta
                    print("Creating fasta: " + fasta_file + "...")
                    plasmid_length = extractfastaplasmids(gbfile, fasta_file,
                                                          plasmid_length)
                else:
                    print("Fasta Found! No conversion needed")
                    # if there was a previous gb->fasta conversion
                    fasta_path = os.path.join(dirname, "fasta",
                                              os.path.splitext(filename)[0]) + \
                                 ".fasta"
                    plasmid_length, fasta_entries = sequencelengthfromfasta(
                        fasta_file, plasmid_length, fasta_path)
                    fasta_file = fasta_path + ".temp"
                    dblist.append(fasta_file)
            elif any(x in filename for x in [".fas", ".fasta", ".fna",
                                             ".fsa", ".fa"]):
                # if it is a fasta file extension
                folderexist(os.path.join(dirname + "fasta"))
                pct += 1
                print("Plasmid file (.fasta) found: " + filename)
                print("#:" + str(pct))
                fasta_file = os.path.join(dirname,
                                          os.path.splitext(filename)[0]) + \
                             os.path.splitext(filename)[1]
                fasta_path = os.path.join(dirname, "fasta/",
                                          os.path.splitext(filename)[0]) + \
                             ".fasta"
                print("Fasta Found! No conversion needed (.fasta)")
                # copyfile(fasta_file, fasta_path)
                # if there was a previous .gb->.fasta conversion
                plasmid_length, fasta_entries = sequencelengthfromfasta(
                    fasta_file, plasmid_length, fasta_path)
                count_entries += fasta_entries
                fasta_file = fasta_path + ".temp"
                dblist.append(fasta_file)
            else:
                print("File extension {0} not recognized! Are you sure this is"
                      " a fasta or gb file? Extensions autorized are .gb, .fas,"
                      " .fasta, .fna, .fsa or .fa".format(
                    os.path.splitext(filename)[
                        1]))

    print("==============================================================="
          "======")
    print("\n")

    # CONCATENATES ALL PLASMID FASTA INTO A SINGLE DB #
    # This avoids the mapping of reads to a "worst matching" plasmid#
    # Instead reads will now be matched against the best scoring plasmid in
    # the entire db#
    main_db = fastaconcatenation(dblist, output_name, plasmids_path)
    # crates a length_json file that is used when index file is used instead
    #  of fasta files.
    length_json = open(output_name + "_length.json", "w")
    length_json.write(json.dumps(plasmid_length))

    return main_db, count_entries

def mapper(pair, idx_file, reads_file, threads, max_k, sam_file, maindb_path,
           trim5, indexes):
    '''
    The function that actually runs bowtie2 and samtools after creating the
    bowtie2 index file and processing the input fastas
    Parameters
    ----------
    pair: bool
        A boolean which tells if paired end is provided in reads folders. If
        not defined then single read file mode is provided to bowtie2
    idx_file: str
        A string with the file path to bowtie2 index files
    reads_file: str
        A string with the path to the location of reads
    threads: str
        The number of threads to be used by bowtie2 and samtools
    max_k: str
        the number of maximum possible matches for a sequence (-k parameter
        for bowtie2
    sam_file: str
        The path to sam_file
    maindb_path: str
        The path to the main fasta file (concatenated plasmid fasta file).
        This is used to either create .fai index for samtools or to indicate
        the path to it
    trim5: str
        indicates the number of bases that should be trimmed from the left of
        the sequence. It is a parameter for bowtie2 (-5).
    indexes: bool
        indicates whether indexes file are provided with the -idx option or not.
        This is useful for the script to know if it needs to create a new
        samtools .fai index or if it can use an existing one.

    Returns
    -------
    depth_file: str
        This is a string that stores the path to depth file that will be used
        afterwards for .txt and .json outputs for that read set

    '''
    cprint("\n=== Running bowtie2 ===\n", "green", attrs=["bold"])
    if pair == True:
        btc = ["bowtie2", "-x", idx_file, "-1", reads_file[0], "-2",
              reads_file[1], "-p", threads, "-k", max_k, "-5", trim5, "-S",
              sam_file]
    else:
        btc = ["bowtie2", "-x", idx_file, "-U", reads_file, "-p",
              threads, "-k", max_k, "-5", trim5, "-S", sam_file]
    print("1) " + " ".join(btc))
    proc1 = Popen(btc, stdout = PIPE, stderr = PIPE)
    proc1.wait()
    #proc1 = Popen(btc, stdout=PIPE, stderr=PIPE, shell=True)
    out, err = proc1.communicate()
    #regex_match = re.search("[\d]{1}[.]{1}[\d]{2}% overall alignment rate",
    # err)
    #try:
    #    alignment_rate = regex_match.group(0).split("%")[0]
    #except AttributeError:
    #    print(err)
    #    print("\nWARNING: bowtie2-build found matching file types and escaped "
    #          "building "
    #          "new index, however the specified file name does not match "
    #          "bowtie index. Try renaming the output '-o' option to match "
    #          "that of the bowtie2 idx files.\n")
    #if alignment_rate > 0:
    cprint("\n=== Running samtools ===\n", "green", attrs=["bold"])
    if indexes == False:
        print("2) " + "samtools faidx " + maindb_path)
        proc2 = Popen(["samtools", "faidx", maindb_path],
                      stdout = PIPE,
                      stderr = PIPE)
        proc2.wait()
    else:
        print("2) " + maindb_path + ".fai found! Proceeding to next step")
    #call('samtools faidx ' + maindb_path, shell=True)
    bam_file = sam_file[:-3] + "bam"
    print("3) " + "samtools view -b -S -t " + maindb_path + ".fai" +
          " -@ " + threads + " -o " + bam_file + " " + sam_file)
    samtools_view_cmd = [
        "samtools",
        "view",
        "-b",
        "-S",
        "-t",
        maindb_path + ".fai",
        "-@",
        threads,
        "-o",
        bam_file,
        sam_file
    ]
    proc3 = Popen(samtools_view_cmd, stdout = PIPE, stderr = PIPE)
    proc3.wait()
    #call('samtools view -b -S -t ' + maindb_path + '.fai' +
    #     ' -@ ' + threads + ' -o ' + bam_file + ' ' + sam_file,
    #     shell=True)
    sorted_bam_file = bam_file[:-3] + "sorted.bam"
    print("4) " + "samtools sort" + " -@ " + threads + " -o " +
          sorted_bam_file + " " + bam_file)
    samtools_sort_cmd = [
        "samtools",
        "sort",
        "-@",
        threads,
        "-o",
        sorted_bam_file,
        bam_file
    ]
    proc4 = Popen(samtools_sort_cmd, stdout = PIPE, stderr = PIPE)
    proc4.wait()
    #call("samtools sort" + " -@ " + threads + " -o " +
    #     sorted_bam_file + " " + bam_file, shell=True)
    print("5) " + "samtools index " + sorted_bam_file)
    samtools_index_cmd = [
        "samtools",
        "index",
        sorted_bam_file
    ]
    proc5 = Popen(samtools_index_cmd, stdout = PIPE, stderr = PIPE)
    proc5.wait()
    #call("samtools index " + sorted_bam_file, shell=True)
    print("6) " + "samtools depth " + sorted_bam_file)
    depth_file = sorted_bam_file + "_depth.txt"
    print("Creating coverage Depth File: " + depth_file)
    #samtools_depth_cmd = [
    #    "samtools",
    #    "depth",
    #    sorted_bam_file,
    #    ">",
    #    depth_file
    #]
    samtools_depth_cmd = "samtools depth {} > {}".format(sorted_bam_file,
                                                         depth_file)
    proc6 = Popen(samtools_depth_cmd, stdout=PIPE, stderr=PIPE, shell=True)
    proc6.wait()
    return depth_file

def main():
    parser = argparse.ArgumentParser(
        description="Outputs a coverage percentage for each Plasmid gbk in "
                    "PlasmidDir using the reads presented in the directory "
                    "structure in ReadsDir")
    parser.add_argument("-p", "--plasmid", dest="plasmid_dir", action=FullPaths,
                        help="Provide the path to the directory containing "
                             "plasmid fastas")
    parser.add_argument("-idx", "--indexes_folder", dest="indexes", action=FullPaths,
                        help="Provide the path to indexes folders")
    parser.add_argument("-r", "--read", dest="read_dir", required=True,
                        action=FullPaths,
                        help="Provide the path to the directory containing "
                             "reads fastq files")
    parser.add_argument("-t", "--threads", dest="threads", default="1",
                        help="Specify the number of threads to be used by "
                             "bowtie2")
    parser.add_argument("-2", "--pair", dest="paired", action="store_true",
                        help="Use this option if you have paired end reads. "
                             "Paired end reads just have to be inside the "
                             "same folder, without any other files inside it. "
                             "Structure to read files should be something "
                             "like: ./reads/read_folder/<with pairs inside "
                             "it>.")
    parser.add_argument("-k", "--max_align", dest="max_align",
                        help="Specify the maximum number of alignments possible "
                             "for each read. This options changes -k parameter "
                             "of Bowtie2. By default this script will set -k to "
                             "the number of fastas in reference directory (e.g. "
                             "if you have 3 reference sequences the number of "
                             "max_align allowed will automatically be set to 3.")
    parser.add_argument("-5", "--trim5", dest="trim5", default="0",
                        help="bowtie2 option: Trim <int> bases from 5' (left)"
                             "end of each read before alignment (default: 0).")
    parser.add_argument("-o", "--output", dest="output_name",
                        default="plasmid_db_out",
                        help="Specify the output name you wish. No need for "
                             "file extension!")
    parser.add_argument("-c", "--cutoff", dest="cutoff_number", default="0.00",
                        help="Specify the cutoff for percentage of plasmid "
                             "coverage that reads must have to be in the "
                             "output. This should be a number between 0.00-1.00")

    # parser.add_argument('--only-plasmid', dest='only_plasmid',
    # action='store_true', help='If you just want to have the result of plasmid sequences rather than whole reads (chromosomal + plasmid reads).') # Still not implemented
    args = parser.parse_args()

    # Lists and dictionaries
    strain_list = []
    dblist = []
    # args
    plasmids_dir = args.plasmid_dir
    reads_dir = args.read_dir
    index_path = args.indexes

    # check the format of input directories to -r and -p options
    if args.plasmid_dir:
        if not plasmids_dir.endswith("/"):
            plasmids_dir += "/"
    if not reads_dir.endswith("/"):
        reads_dir += "/"
    if not index_path.endswith("/"):
        index_path += "/"

    # Process plasmids references into a single fasta

    if not index_path:
        plasmid_length = {}
        maindb, count_entries = plasmidprocessing(dblist, plasmids_dir,
                                                  plasmid_length, args.output_name)

        maindb_path = os.path.join(plasmids_dir + "fasta/" + maindb)

        # Deletes temporary fastas created during plasmidprocessing function
        deltemp(os.path.join(plasmids_dir + "fasta/"))

    # Create Bowtie Idx files for plasmid references
        idx_file = createbowtieidx(maindb, plasmids_dir, args.threads)
        indexes = False
    else:
        idx_file = os.path.join(index_path + "bowtie2idx/" +
                                "bowtie2.idx")
        #  fasta
        maindb_path = os.path.join(index_path + "fasta/" + "samtools.fasta")
        indexes = True
        # import json_length from json file
        script_dir = os.path.dirname(__file__)
        rel_path = "/json/reads_sample_result_length.json"
        length_import = open(script_dir + rel_path)
        plasmid_length = json.load(length_import)
        count_entries = len(plasmid_length)

    # READS#
    output_txt = open(args.output_name + ".txt", "w")

    master_keys = []
    trace_list = []
    counter = 0  # counter used to control output file
    file_reset = False
    for dirname, dirnames, filenames in os.walk(reads_dir):
        for subdirname in dirnames:

            for dirname2, dirnames2, filenames2 in os.walk(os.path.join(dirname,
                                                                        subdirname)):
                for filename in filenames2:
                    if filename.find('fastq') != -1 or filename.find('fq') !=\
                            -1:
                        fn = filename.split('.')[0]
                        strain_list.append(fn)
                        print("\n")
                        print("+++++++++++++++++++++++++++++++++++++++++++++++"
                              "++++++++++++++++++++++++++++++++++++++++")
                        print("\n")
                        print("Filename: " + filename)
                        print("\n")
                        datetime.fromtimestamp(time()).strftime(
                            '%Y-%m-%d %H:%M:%S')
                        # print "Mapping "+ filename+" vs "+ maindb_path
                        sam_file = dirname2 + '/' + args.output_name + '_' + \
                                   subdirname + '.sam'
                        reads_file = os.path.join(dirname2, filename)
                        threads = args.threads
                        max_k = alignmaxnumber(args.max_align, count_entries)
                        if args.paired \
                                and file_reset != True:
                            reads_list = [reads_file]
                            file_reset = True
                        elif args.paired \
                                and file_reset == True:
                            reads_list.append(reads_file)
                            depth_file = mapper(args.paired, idx_file,
                                                reads_list, threads, max_k,
                                                sam_file, maindb_path,
                                                args.trim5, indexes)
                            reads_list = []
                            file_reset = False
                        else:
                            depth_file = mapper(args.paired, idx_file,
                                                reads_file, threads, max_k,
                                                sam_file, maindb_path,
                                                args.trim5, indexes)

                # Compute descriptive statistics and prints to tabular txt file

                #try:
                percentage_basescovered, mean, metadata = depthfilereader(
                    depth_file, plasmid_length)
                # print(percentage_basescovered)
                # print(mean)
                # print(metadata)
                sorted_perccoverage_dic = sorted(
                    percentage_basescovered.items(),
                    key=operator.itemgetter(1),
                    reverse=True)
                tmp_list_k, list_all_k, list_all_v = [], [], []
                if 0 <= float(args.cutoff_number) <= 1:
                    # outputs a json file per input file
                    output_json = open(args.output_name + fn + ".json",
                                       "w")  # new output json
                    json_dict = {}

                    if counter == 0:
                        output_txt.write(
                            "NOTE: outputted results for plasmids with "
                            "more than " + args.cutoff_number + " mapping"
                                                                " coverage.\n\n")
                        counter = 1
                    for k, v in sorted_perccoverage_dic:
                        if v >= float(args.cutoff_number):
                            tmp_list_k.append(k)
                            #tmp_list_v.append(v)
                            json_dict[k] = v
                        if k not in master_keys:
                            master_keys.append(k)
                        list_all_v.append(v)
                        list_all_k.append(k)

                    # COVERAGE PERCENTAGE #
                    output_txt.write(
                        "Read name: " + fn + "\nReference sequence\t"
                                             "Coverage percentage\t"
                                             "Mean mapping depth\t"
                                             "Sequence length\t"
                                             "Species name\t"
                                             "Plasmid name\n")
                    # count_x=0
                    for element in tmp_list_k:
                        output_txt.write(
                            str(element) + "\t")  # outputs ref sequence
                        output_txt.write(
                            str(percentage_basescovered[
                                    element]) + "\t")  # outputs coverage percentage
                        output_txt.write(str(mean[element]) + "\t")
                        # outputs mean mapping depth
                        output_txt.write("{}\t{}\t{}\n".format(str(metadata[
                                                                       element][2]), " ".join(str(metadata[element][
                                                                                                      0]).split("_")), str(metadata[element][1])))

                    # count_x += 1
                    trace = go.Bar(x=list_all_k, y=list_all_v, name=fn)
                    trace_list.append(trace)
                output_txt.write("\n")
                output_json.write(json.dumps(json_dict))
                output_json.close()
                #except NameError:
                #    depth_file = None
                #    print("error: samtools depth file not correct -> ",
                #          depth_file)

    # Graphical outputs #
    bar_plot(trace_list, float(args.cutoff_number), master_keys,
             args.output_name)
    output_txt.close()

if __name__ == "__main__":
    main()
