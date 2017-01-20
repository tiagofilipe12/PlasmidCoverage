#!/usr/bin/env python

## Last update: 20/1/2017
## Author: T.F. Jesus
## This version runs with bowtie2 build 2.2.9, with multithreading for bowtie2-build

import argparse
import os
import re
import operator
import plotly
import plotly.plotly as py
import plotly.graph_objs as go
from subprocess import Popen, PIPE, call
from Bio import SeqIO
from time import time
from datetime import datetime
from shutil import copyfile


#TODO: 
#1) Set a minnimun read threshold to avoid computation on less than x reads mapped (depth_file?)
#4) delete unsorted BAM files
#5) Create directory structure to store all comparisons.

def alignmaxnumber(max_align, dblist):
	if max_align:
		k_value=max_align
	else:
		k_value = str(len(dblist))
	return k_value

def FolderExist(directory):
	if not directory.endswith("/"):
		directory = directory + "/"
	if not os.path.exists(os.path.join(directory)):		
		os.makedirs(os.path.join(directory))
		print os.path.join(directory) + " does not exist. One will be created..."
	else:
		print os.path.join(directory) + " exists!"	

def fastadict(fasta_file):
	if_handle=open(fasta_file,'r')
	x = 0
	fasta_dic = {}
	problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/"]

	for line in if_handle:
		if len(line) > 0: 
			line = line.splitlines()[0]  
		if x == 0 and line.startswith(">"):
			sequence =""
			sequence_list =[]
			PlasmidName = line[1:]
			for char in problematic_characters:
					PlasmidName = PlasmidName.replace(char, '_')
			x+=1
		elif x == 0 and not line.startswith(">"):
			print "Is this a fasta file? " + fasta_file
			raise SystemExit
		elif x >=1 and line.startswith(">"):
			fasta_dic[PlasmidName] = sequence_list	#appends last sequence to be parsed before new structure for sequence
			sequence_list =[]
			sequence =""
			PlasmidName = line[1:]
			for char in problematic_characters:
					PlasmidName = PlasmidName.replace(char, '_')
			x+=1
		else:
			sequence_list.append(line)
	fasta_dic[PlasmidName] = sequence_list	#appends last sequence on the fasta
	if_handle.close()
	return fasta_dic

def SequenceLengthFromFasta(fasta_file,plasmid_length,fasta_path):
	fasta_dic = fastadict(fasta_file)
	print len(fasta_dic.keys())
	print len(fasta_dic.values())
	out_handle = open(os.path.join(fasta_path + ".temp"), 'w')
	for key in fasta_dic:
		plasmid_length[key]=sum(len(s) for s in fasta_dic[key])
		out_handle.write('>' + key + '\n' + ''.join(fasta_dic[key]) + '\n')
	out_handle.close()
	return plasmid_length 

def ExtractFastaPlasmids(gbkfile,fastafile,plasmid_length):
	if_handle=open(gbkfile,'r')
	gbkdata=SeqIO.read(if_handle, "genbank")
	problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/"]
	plasmid_name=gbkdata.description[:-1]
	for char in problematic_characters:
		plasmid_name=plasmid_name.replace(char, '_')
	sequence = str(gbkdata.seq)
	if not sequence.endswith("\n"):
		sequence = sequence + "\n"
	out_handle=open(fastafile,'w')
	out_handle.write('>'+plasmid_name+'\n')
	out_handle.write(sequence)
	if_handle.close() 
	out_handle.close()	
	plasmid_length[plasmid_name]=len(str(gbkdata.seq))
	print plasmid_length.values(), " values each gb file"
	print(" Wrote fasta file: "+ fastafile)
	return plasmid_length
	

def CreateBowtieIdx(filename):
	check =""
	dirname = args.plasmid_dir
	if not os.path.exists(os.path.join(dirname + "bowtie2idx")):		
		os.makedirs(os.path.join(dirname + "bowtie2idx"))
	else:
		pass	
	idx_file=os.path.join(dirname,'bowtie2idx',os.path.splitext(filename)[0])+'.idx'
	fasta_file = os.path.join(dirname,'fasta',os.path.splitext(filename)[0])+'.fasta'
	list_idxfiles = [idx_file + ".1.bt2", idx_file + ".2.bt2", idx_file + ".3.bt2", idx_file + ".4.bt2", idx_file + ".rev.2.bt2", idx_file + ".rev.1.bt2"]
	for idx in list_idxfiles:
		if not os.path.isfile(idx): 
			check = "yes"
	if "yes" in check:
		print("Creating " + idx_file)
			#Create bowtie index 			
		call('bowtie2-build -q '+ fasta_file +' --threads ' +args.threads+ ' ' +idx_file, shell=True)		#convert to popen
	else:
		print idx_file + " already exists!"

	return idx_file

def FastaConcatenation(dblist):
	print(dblist)
	main_filename = args.output_name + ".fasta"
	dirname = os.path.join(args.plasmid_dir, "fasta", "")
	if os.path.isfile(dirname + main_filename):
		print args.output_name + ".fasta already exists. Overriding file..."	
	print "Saving to: " + main_filename
	concat=Popen("cat " + ' '.join(dblist) + "> " + dirname + main_filename, stdout = PIPE, stderr = PIPE, shell=True)
	stdout,stderr= concat.communicate()
	print stderr
	return main_filename

## function to delete temporary fasta files
def deltemp(directory):
	files = os.listdir(directory)
	print "Deleting temporary fasta files in: " + directory
	for f in files:
		if f.endswith('.temp'):
			os.remove(os.path.join(directory,f))

## Creates an alignment for each .sam file (per reads file) and outputs a dictionary with keys = reads and values = 
## Still not implemented
#def SamDictMultipleHits(samfile):
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
		Reference = tab_split[0].strip()
		Position = tab_split[1]
		NumReadsAlign = float(tab_split[2].rstrip("\n"))
		if Reference not in depth_dic_coverage:
			depth_dic_coverage[Reference] = {} 
		depth_dic_coverage[Reference][Position] = NumReadsAlign 
	Percentage_BasesCovered, Mean = {}, {}
	for ref in depth_dic_coverage:
		Percentage_BasesCovered[ref] = float(len(depth_dic_coverage[ref]))/float(plasmid_length[ref])
		Mean[ref] = sum(depth_dic_coverage[ref].values())/len(depth_dic_coverage[ref])
	return Percentage_BasesCovered, Mean 

def bar_plot(trace_list, cutoff, number_plasmid):
	trace_line = go.Scatter(x=number_plasmid, y=[cutoff]*len(number_plasmid), mode = 'lines', name='cut-off', marker=dict(color= 'rgb(255, 0, 0)')) 
	trace_list.append(trace_line)
	data = trace_list
	layout = go.Layout(barmode='group', yaxis= dict(title='Percentage of plasmid length covered'))#, shapes=[{'type':'line', 'x0' = 0,  'y0'= cutoff, 'x1': len(number_plasmid), 'y1': cutoff,},],)
	fig = go.Figure(data=data, layout=layout)
	plot_url = plotly.offline.plot(fig, filename='graph_new.html',auto_open=False)	
	

############# PLASMIDS ##################
def PlasmidProcessing(dblist,plasmids_path,plasmid_length):
	print "========================================================================="
	print "Processing Plasmids in "+ plasmids_path
	pct=0;
	for dirname, dirnames, filenames in os.walk(plasmids_path, topdown=True): 
		dirnames[:] = [d for d in dirnames if d not in ["bowtie2idx", "fasta"]]
		for filename in filenames:
			#if it is a genebank file
			if filename.endswith('.gb'):
				pct+=1
				print "Plasmid (.gb) file found: "+ filename
				print "#:"+str(pct)
				gbfile = os.path.join(dirname, filename)
				FolderExist(os.path.join(dirname + "fasta"))			
				fasta_file = os.path.join(dirname,'fasta',filename[:-len('.gb')])+'.fasta'
				if not(os.path.exists(fasta_file)):
					print "Fasta file not file. Converting..."
					#If fasta not found Transform .gb to .fasta
					print "Creating fasta: " + fasta_file + "..."
					plasmid_length=ExtractFastaPlasmids(gbfile,fasta_file,plasmid_length)
				else:
					print "Fasta Found! No conversion needed"
					# if there was a previous gb->fasta conversion
					fasta_path = os.path.join(dirname,'fasta',os.path.splitext(filename)[0])+'.fasta'
					plasmid_length=SequenceLengthFromFasta(fasta_file,plasmid_length,fasta_path)
					fasta_file = fasta_path + ".temp"
					dblist.append(fasta_file)
			elif any (x in filename for x in [".fas",".fasta",".fna",".fsa", ".fa"]):
				#if it is a fasta file extension
				FolderExist(os.path.join(dirname + "fasta"))	
				pct+=1
				print "Plasmid file (.fasta) found: "+ filename
				print "#:"+str(pct)
				fasta_file = os.path.join(dirname,os.path.splitext(filename)[0]) + os.path.splitext(filename)[1]
				fasta_path = os.path.join(dirname,'fasta/',os.path.splitext(filename)[0])+'.fasta'
				print "Fasta Found! No conversion needed (.fasta)"
				#copyfile(fasta_file, fasta_path)
					# if there was a previous .gb->.fasta conversion
				plasmid_length=SequenceLengthFromFasta(fasta_file,plasmid_length,fasta_path)
				fasta_file = os.path.join(fasta_path + ".temp")
				dblist.append(fasta_file)
			else:
				print "File extension " +os.path.splitext(filename)[1]+ " not recognized! Are you sure this is a fasta or gb file? Extensions autorized are .gb, .fas, .fasta, .fna, .fsa or .fa"
				
	print "========================================================================="
	print 

## CONCATENATES ALL PLASMID FASTA INTO A SINGLE DB ######################################
## This avoids the mapping of reads to a "worst matching" plasmid########################
## Instead reads will now be matched against the best scoring plasmid in the entire db###
	main_db = FastaConcatenation(dblist) 
	return main_db

def mapper(idx_file,reads_file,threads,max_k, sam_file,maindb_path):
	btc ='bowtie2 -x '+idx_file+' -U '+reads_file+' -p ' +threads+ ' -k '+ max_k + ' -5 15 -S '+ sam_file
	print "1) " + btc
	proc1=Popen(btc, stdout = PIPE, stderr = PIPE, shell=True)
	out,err= proc1.communicate()
	print err
	regex_match=re.search('[\d]{1}[.]{1}[\d]{2}% overall alignment rate',err)
	alignment_rate=regex_match.group(0).split('%')[0]
	if alignment_rate>0.00:
		print "2) " + 'samtools faidx '+maindb_path
		call('samtools faidx '+maindb_path, shell=True)
		bam_file = sam_file[:-3]+'bam'
		print("3) " + 'samtools view -b -S -t '+maindb_path+'.fai'+' -@ ' +threads+' -o '+bam_file+' '+sam_file)
		call('samtools view -b -S -t '+maindb_path+'.fai' +' -@ ' +threads+' -o '+bam_file+' '+sam_file, shell=True)
		sorted_bam_file = bam_file[:-3]+'sorted.bam'
		print("4) "+ 'samtools sort'+ ' -@ ' +threads+' -o '+sorted_bam_file+ ' ' + bam_file)
		call('samtools sort'+ ' -@ ' +threads+' -o '+sorted_bam_file+ ' ' + bam_file, shell=True)
		print("5)" + 'samtools index '+sorted_bam_file)
		call('samtools index '+sorted_bam_file, shell=True)
		print("6) " + 'samtools depth '+sorted_bam_file)
		depth_file = sorted_bam_file+'_depth.txt'
		print("Creating coverage Depth File: " + depth_file)
		proc2=Popen('samtools depth '+sorted_bam_file + ' > '+ depth_file, stdout = PIPE, stderr = PIPE, shell=True)
		out2,err2 = proc2.communicate()
		return depth_file


## Argparser arguments

parser = argparse.ArgumentParser(description="Outputs a coverage percentage for each Plasmid gbk in PlasmidDir using the reads presented in the directory structure in ReadsDir")
parser.add_argument('-p','--plasmid', dest='plasmid_dir', required=True, help='Provide the path to the directory containing plasmid fastas')
parser.add_argument('-r','--read', dest='read_dir', required=True, help='Provide the path to the directory containing reads fastas')
parser.add_argument('-t', '--threads', dest='threads', default="1", help="Specify the number of threads to be used by bowtie2")
parser.add_argument('-k', "--max_align", dest="max_align", help="Specify the maximum number of alignments possible for each read. This options changes -k parameter of Bowtie2. By default this script will set -k to the number of fastas in reference directory (e.g. if you have 3 reference sequences the number of max_align allowed will automatically be set to 3.")
parser.add_argument('-o','--output', dest='output_name', default="plasmid_db_out", help='Specify the output name you wish. No need for file extension!')
parser.add_argument('-c','--cutoff', dest='cutoff_number', default = "0.00", help='Specify the cutoff for percentage of plasmid coverage that reads must have to be in the output. This should be a number between 0.00-1.00')
#parser.add_argument('-g', '--graphs', dest='graphical', help='This option enables the output of graphical visualization of the % coverage in each plasmid. This options is intended to provide the user a better criteria for defining the optimal cut-off value (-c option)')
args = parser.parse_args()

## Lists and dictionaries

plasmid_length={}
strain_list=[]
pidx2name={}
dblist=[]
sam_dict={}


## Process plasmids references into a single fasta

maindb = PlasmidProcessing(dblist,args.plasmid_dir,plasmid_length)
maindb_path = os.path.join(args.plasmid_dir + "fasta/" + maindb)

## Deletes temporary fastas created during PlasmidProcessing function
deltemp(os.path.join(args.plasmid_dir + "fasta/"))

##Create Bowtie Idx files for plasmid references
idx_file=CreateBowtieIdx(maindb)

### READS#########################
output_txt = open(args.output_name +".txt", "w")
master_keys = []
trace_list = []
counter=0 # counter used to control output file
for dirname, dirnames, filenames in os.walk(args.read_dir):
	for subdirname in dirnames:
		for dirname2, dirnames2, filenames2 in os.walk(os.path.join(dirname,subdirname)):
			for filename in filenames2:
				if filename.find('fastq')!=-1:
					fn = filename.split('.')[0]
					strain_list.append(fn)
					print 
					print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
					print
					print "Filename :"+ filename
					print datetime.fromtimestamp(time()).strftime('%Y-%m-%d %H:%M:%S')
					print "Mapping "+ filename+" vs "+ maindb_path
					sam_file= dirname2+'/'+args.output_name+'_'+subdirname+'.sam'
					reads_file=os.path.join(dirname2,filename)
					threads = args.threads
					max_k=alignmaxnumber(args.max_align,dblist)
					depth_file=mapper(idx_file,reads_file,threads,max_k, sam_file,maindb_path)
					
		
			## Compute descritptive statistics and prints to tabular txt file
			
			Percentage_BasesCovered, Mean = DepthFileReader(depth_file, plasmid_length)
			sorted_percCoverage_dic = sorted(Percentage_BasesCovered.items(), key=operator.itemgetter(1), reverse=True)
			tmp_list_k = []
			tmp_list_v = []
			list_all_k = []
			list_all_v = []
			if 0.00<=float(args.cutoff_number)<=1.00: 
			## Filtered output
				if counter == 0:
					output_txt.write("NOTE: outputed results for plasmids with more than "+ args.cutoff_number + " coverage.\n\n")
					counter=1
				for k,v in sorted_percCoverage_dic:
					if v >= float(args.cutoff_number):
						tmp_list_k.append(k)
						tmp_list_v.append(v)
					if k not in master_keys:
						master_keys.append(k)
					list_all_v.append(v)
					list_all_k.append(k)
				## COVERAGE PERCENTAGE ##
				output_txt.write(fn + "\t" + ("\t").join(tmp_list_k) + "\nCoverage Percentage\t")
				for element in tmp_list_v:
					output_txt.write(str(element) +"\t")
				trace = go.Bar(x=list_all_k, y=list_all_v, name=fn)
				trace_list.append(trace)
				## MEAN ##
				output_txt.write("\nMean mapping depth\t")
				for element in tmp_list_k:
					output_txt.write(str(Mean[element]) + "\t")

				output_txt.write("\n")
			#Standard output
			else:
				print "cutoff value out of bounds. cutoff value must be between 0 and 1 since it is a probability"

output_txt.close()

### Graphical outputs ###

bar_plot(trace_list, float(args.cutoff_number), master_keys)