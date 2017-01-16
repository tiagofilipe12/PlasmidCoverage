#!/usr/bin/env python

## Last update: 14/1/2017

#import sys
import argparse
import os
import re
#import numpy as np
#import prettytable
#import matplotlib.pyplot as plt
from subprocess import Popen, PIPE, call
#from collections import defaultdict
from Bio import SeqIO
#from shutil import copyfile
from time import time
from datetime import datetime


#TODO: 
#1) Set a minnimun read threshold to avoid computation on less than x reads mapped
#3)correct check for idx files
#4) delete unsorted BAM files
#5) Create directory structure to store all comparisons.

#except IndexError:
#	print "Usage: PlasmidCoverage.py <PlasmidDir> <ReadsDir>"
#	print "Outputs a coverage percentage for each Plasmid gbk in PlasmidDir using the reads presented in the directory structure in ReadsDir"
#	print "jcarrico@fm.ul.pt - 23/10/2013"
#	raise SystemExit

parser = argparse.ArgumentParser(description="Outputs a coverage percentage for each Plasmid gbk in PlasmidDir using the reads presented in the directory structure in ReadsDir")
parser.add_argument('-p','--plasmid', dest='plasmid_dir', required=True, help='Provide the path to the directory containing plasmid fastas')
parser.add_argument('-r','--read', dest='read_dir', required=True, help='Provide the path to the directory containing reads fastas')
parser.add_argument('-t', '--threads', dest='threads', default="1", help="Specify the number of threads to be used by bowtie2")
#add an option to filter low coverage reads from samtools depth
#argument if depth files already exist for the alignment
parser.add_argument('-k', "--max_align", dest="max_align", help="Specify the maximum number of alignments possible for each read. This options changes -k parameter of Bowtie2. By default this script will set -k to the number of fastas in reference directory (e.g. if you have 3 reference sequences the number of max_align allowed will automatically be set to 3.")
#Define the max number of alingments to the number of sequences
parser.add_argument('-o','--output', dest='output_name', default="plasmid_db_out", help='Specify the output name you wish. No need for file extension!')
args = parser.parse_args()

plasmid_length={}
#output_values=defaultdict(dict)
strain_list=[]
pidx2name={}
dblist=[]
sam_dict={}

def alignmaxnumber(max_align, dblist):
	if max_align:
		k_value=max_align
	else:
		k_value = str(len(dblist))
	return k_value

def FolderExist(directory):
#	if not directory.endswith("/"):
#		directory = directory + "/"
	if not os.path.exists(os.path.join(directory)):		
		os.makedirs(os.path.join(directory))
		print os.path.join(directory) + " does not exist. One will be created..."
	else:
		print os.path.join(directory) + " exists!"	

def FastaDict(fasta_file):
	if_handle=open(fasta_file,'r')
	x = 0
	sequence = ""
	fasta_dic = {}
	problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/"]

	for line in if_handle:
		if len(line) > 0: 
			line = line.splitlines()[0]  
		if x >= 0 and line.startswith(">"):
			PlasmidName = line[1:]
			for char in problematic_characters:
					PlasmidName = PlasmidName.replace(char, '_')
			x+=1
		elif x == 0 and not line.startswith(">"):
			print "Is this a fasta file? " + fasta_file
			raise SystemExit
		else:
			print("Passing sequence itself...")
			sequence += line
			fasta_dic[PlasmidName] = sequence
	if_handle.close()
	return fasta_dic

def SequenceLengthFromFasta(fasta_file,plasmid_length,fasta_path):
	fasta_dic = FastaDict(fasta_file)
	out_handle = open(fasta_path + ".temp", 'w')
	for key in fasta_dic:
		plasmid_length[key]=len(fasta_dic[key])
		out_handle.write('>' + key + '\n' + fasta_dic[key] + '\n')
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
	return plasmid_length
	print(" Wrote fasta file: "+ fastafile)

def CreateBowtieIdx(filename):
	dirname = args.plasmid_dir
	if not os.path.exists(os.path.join(dirname + "bowtie2idx")):		
		os.makedirs(os.path.join(dirname + "bowtie2idx"))
	else:
		pass	
	idx_file=os.path.join(dirname,'bowtie2idx',os.path.splitext(filename)[0])+'.idx'
	fasta_file = os.path.join(dirname,'fasta',os.path.splitext(filename)[0])+'.fasta'
	if not(os.path.exists(idx_file)):
		print("Creating " + idx_file)
		#Create bowtie index 			
		call('bowtie2-build -q '+ fasta_file +' '+idx_file, shell=True)		#convert to popen
	return idx_file

def FastaConcatenation(dblist):
	output_filename = args.output_name + ".fasta"
	dirname = os.path.join(args.plasmid_dir, "fasta", "")
	print(dirname)
	if os.path.isfile(dirname + output_filename):
		print args.output_name + ".fasta already exists. Overriding file..."
	main_filename = args.output_name + ".fasta"
	print "Saving to: " + main_filename
	concat=Popen("cat " + ' '.join(dblist) + "> " + dirname + main_filename, stdout = PIPE, stderr = PIPE, shell=True)
	stdout,stderr= concat.communicate()
	print stderr
	return main_filename

## function to delete temporary fasta files
def DelTemp(directory):
	files = os.listdir(directory)
	print "Deleting temporary fasta files..."
	for f in files:
		if f.endswith('.temp'):
			os.remove(os.path.join(directory,f))

## Creates an alignment for each .sam file (per reads file) and outputs a dictionary with keys = reads and values = 
## Still not implemented
def SamDictMultipleHits(samfile):
	input_sam = open(samfile,"r+")
	for line in input_sam:
		if not line.startswith("@"):
			tab_split = line.split("\t")
			Read = tab_split[0]
			Flag = tab_split[1]
			Ref = tab_split[2]
			Leftmost_pos = tab_split[3]
			Seq_lenght = tab_split[9]
			if Ref != "*":
				sam_dict[Read]= [Flag, Ref, Leftmost_pos, Seq_lenght]
	return sam_dict

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
#		Mean[ref] = np.mean(depth_dic_numreads[ref])
#		Median[ref] = np.median(depth_dic_numreads[ref])
#		Minimum[ref] = np.min(depth_dic_numreads[ref])
#		Maximum[ref] = np.max(depth_dic_numreads[ref])
	return Percentage_BasesCovered, Mean #, Median, Minimum, Maximum



############# PLASMIDS ##################
def PlasmidProcessing(dblist,plasmids_path,plasmid_length):
	print "========================================================================="
	print "Processing Plasmids in "+ plasmids_path
	pct=0;
	for dirname, dirnames, filenames in os.walk(plasmids_path, topdown=True): 
		dirnames[:] = [d for d in dirnames if d not in ["bowtie2idx", "fasta"]]
		for filename in filenames:
			#if it is a genebank file
			if '.gb' in filename:
				pct+=1
				print "Plasmid (.gb) file found: "+ filename
				print "#:"+str(pct)
				gbfile = os.path.join(dirname, filename)
				FolderExist(os.path.join(dirname + "fasta"))			
				fasta_file = os.path.join(dirname,'fasta/',filename[:-len('.gb')])+'.fasta'
				if not(os.path.exists(fasta_file)):
					print "Fasta file not file. Converting..."
					#If fasta not found Transform .gb to .fasta
					print "Creating fasta: " + fasta_file + "..."
					plasmid_length=ExtractFastaPlasmids(gbfile,fasta_file,plasmid_length)
				else:
					print "Fasta Found! No conversion needed"
					# if there was a previous gb->fasta conversion
					fasta_path = os.path.join(dirname,'fasta/',os.path.splitext(filename)[0])+'.fasta'
					plasmid_length=SequenceLengthFromFasta(fasta_file,plasmid_length,fasta_path)
					fasta_file = fasta_path + ".temp"
			elif any (x in filename for x in [".fas",".fasta",".fna",".fsa", ".fa"]):
				#if it is a fasta file extension
				FolderExist(dirname + "fasta/")	
				pct+=1
				print "Plasmid file (.fasta) found: "+ filename
				print "#:"+str(pct)
				fasta_file = os.path.join(dirname,os.path.splitext(filename)[0])+'.fasta'
				fasta_path = os.path.join(dirname,'fasta/',os.path.splitext(filename)[0])+'.fasta'
				print "Fasta Found! No conversion needed (.fasta)"
					# if there was a previous .gb->.fasta conversion
				plasmid_length=SequenceLengthFromFasta(fasta_file,plasmid_length,fasta_path)
				fasta_file = fasta_path + ".temp"
			else:
				print "File extension not recognized! Are you sure this is a fasta or gb file? Extensions autorized are .gb, .fas, .fasta, .fna, .fsa or .fa"
				
			dblist.append(fasta_file)


	print "========================================================================="
	print 

## CONCATENATES ALL PLASMID FASTA INTO A SINGLE DB ######################################
## This avoids the mapping of reads to a "worst matching" plasmid########################
## Instead reads will now be matched against the best scoring plasmid in the entire db###
	main_db = FastaConcatenation(dblist) 
	return main_db



## Process plasmids references into a single fasta

maindb = PlasmidProcessing(dblist,args.plasmid_dir,plasmid_length)
maindb_path = args.plasmid_dir + "fasta/" + maindb

## Deletes temporary fastas created during PlasmidProcessing function
DelTemp(args.plasmid_dir + "fasta/")

##Create Bowtie Idx files for plasmid references

idx_file=CreateBowtieIdx(maindb)

### READS#########################
output_txt = open(args.output_name +".txt", "w")
for dirname, dirnames, filenames in os.walk(args.read_dir):
	print(dirnames)
	for subdirname in dirnames:
		#print(subdirname) If fastq files are not in subdirectories of main directory the script won't do anything... this must be handled
		for dirname2, dirnames2, filenames2 in os.walk(os.path.join(dirname,subdirname)):
			#print filenames2
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
					btc ='bowtie2 -x '+idx_file+' -U '+reads_file+' -p ' +args.threads+ ' -k '+ alignmaxnumber(args.max_align,dblist) + ' -5 15 -S '+ sam_file
					print "1) " + btc
					proc1=Popen(btc, stdout = PIPE, stderr = PIPE, shell=True)
					out,err= proc1.communicate()
					print err
					sam_parser = SamDictMultipleHits(sam_file)
					regex_match=re.search('[\d]{1}[.]{1}[\d]{2}% overall alignment rate',err)
					alignment_rate=regex_match.group(0).split('%')[0]
					print alignment_rate

					if alignment_rate>0.00:
						print "2) " + 'samtools faidx '+maindb_path
						call('samtools faidx '+maindb_path, shell=True)

						bam_file = sam_file[:-3]+'bam'
						print("3) " + 'samtools view -b -S -t '+maindb_path+'.fai'+' -@ ' +args.threads+' -o '+bam_file+' '+sam_file)
						call('samtools view -b -S -t '+maindb_path+'.fai' +' -@ ' +args.threads+' -o '+bam_file+' '+sam_file, shell=True)
						SamDictMultipleHits(sam_file)

						sorted_bam_file = bam_file[:-3]+'sorted.bam'
						print("4) "+ 'samtools sort'+ ' -@ ' +args.threads+' -o '+sorted_bam_file+ ' ' + bam_file)
						call('samtools sort'+ ' -@ ' +args.threads+' -o '+sorted_bam_file+ ' ' + bam_file, shell=True)

						print("5)" + 'samtools index '+sorted_bam_file)
						call('samtools index '+sorted_bam_file, shell=True)

						print("6) " + 'samtools depth '+sorted_bam_file)
						depth_file = sorted_bam_file+'_depth.txt'
						print("Creating coverage Depth File: " + depth_file)
						proc2=Popen('samtools depth '+sorted_bam_file + ' > '+ depth_file, stdout = PIPE, stderr = PIPE, shell=True)
						print("done")
						out2,err2 = proc2.communicate()
		
		## Compute descritptive statistics and prints to tabular txt file
			
			Percentage_BasesCovered, Mean = DepthFileReader(depth_file, plasmid_length)
			output_txt.write(fn + "\t" + ("\t").join(Percentage_BasesCovered) + "\nCoverage Percentage\t")
			#COVERAGE PERCENTAGE
			for element in Percentage_BasesCovered:
				#print Percentage_BasesCovered, "Percentage"
				output_txt.write(str(Percentage_BasesCovered[element]) +"\t")
			#MEAN#
			output_txt.write("\nMean mapping depth\t")
			for element in Percentage_BasesCovered:
				#print Percentage_BasesCovered, "Mean" 
				output_txt.write(str(Mean[element]) + "\t")
			#MEDIAN#
#			output_txt.write("\nMedian mapping depth\t")
#			for element in Reference_list:
#				output_txt.write(str(Median[element]) + "\t")
			#MAXIMUM#
#			output_txt.write("\nMaximum mapping depth\t")
#			for element in Reference_list:
#				output_txt.write(str(Maximum[element]) + "\t")
			#MINIMUM#
#			output_txt.write("\nMinimum mapping depth\t")
#			for element in Reference_list:
#				output_txt.write(str(Minimum[element]) + "\t")
			output_txt.write("\n")

output_txt.close()