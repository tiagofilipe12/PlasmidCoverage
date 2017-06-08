#!/usr/bin/env python

## Last update: 1/2/2017
## Author: T.F. Jesus
## This version runs with bowtie2 build 2.2.9, seqtk (1.2-r95-dirty) and samtools, with multithreading for bowtie2-build
## This script filters reads, given a plasmid database, allowing for example to separate both chromossomal reads from plasmid reads.

import argparse
from subprocess import Popen, PIPE
import os
import re

## Function to fix several issues that fasta header names can have with some programs 
def header_fix(input_header):
	problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/","[","]",":","{","}"]
	for char in problematic_characters:
		input_header=input_header.replace(char, '_')
	return input_header

## Function to create a master fasta file from several fasta databases. One fasta is enought though
def master_fasta(fastas, output_tag):
	master_fasta = open("master_fasta_" + output_tag + ".fas", "w")
	for filename in fastas:
		path=os.path.dirname(os.path.abspath(filename))
		fasta = open(filename,"r")
		for line in fasta:
			if line.startswith(">"):
				line = header_fix(line)
			master_fasta.write(line)
	master_fasta.close()
	return os.path.join(path, "master_fasta_" + output_tag + ".fas")

## Create bowtie index
def createbowtieidx(filename, threads):
	check =""
	idx_file=os.path.splitext(filename)[0]+'.idx'
	print idx_file
	fasta_file = os.path.splitext(filename)[0]+'.fas'
	print fasta_file
	list_idxfiles = [idx_file + ".1.bt2", idx_file + ".2.bt2", idx_file + ".3.bt2", idx_file + ".4.bt2", idx_file + ".rev.2.bt2", idx_file + ".rev.1.bt2"]
	for idx in list_idxfiles:
		if not os.path.isfile(idx): 
			check = "yes"
	if "yes" in check:
		print("Creating " + idx_file)
		#Create bowtie index 			
		p=Popen('bowtie2-build -q '+ fasta_file +' --threads ' +threads+ ' ' +idx_file, stdout = PIPE, stderr = PIPE, shell=True)
		p.wait()
	else:
		print idx_file + " already exists!"
	return idx_file

##Runs bowtie and samtools to retrieve either the plasmid filtered sequences or chromossomal sequences from reads making a unique file with plasmid or chromossomes reads only
def mapper(idx_file,read,threads,main_file, output_name, unmapped):
	## specify bam and sam file names
	print read
	sam_file = os.path.basename(read).split(".")[0] +"_"+ output_name + ".sam"
	bam_file = sam_file[:-3]+'bam'
	## Runs the three commands necessary to have only unmapped reads
	btc ='bowtie2 -x '+idx_file+' -U '+read+' -p ' +threads+ ' -5 15 -S '+ sam_file
	print "1) " + btc
	proc1=Popen(btc, stdout = PIPE, stderr = PIPE, shell=True)
	proc1.wait()
	err= proc1.communicate()
	print err
	regex_match=re.search('[\d]{1}[.]{1}[\d]{2}% overall alignment rate',err)
	alignment_rate=regex_match.group(0).split('%')[0]
	if alignment_rate>0.00:
		sf = "samtools faidx "+main_file
		print "2) "+sf
		proc2=Popen(sf, stdout = PIPE, stderr = PIPE, shell=True)
		proc2.wait()
		out,err= proc2.communicate()
		print err
		if unmapped:
			sv = 'samtools view -b -f 4 ' + sam_file + ' > '+bam_file 		## all reads but the ones in the database provided
		else:
			sv = 'samtools view -b -F 4 ' + sam_file + ' > '+bam_file		## all reads mapped against the database provided
		print "3) "+sv
		proc3=Popen(sv, stdout = PIPE, stderr = PIPE, shell=True)
		proc3.wait()
		out,err= proc3.communicate()
		print err

		## convert file to fasta format

		sbf = "samtools bam2fq " + bam_file+" | seqtk seq -A > "+sam_file[:-3] + "fas"
		proc4 = Popen(sbf, stdout = PIPE, stderr = PIPE, shell=True)
		proc4.wait()
		out,err= proc4.communicate()
		print err
		print "Saved reads to: " + str(sam_file[:-3]) + "fas"


def main():
	parser = argparse.ArgumentParser(description="")
	parser.add_argument('-p','--plasmid', dest='plasmid', nargs='+', required=True, help='Provide the plasmid fastas')
	parser.add_argument('-r','--read', dest='reads', nargs='+', required=True, help='Provide the path to the directory containing reads fastas')
	parser.add_argument('-t', '--threads', dest='threads', default="1", help="Specify the number of threads to be used by bowtie2")
	parser.add_argument('-o','--output', dest='output_name', required=True, help='Specify the output name you wish. No need for file extension!')
	parser.add_argument('-unmap','--unmapped', dest='unmapped_reads', action='store_true', help='By default this script attempts to save sequences available in the provided read files. If you want to save the reads that do not belong to plasmids, use this option.')
	args = parser.parse_args()

	if args.threads is None:
		threads = "1"
	else:
		threads = args.threads
	fastas = []
	for filename in args.plasmid:
		if any (x in filename for x in [".fas",".fasta",".fna",".fsa", ".fa"]):
			fastas.append(filename)
	for read in args.reads:
		main_file=master_fasta(fastas, args.output_name)
		idx_file=createbowtieidx(main_file, threads)
		mapper(idx_file,read,threads,main_file, args.output_name, args.unmapped_reads)

if __name__ == "__main__":
	main()