import sys
import os
import re
import numpy as np
import prettytable
#import matplotlib.pyplot as plt
from subprocess import Popen, PIPE, call
from collections import defaultdict
from Bio import SeqIO
from shutil import copyfile
from time import time
from datetime import datetime

try:
	plasmid_dir = sys.argv[1]
	read_dir = sys.argv[2]

except IndexError:
	print "Usage: PlasmidCoverage.py <PlasmidDir> <ReadsDir>"
	print "Outputs a coverage percentage for each Plasmid gbk in PlasmidDir using the reads presented in the directory structure in ReadsDir"
	print "jcarrico@fm.ul.pt - 23/10/2013"
	raise SystemExit

plasmid_idx_list=[]
plasmid_length={}
output_values=defaultdict(dict)
strain_list=[]
plasmid_list=[]
pidx2name={}

def ExtractPlasmidNameFromFasta(fasta_file):
	print fasta_file
	if_handle=open(fasta_file,'r')
	fastadata=SeqIO.read(if_handle,"fasta")
	pieces = fastadata.description.split('|')
	name=pieces[len(pieces)-1].lstrip()
	name_pieces = name.split(' ')
	rem=re.search('plasmid ([\w-]+)',name)
	PlasmidName=name_pieces[0][0]+'_'+name_pieces[1]+'_'+rem.groups()[0]
	if_handle.close()
	return PlasmidName

def ExtractFastaPlasmids(gbkfile,fastafile,plasmid_length):
	if_handle=open(gbkfile,'r')
	gbkdata=SeqIO.read(if_handle, "genbank")
	out_handle=open(fasta_file,'w')
	out_handle.write('>'+gbkdata.description+'\n')
	out_handle.write(gbkdata.seq.tostring())
	if_handle.close() 
	out_handle.close()
	plasmid_name=ExtractPlasmidNameFromFasta(fasta_file)
	plasmid_length[plasmid_name]=len(gbkdata.seq.tostring())
	return plasmid_length
	print " Wrote fasta file: "+fastafile

def SequenceLengthFromFasta(fasta_file,plasmid_length):
	if_handle=open(fasta_file,'r')
	fastadata=SeqIO.read(if_handle,"fasta")
	plasmid_name=ExtractPlasmidNameFromFasta(fasta_file)
	plasmid_length[plasmid_name]=len(fastadata.seq.tostring())
	if_handle.close()
	return plasmid_length 

def CreateBowtieIdx(filename,pidx2name):
	idx_file=os.path.join(dirname,'bowtie2idx/',os.path.splitext(filename)[0])+'.idx'
	fasta_file = os.path.join(dirname,'fasta/',os.path.splitext(filename)[0])+'.fasta'
	if not(os.path.exists(idx_file)):
		print "Creating " + idx_file 
		#Create bowtie index 			
		call('bowtie2-build -q '+ fasta_file +' '+idx_file, shell=True)
		plasmid_name=ExtractPlasmidNameFromFasta(fasta_file)
		pidx2name[idx_file]=plasmid_name
	return pidx2name

#TODO: 
#1) Set a minnimun read threshold to avoid computation on less than x reads mapped
#2) use the -p option set for the number of cores DONE: Set for 3 cores in Laptop. Change to 6 in dawkins.
#3)correct check for idx files
#4) delete unsorted BAM files
#5) Create directory structure to store all comparisons.


############# PLASMIDS ##################
print "========================================================================="
print "Processing Plasmids in "+ plasmid_dir
pct=0;
for dirname, dirnames, filenames in os.walk(plasmid_dir):
	for filename in filenames:
		#if it is a genebank file
		if filename.find('gb')!=-1:
			pct+=1
			print "Plasmid file found:"+ filename
			print "#:"+str(pct)
			gbfile = os.path.join(dirname, filename)
			fasta_file = os.path.join(dirname,'fasta/',filename[:-len('.gb')])+'.fasta'

			if not(os.path.exists(fasta_file)):
				print "Fasta file not file. Converting..."
				#If fasta not found Transform gbk to fasta
				print "Creating fasta: " + fasta_file + "..."
				plasmid_length=ExtractFastaPlasmids(gbfile,fasta_file,plasmid_length)
				#plasmid_idx_list=CreateBowtieIdx(filename,plasmid_idx_list,pidx2name)

			else:
				print "Fasta Found! No conversion needed"
				# if there was a previous gb->fasta conversion
				plasmid_length=SequenceLengthFromFasta(fasta_file,plasmid_length)
				#plasmid_idx_list=CreateBowtieIdx(filename,plasmid_idx_list,pidx2name)
			print
			pidx2name=CreateBowtieIdx(filename,pidx2name) 
		elif filename.find('fna')!=-1:
			pct+=1
			print "Plasmid file found:"+ filename
			print "#:"+str(pct)

			fasta_file = os.path.join(dirname,'fasta/',os.path.splitext(filename)[0])+'.fasta'

			if not(os.path.exists(fasta_file)):
				# If it is not a genebank it needs to be a fasta format
				original_fasta_filename = os.path.join(dirname, filename)
	
				copyfile(original_fasta_filename,fasta_file)

				plasmid_length=SequenceLengthFromFasta(fasta_file,plasmid_length)
				#plasmid_idx_list=CreateBowtieIdx(filename,plasmid_idx_list,pidx2name)
			else:
				print "Fasta Found! No conversion needed"
				# if there was a previous gb->fasta conversion
				plasmid_length=SequenceLengthFromFasta(fasta_file,plasmid_length)
				#plasmid_idx_list=CreateBowtieIdx(filename,plasmid_idx_list,pidx2name)
			print
			pidx2name=CreateBowtieIdx(filename,pidx2name)

plasmid_name_list=pidx2name.values()
plasmid_idx_list=pidx2name.keys()
print "========================================================================="
print 

### READS#########################
for dirname, dirnames, filenames in os.walk(read_dir):
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
					for plasmid_idx in plasmid_idx_list:

						plasmid_name=pidx2name[plasmid_idx]
						print
						print "######################################"
						print datetime.fromtimestamp(time()).strftime('%Y-%m-%d %H:%M:%S')
						print "Mapping "+ filename+" vs "+ plasmid_name
						sam_file= dirname2+'/'+plasmid_idx.split('/')[1][:-4]+'_'+subdirname+'.sam'
						reads_file=os.path.join(dirname2,filename)
						btc ='bowtie2 -x '+plasmid_idx+' -U '+reads_file+' -p 3 -5 15 -S '+sam_file
						print "1) " + btc
						proc1=Popen(btc, stdout = PIPE, stderr = PIPE, shell=True)
						out,err= proc1.communicate()
						print err
						regex_match=re.search('[\d]{1}[.]{1}[\d]{2}% overall alignment rate',err)
						alignment_rate=regex_match.group(0).split('%')[0]
						print alignment_rate

						output_values[fn,plasmid_name]['alignment_rate']=alignment_rate

						if err.find('0.00% overall alignment rate')<0:
							print "2) " + 'samtools faidx '+fasta_file
							call('samtools faidx '+fasta_file, shell=True)

							bam_file = sam_file[:-3]+'bam'
							print "3) " + 'samtools view -b -S -t '+fasta_file+'.fai -o '+bam_file+' '+sam_file
							call('samtools view -b -S -t '+fasta_file+'.fai -o '+bam_file+' '+sam_file, shell=True)
							call('rm -f '+sam_file,shell=True)

							sorted_bam_file = bam_file[:-3]+'sorted'
							print "4) "+ 'samtools sort '+ bam_file +' '+sorted_bam_file
							call('samtools sort '+ bam_file +' '+sorted_bam_file, shell=True)

							print "5)" + 'samtools index '+sorted_bam_file+'.bam'
							call('samtools index '+sorted_bam_file+'.bam', shell=True)

							print "6) " + 'samtools depth '+sorted_bam_file+'.bam'
							depth_file = sorted_bam_file+'_depth.txt'
							print "Creating coverage Depth File: " + depth_file
							proc2=Popen('samtools depth '+sorted_bam_file+'.bam >'+ depth_file, stdout = PIPE, stderr = PIPE, shell=True)
							print "done"
							out2,err2 = proc2.communicate()
							#print "--out2--"
							#print out2
							#print "--out2--"
							#print err2
							depth_info=np.loadtxt(depth_file,usecols=(1,2))
							output_values[fn,plasmid_name]['n_bases_covered']=len(depth_info)
							output_values[fn,plasmid_name]['mean_depth']=np.mean(depth_info[:,1])
							output_values[fn,plasmid_name]['median_depth']=np.median(depth_info[:,1])
							output_values[fn,plasmid_name]['min_depth']=np.min(depth_info[:,1])
							output_values[fn,plasmid_name]['max_depth']=np.max(depth_info[:,1])
							output_values[fn,plasmid_name]['plasmid_length']=plasmid_length[plasmid_name]

						else:
							print "No reads were aligned"
							output_values[fn,plasmid_name]['n_bases_covered']='na'
							output_values[fn,plasmid_name]['alignment_rate']='na'
							output_values[fn,plasmid_name]['mean_depth']='na'
							output_values[fn,plasmid_name]['median_depth']='na'
							output_values[fn,plasmid_name]['min_depth']='na'
							output_values[fn,plasmid_name]['max_depth']='na'
							output_values[fn,plasmid_name]['plasmid_length']=plasmid_length[plasmid_name]
						print "######################################"
						print


results_html_file="Results.html"
fh=open(results_html_file,'w')

first_row_header=['----']

for plasmid in plasmid_name_list:
	first_row_header+=[plasmid]

Result_table_cov_perc=prettytable.PrettyTable(first_row_header)

for strain in strain_list:
	row=[ strain ]
	for plasmid in plasmid_name_list:
		#print "~.....~"
		if output_values[strain,plasmid]['n_bases_covered']=='na':
			percent_mapped='na'
			row+=[percent_mapped]
		else:
			percent_mapped=float(output_values[strain,plasmid]['n_bases_covered'])/float(output_values[strain,plasmid]['plasmid_length'])
			row+=["%2.3f" % percent_mapped]
	Result_table_cov_perc.add_row(row)
print "==============>   Coverage Percentage Table <==============" 
print Result_table_cov_perc
fh.write("<H1>Coverage Percentage Table</H1>")
fh.write(Result_table_cov_perc.get_html_string())

Result_table_mean_depth=prettytable.PrettyTable(first_row_header)

for strain in strain_list:
	row=[ strain ]
	for plasmid in plasmid_name_list:
		if output_values[strain,plasmid]['mean_depth']=='na':
			mean_depth='na'
			row+=[mean_depth]
		else:
			mean_depth=float(output_values[strain,plasmid]['mean_depth'])
			row+=["%2.3f" % mean_depth]
	Result_table_mean_depth.add_row(row)
print "==============>   Mean mapping depth Table <=============="
print Result_table_mean_depth

fh.write("<H1>Mean mapping depth Table</H1>")
fh.write(Result_table_mean_depth.get_html_string())

Result_table_median_depth=prettytable.PrettyTable(first_row_header)

for strain in strain_list:
	row=[ strain ]
	for plasmid in plasmid_name_list:
		if output_values[strain,plasmid]['median_depth']=='na':
			median_depth='na'
			row+=[median_depth]
		else:
			median_depth=float(output_values[strain,plasmid]['median_depth'])
			row+=["%2.3f" % median_depth]
	Result_table_median_depth.add_row(row)
print "==============>   Median mapping depth Table <=============="
print Result_table_median_depth
fh.write("<H1>Median mapping depth Table</H1>")
fh.write(Result_table_median_depth.get_html_string())

fh.close()

