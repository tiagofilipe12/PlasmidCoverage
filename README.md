#PlasmidCoverage_Sdb.py
Initial tests in finding plasmid coverages from reads and inferring plasmid contents

For now, this developing branch of this repo provide two versions of the same script:

* Version *PlasmidCoverage.py* maps against only one plasmid sequence.

* Version *PlasmidCoverage_Sdb.py* allows users to compile multiple fastas into one or having one multi-fasta file to calculate plasmid coverage. 

In both versions, .gb files are also supported. These files will be converted into fasta files in a subdirectory named "fasta" within the directory the user specifies for plasmid.

###Options for PlasmidCoverage_Sdb.py:

**'-p'**,**'--plasmid'**, dest='plasmid_dir', required=True, help='Provide the path to the directory containing plasmid fastas.'

**'-r'**,**'--read'**, dest='read_dir', required=True, help='Provide the path to the directory containing reads fastas.'

**'-t'**, **'--threads'**, dest='threads', default="1", help="Specify the number of threads to be used by bowtie2."

**'-k'**, **"--max_align"**, dest="max_align", help="Specify the maximum number of alignments possible for each read. This option changes -k parameter of Bowtie2. By default this script will set -k to the number of fasta files in reference directory (e.g. if you have 3 reference sequences the number of max_align allowed will automatically be set to 3). So, if you have a single multi-fasta -k will be set to 1."

**'-o'**,**'--output'**, dest='output_name', default="plasmid_db_out", help='Specify the output name you wish. There is no need of file extension.')

**'-c'**,**'--cutoff'**, dest='cutoff_number', help='Specify the cutoff for percentage of plasmid coverage that reads must have to be in the output. This should be a number between 0-1.')

---

####- Acessory tools

**plasmid_or_not.py**

Currently this tool is separated from the main script (PlasmidCoverage_Sdb.py). This script is intended to separate reads from plasmids (or other type of sequences that have a database in fasta format) from reads for mixed libraries of chromosomal + plasmid reads.

######Options for plasmid_or_not.py

**'-p'**,**'--plasmid'**, dest='plasmid', nargs='+', required=True, help='Provide the plasmid fastas'

**'-r'**,**'--read'**, dest='reads', nargs='+', required=True, help='Provide the path to the directory containing reads fastas'

**'-t'**, **'--threads'**, dest='threads', default="1", help="Specify the number of threads to be used by bowtie2"

**'-o'**,**'--output'**, dest='output_name', required=True, help='Specify the output name you wish. No need for file extension! Output file will be a fasta.'

**'-unmap'**,**'--unmapped'**, dest='unmapped_reads', action='store_true', help='By default this script attempts to save sequences available in the provided read files. If you want to save the reads that do not belong to plasmids, use this option.'



