# PlasmidCoverage_Sdb.py
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

