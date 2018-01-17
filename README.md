# PlasmidCoverage

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/4ecf8dfe775746f4bc5f3d154a7207df)](https://www.codacy.com/app/tiagofilipe12/PlasmidCoverage?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=tiagofilipe12/PlasmidCoverage&amp;utm_campaign=Badge_Grade)

## About

PlasmidUNCover.py is a script that is intended to allow users to map reads 
against a plasmid database (provided in [indexes folder](https://github.com/tiagofilipe12/PlasmidCoverage/releases/download/v1.0.0/indexes.tar.gz))
. Then, the resulting accession numbers that have a percentage of covered 
base pairs higher than the specified cutoff value will be reported into a 
.json file which may be imported by [Plasmid Atlas](http://www.patlas.site).

#### Other uses

It is also possible to map against any reference fasta, since it will 
automatically recognize the option provided and construct the required index 
files for that reference you want.

## How to install

* First install [bowtie2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/) (v2.2.9 or higher) and 
[samtools](https://sourceforge.net/projects/samtools/files/samtools/) (v.1.3.1 or 
higher).

### Using pipy (recommended)

* `pip3 install plasmiduncover`

* Download [indexes folder](https://github.com/tiagofilipe12/PlasmidCoverage/releases/download/v1.0.0/indexes.tar.gz)
and uncompress it.

### Using github release

* Second, download the [latest release](https://github.com/tiagofilipe12/PlasmidCoverage/releases/tag/v1.0.4) of this script
(don't forget to download [indexes folder](https://github.com/tiagofilipe12/PlasmidCoverage/releases/download/v1.0.4/indexes.tar.gz)).

* `pip3 install -r requirements.txt`

* Uncompress the `indexes.tar.gz` file.

## Example run

### If you installed using git clone
`PlasmidUNCover.py -idx <path/to/indexes_folder> -r 
<path/to/reads_folder> -t 1 -o test -c 0.6`

Note that each read or pair of reads should be inside its own folder within 
the `path/to/reads_folder/`. This allows the user to run multiple samples at 
once, since `PlasmidUNCoverage.py` will crawl this directory to search for 
directories with samples.

For advanced users, you may use `PlasmidUNCoverage.py -h` for additional 
options.

## Dependencies

Note: ignore this if you have read the [How to install](#how-to-install) 
section.

* plotly - ```pip install plotly```
* termcolor
* bowtie2 (tested for version 2.2.9)
* samtools (tested for version 1.3.1)

You can now simply: ```pip3 install -r requirements.txt```

## Important note regarding the parsing of argument -r

You should provide the path to the directory containing the directories with the reads (.fastq files). For instance if you have something like `~/Reads/sample1/example.fastq` and `~/Reads/sample2/example2.fastq`, you should provide `~/Reads` as the argument for `-r` option.

### Options for PlasmidUNCoverage.py:

```
-p, --plasmid - Provide the path to the directory containing plasmid fastas.

-idx, --bowtie-index - Provide the path to bowtie index file

-r, --read - Provide the path to the directory containing reads fastas.

-t, --threads - Specify the number of threads to be used by bowtie2.

-2, --pair - Use this 
option if you have paired end reads. Paired end reads just have to be inside 
the same folder, without any other files inside it. Structure to read files 
should be something like: ./reads/read_folder/<with pairs inside it>.)

-k, --max_align - Specify the maximum number of alignments possible for each 
read. This option changes -k parameter of Bowtie2. By default this script will 
set -k to the number of fasta files in reference directory (e.g. if you have 3 
reference sequences the number of max_align allowed will automatically be set 
to 3). So, if you have a single multi-fasta -k will be set to 1.

-o, --output - Specify the output name you wish. There is no need of file 
extension.)

-c, --cutoff - Specify the cutoff for percentage of plasmid coverage that reads
 must have to be in the output. This should be a number between 0-1.')
```

---

#### - Accessory tools

**plasmid_or_not.py**

Currently this tool is separated from the main script (PlasmidCoverage_Sdb.py). This script is intended to separate reads from plasmids (or other type of sequences that have a database in fasta format) from reads for mixed libraries of chromosomal + plasmid reads.

###### Options for plasmid_or_not.py
```
-p, --plasmid - Provide the plasmid fastas

-r, --read - Provide the path to the directory containing reads fastas

-t, --threads - Specify the number of threads to be used by bowtie2

-o, --output - Specify the output name you wish. No need for file extension! 
Output file will be a fasta.

-unmap, --unmapped - By default this script attempts to save sequences 
available in the provided read files. If you want to save the reads that do not 
belong to plasmids, use this option.
```

**diffs_json.py**

Currently this tool is separated from the main script (PlasmidCoverage_Sdb.py). This script compares two json files with coverage percentage per gi (json files retrieved by PlasmidCoverage_Sdb.py.

##### Options for diffs_json.py
```
-i, --input_jsons - Provide the input json files of interest.

-c, --cutoff - Provide the cutoff value used for the minimum coverage allowed
 to be outputed to json files. This value should be equal in both files to compare.
```

