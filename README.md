# minitax
Taxonomic assignement of long metagenomic reads with minimap2

This tool uses the minimap2 aligner to map metagenomic reads to a selected database. 
The database can be a collection of genomes from NCBI for WGS reads, or for 16S reads the databse of EMU (https://gitlab.com/treangenlab/emu) or SILVA (https://www.arb-silva.de/).
Afterwards it tries to find the best taxonomic assignemnt for each read based on mapping quality scores (mapq) and potentially, platform-specifc CIGAR-scoring matrices. 
If there are still multiple alignments with same mapq values, minitax processes the alignments with either of the following three methods:
1. SpeciesEstimate: This methd does not counts the source read for each. It will use every alignment of reads that were kept after MAPQ filtering and normalize the counts based on the number alignments.
3. LCA: This will provide the Lowest Common Ancestor (tax.identity) for each read that were kept after MAPQ filtering.
4. BestAln: This will keep the best matching alignment for each read, based on cigar_scores
Finally it summarises the counts for the lowest taxonomic rank (species by default)

## Installing
the program can be downloaded from github using  
*git clone https://github.com/Balays/minitax.git*


## Configuration file:
The configuration file is tab-separated file, and should contain the follwing information and should look like this:

argument	value	step	description  
platform	ONT	both	Either: 'Illumina', 'PacBio' or 'ONT'  
db	'all_NCBI_genomes'	both	options: 'all_NCBI_genomes' or 'EMUdb'  
db.dir	/mnt/d/data/databases/all_NCBI_genomes	both	absolute path of database home directory  
project	project	optional	project identifier  
...  
*A sample configuration file is provided*

## Obtaining the database
*gunzip NCBI_genome_collection/all_NCBI_genomes_sequence_lengths.zip*
*minitax/download_ncbi_accessions.sh NCBI_genome_collection/all_NCBI_genomes_sequence_lengths.tsv*

## Running the tool:
Mappping part  
*minitax/minitax.sh minitax_config_allNCBI.txt*

Finding the best taxonomic assignement for each read  
*minitax/minitax.complete.R minitax_config_NCBI.txt*

## Outputs:
The outputs include a .tsv file containing the counts for each sample
And an .rds file containing a phyloseq-object

## R-pcakage Dependencies:
Rsamtools
readr
tidyr
dplyr
GenomicAlignments
data.table
stringi
stringr
phyloseq
future.apply



