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
git clone https://github.com/Balays/minitax.git


## Configuration file:
The configuration file is tab-separated file, and should contain the follwing information and should look like this:

argument	value	step	description
platform	ONT	both	Either: 'Illumina', 'PacBio' or 'ONT'
db	'UHGV'	both	options: 'all_NCBI_genomes' or 'EMUdb' or 'UHGV'
db.dir	/mnt/d/data/databases/all_NCBI_genomes	both	absolute path of database home directory
project	Toti_UHGV	optional	project identifier
Vregion	V1_V2	optional	optional
indir	./fastq	minimap2	absolute path to directory of fastq files
outdir	NA	both	will be composed from platform, db, and project arguments, this will be the output directory (will be created and files here will be ovrewritten)
debug	FALSE	both	if there was a problem, try to show in which step
mm2_path	/mnt/d/ubuntu/programs/mm2-fast/minimap2	minimap2	absolute path of minimap2
mm2_index	all_NCBI_genomes.idx	minimap2	filename of minimap2 index, relative to db.dir
mm2_ref	NA	minimap2	filename of minimap2 reference database (.fasta), if there's no index yet (relative to db.dir)
fastq_pair_pattern	_L001	minimap2	In case of Illumina paired-end seq what is the common pattern for the pairs? Usually it is L0001, and what is after 
fastq_suffix	_001.fastq.gz	minimap2	extension of fastq files
Nsec	20	minimap2	number of secondary alignments to keep in minimap2
nproc	42	both	number of cores to use
minitax.dir	/mnt/e/my.R.packages/minitax	minitax	absolute path of minitax home directory
misc.dir	/mnt/e/my.R.packages/Rlyeh-main	minitax	absolute path of minitax home directory
metadata	metadata.tsv	minitax	tab-separated metadata file, the 'sample' coulmn should be the same as the the bamfiles (without the 'pattern' argumetn, i.e. file suffix)
keep.highest.mapq.aln.only	T	minitax	in case of multi-mapping reads, keep the one with the highest MAPQ value only
crop.na.tax	F	minitax	???
multicore	T	minitax	Run minitax in multicore?
saveRAM	F	minitax	Frequently save intermediate files and load them back in, or be as fast as possible?
mapq.filt	NA	minitax	If other than NA, filter out alignments with MAPQ values lower than this
outputs	bam.sum		
methods	BestAln;SpeciesEstimate	minitax	
pardir	NA	minitax	directory of bamfiles
pattern	.bam	minitax	pattern of bamfiles (for subseting)
CIGAR_points	match_score = 1; mismatch_score = -3; insertion_score = -2; deletion_score = -2; gap_opening_penalty = -2; gap_extension_penalty = -1	minitax	rules for scoring the alignments
best.mapq	T	minitax find the alignment with the highest mapq for each read before analysing CIGAR scores

## Obtaining the database
gunzip NCBI_genome_collection/all_NCBI_genomes_sequence_lengths.zip
minitax/download_ncbi_accessions.sh NCBI_genome_collection/all_NCBI_genomes_sequence_lengths.tsv

## Running the tool:
Mappping part  
minitax/minitax.sh minitax_config_allNCBI.txt

Finding the best taxonomic assignement for each read  
minitax/minitax.complete.R minitax_config_NCBI.txt

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



