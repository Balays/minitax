
# đť‘šďĽ©ďĽ®ďĽ©ďĽ´Î›Î§
Versatile Taxonomic Assignment Tool for Metagenomic Reads Using minimap2

Minitax is a taxonomic assignment tool designed for robust profiling across diverse sequencing platforms, including Oxford Nanopore (ONT), PacBio, and Illumina, as well as different library types like metagenomic whole genome sequencing (mWGS) and 16S rRNA gene sequencing. It utilizes minimap2 for initial alignment, followed by sophisticated post-alignment processing to ensure accurate taxonomic assignments.

## How minitax Works
### alignment with minimap2
Minitax begins by aligning reads using minimap2 with platform-specific parameter settings, ensuring the best performance for each sequencing platform.
The default parameters settings are as follows:
```
Platform	Match Score	Mismatch Score	Insertion Score	Deletion Score	Gap Opening Penalty	Gap Extension Penalty	Description
Illumina	2	-4	-3	-3	-4	-2	*Optimized for high-accuracy, short reads. Higher penalties for mismatches and indels to reflect the platformâ€™s low error rate.*
ONT	1	-3	-2	-2	-2	-1	*Adjusted for longer reads with higher error rates. More lenient penalties to accommodate frequent indels and mismatches.*
PacBio	2	-3	-3	-3	-3	-2	*Balanced settings for long, high-fidelity reads (e.g., HiFi mode). Moderate penalties for indels to support accurate alignment in repetitive regions.*
```
### Databases
Minitax supports a variety of databases, including a comprehensive genome collection from NCBI (approximately 16,000 genomes) for WGS reads,
and EMUdb  (https://gitlab.com/treangenlab/emu) or SILVA (https://www.arb-silva.de/) for 16S gene sequencing data.

### Post-alignment processing
The alignment data is imported into R using Rsamtools and merged with database information using data.table for computational efficiency, a design that supports the processing of large datasets.
The alignments may be optionally filtered for lower MAPQ scores (eg. 0-59) in order to decrease the number of false positive hits. From th filtered alignemnts, the software then determines the best for each read based on MAPQ values and the above mentioned user-controllable, platform-optimized CIGAR scoring schemes. For reads with multiple high-confidence alignments (with the same MAPQ and CIGAR scores), minitax provides four refinement methods:
1. Random Alignment (RandAln): Selects a random alignment for reads with multiple alignments that have identical MAPQ and CIGAR scores (faster).
2. Best Alignment (BestAln): Selects the alignment with the best MAPQ and CIGAR score for the most precise taxonomic assignment.
3. Species Estimation (SpeciesEstimate): Uses all alignments to estimate species-level abundances by normalizing the read counts. It will use every alignment of reads (that were kept after MAPQ filtering) and normalize the counts based on the number alignments.

### Summarization and outputs
After determining the best taxonomic assignment for each read, it summarizes read counts at the chosen taxonomic rank (e.g., species), providing outputs in both .tsv and .rds formats (the latter as a phyloseq object for downstream analysis).


## Installing
the program can be downloaded from github using
```
git clone https://github.com/Balays/minitax.git
```


## Configuration file
The configuration file is tab-separated file, and should contain the follwing information and should look like this:
```
argument	value	step	description
platform	ONT	both	Either: 'Illumina', 'PacBio' or 'ONT'
db	'all_NCBI_genomes'	both	options: 'all_NCBI_genomes' or 'EMUdb'
db.dir	/mnt/d/data/databases/all_NCBI_genomes	both	absolute path of database home directory
project	project	optional	project identifier
...
```
*A sample configuration file is provided*

## Building the NCBI genome collection database
Minitax needs more than a FASTA file. The NCBI database directory must contain
the genome FASTA, minimap2 index, taxonomy tables, sequence lengths, and genome
sizes. Use the bundled builder instead of the old accession-only downloader.

Activate an environment with `Rscript`, `minimap2`, `curl` or `wget`, and the R
package `data.table`. `pigz` is optional but speeds up compression. For example:
```
mamba activate minitax
```

Check the selected assemblies without downloading genomes:
```
./build_NCBI_genome_collection.sh \
  --outdir /mnt/d/data/databases/all_NCBI_genomes \
  --dry-run \
  --max-genomes 10
```
The dry run prints the selected assembly count and total genome size in bp/Gbp,
then writes the selected table to `metadata/assembly_summary.filtered.tsv`.

Build a small test database:
```
./build_NCBI_genome_collection.sh \
  --outdir /mnt/d/data/databases/minitax_NCBI_test \
  --max-genomes 5 \
  --threads 8
```

Build the default RefSeq representative/reference collection:
```
./build_NCBI_genome_collection.sh \
  --outdir /mnt/d/data/databases/all_NCBI_genomes \
  --threads 32
```

The output directory will contain:
```
all_NCBI_genomes.fna.gz
all_NCBI_genomes.idx
NCBI.db.tsv
NCBI.db.uni.tsv
NCBI.db.uni.spec.tsv
NCBI.db.genomesize.tsv
NCBI_genome_collection_seqlengths.txt
```

Then set these config values:
```
db	all_NCBI_genomes
db.dir	/mnt/d/data/databases/all_NCBI_genomes
mm2_index	all_NCBI_genomes.idx
```

## Running minitax
Validate the config before launching a run:
```
./minitax_validate.sh minitax_config.txt
```

You can validate only one step when needed:
```
./minitax_validate.sh minitax_config.txt --step map
./minitax_validate.sh minitax_config.txt --step classify
./minitax_validate.sh minitax_config.txt --step database
```

Mappping part
```
minitax/minitax.sh minitax_config_allNCBI.txt
```

Finding the best taxonomic assignement for each read
```
minitax/minitax.complete.R minitax_config_NCBI.txt
```

## Outputs
The outputs include a .tsv file containing the counts for each sample
And an .rds file containing a phyloseq-object

## R-pcakage Dependencies
Rsamtools
readr
ggplot2
tidyr
dplyr
GenomicAlignments
data.table
stringi
stringr
phyloseq
future.apply
