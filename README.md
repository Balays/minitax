
# 𝑚ＩＮＩＴΛΧ
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
the program can be downloaded from github and installed using mamba
```
# 1. Clone minitax
git clone https://github.com/Balays/minitax.git

cd minitax

# 2. Create the env where you want it to live
mamba env create \
  -p ~/mambaforge/envs/minitax-mm2fast \
  -f envs/minitax-mm2fast.yml

# 3. Activate it
mamba activate /home/mdbio/mamba/envs/minitax-mm2fast

# 3. Install zlib for mm2-fast
mamba install -c conda-forge zlib libzlib

# 4. Make scripts executable, useful after Windows/WSL checkouts
chmod +x minitax.sh minitax_validate.sh scripts/*.sh

# 5. Build and install mm2-fast into this env
bash scripts/install_mm2_fast_into_env.sh \
  --env-prefix "$CONDA_PREFIX"
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

### Reliable MAPQ with large minimap2 indexes

For large genome collections, minimap2 can split the reference into multiple
index parts. Mapping still works, but minimap2's MAPQ is not reliable across
multi-part indexes because each part is scored independently. Minitax uses MAPQ
when selecting best alignments, so large database indexes must be built with an
`-I` value larger than the reference collection.

The bundled builders default to `-I 128G`:
```
scripts/build_minimap2_index.sh \
  --db-dir /mnt/d/data/databases/minitax_bacteria_refseq \
  --threads 8 \
  --index-batch-size 128G \
  --force

scripts/build_minimap2_index.sh \
  --db-dir /mnt/d/data/databases/all_NCBI_genomes \
  --threads 8 \
  --index-batch-size 128G \
  --force
```

The index builder writes to `all_NCBI_genomes.idx.tmp` first and only replaces
the final index after a successful build. It rejects builds whose minimap2 log
shows more than one `loaded/built the index` block.

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
  --threads 32 \
  --index-batch-size 128G
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
mm2_ref	all_NCBI_genomes.fna.gz
mapper_backend	auto
mm2_index_batch	128G
```

## Mapper backend: mm2-fast or Parabricks

By default, `mapper_backend=auto` tries NVIDIA Parabricks minimap2 only when the
local CUDA/Docker/Parabricks test passes. Otherwise it uses `mm2_path`, usually
`mm2-fast` or vanilla `minimap2`.

Install mm2-fast after creating/activating the conda environment:
```
mamba env update -f envs/minitax-mm2fast.yml --prune
mamba activate minitax-mm2fast
scripts/install_mm2_fast_into_env.sh
mm2-fast --version
```

The installer clones and builds mm2-fast from
https://github.com/bwa-mem2/mm2-fast, then creates environment-local
`mm2-fast` and `minimap2-mm2fast` commands. The mm2-fast README describes it as
a drop-in accelerated minimap2 implementation for Linux x86_64 CPUs with AVX2
or AVX512, built with `git clone --recursive ... && make`.

Relevant config keys:
```
mm2_path	mm2-fast
mapper_backend	auto
parabricks_image	nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1
parabricks_num_gpus	1
parabricks_extra_flags	NA
```

Test a system before running a GPU mapping job:
```
scripts/test_CUDA.sh --pull --self-test
```

The detector prints either:
```
recommended_mapper=parabricks
```
or:
```
recommended_mapper=mm2-fast
```

Parabricks minimap2 requires Docker with NVIDIA GPU support, a local or pullable
Parabricks image, and `mm2_ref` pointing to the reference FASTA inside `db.dir`.
The current automatic backend uses Parabricks for supported long-read presets
such as `map-ont`; Illumina short-read mapping stays on `mm2_path`.

NVIDIA Parabricks minimap2 reference:
https://docs.nvidia.com/clara/parabricks/latest/documentation/tooldocs/man_minimap2.html

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

Mapping part:
```
./minitax.sh minitax_config.txt
```

The mapper validates the config first, writes BAM files to `outdir/bam`, and
writes per-sample logs to `outdir/logs`. Useful options:
```
./minitax.sh minitax_config.txt --dry-run
./minitax.sh minitax_config.txt --force
./minitax.sh minitax_config.txt --no-validate
```

If `mm2_index` is missing and `mm2_ref` points to a FASTA, the mapper will build
the minimap2 index before mapping, using `mm2_index_batch` as minimap2's `-I`
value.

## Prebuilt database bundles

Do not commit database binaries to git. Publish full runnable bundles containing
the FASTA, minimap2 index, taxonomy tables, sequence lengths, manifests, and
checksums.

Create split archive parts:
```
scripts/package_minitax_db.sh \
  --db-dir /mnt/d/data/databases/minitax_bacteria_refseq \
  --name minitax_bacteria_refseq_YYYYMMDD.minimap2-I128G \
  --outdir /mnt/d/data/databases/minitax_bacteria_refseq/bundles \
  --part-size 45G
```

Recommended hosting:
- Zenodo: DOI, metadata, README, checksums, and small files.
- Hugging Face dataset: oversized archive parts when the full bundle exceeds
  Zenodo's practical quota.

Download and verify a published bundle:
```
scripts/download_minitax_db_bundle.sh \
  --url-base https://huggingface.co/datasets/ORG/DATASET/resolve/main \
  --bundle minitax_bacteria_refseq_YYYYMMDD.minimap2-I128G \
  --outdir /mnt/d/data/databases/minitax_bacteria_refseq
```

Zenodo size-limit reference:
https://support.zenodo.org/help/en-gb/1-upload-deposit/80-what-are-the-size-limitations-of-zenodo

Finding the best taxonomic assignment for each read:
```
./minitax.complete.R minitax_config.txt
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
