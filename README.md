# 𝑚ＩＮＩＴΛΧ

Versatile taxonomic assignment tool for metagenomic and amplicon reads using minimap2/mm2-fast plus post-alignment taxonomic refinement in R.

Minitax supports ONT, PacBio, and Illumina reads from metagenomic WGS and marker-gene workflows such as full-length 16S rRNA sequencing. Reads are aligned to a user-selected reference database, then minitax imports the alignments, applies optional MAPQ and CIGAR-score refinement, and summarizes taxonomic assignments.

## Main workflow

1. Build or download a compatible database bundle.
2. Prepare a tab-separated `minitax_config.txt`.
3. Validate the config.
4. Run mapping with `minitax`.
5. Run classification/summarisation with `minitax-complete`.

```bash
minitax-validate minitax_config.txt
minitax minitax_config.txt
minitax-complete minitax_config.txt
```

## Supported database types

| `db` value | Intended use | Required database files |
|---|---|---|
| `all_NCBI_genomes` | NCBI genome collection for WGS reads | `NCBI.db.tsv`, `NCBI.db.uni.tsv`, `NCBI.db.genomesize.tsv`, minimap2 index |
| `ncbi_refseq_16S` | NCBI RefSeq TargetedLoci 16S rRNA references | `NCBI.db.tsv`, `NCBI.db.uni.tsv`, `db_data.tsv`, `ncbi_refseq_16s.fna.gz`, minimap2 index |
| `GTDB_SSU` | GTDB representative SSU/16S references | `gtdb_ssu_taxonomy.tsv`, FASTA, minimap2 index |
| `EMUdb` | EMU-compatible 16S database | `taxonomy.tsv`, `species_taxid.fasta` |
| generic/custom | Any FASTA with matching taxonomy table | `db_data.tsv`, FASTA, minimap2 index |

For custom sequence-name joined databases, `db_data.tsv` should contain at least:

```text
seqnames	taxid	superkingdom	phylum	class	order	family	genus	species
```

`seqnames` must match the first whitespace-delimited token in the FASTA headers and the BAM reference names.

## Installing

Clone the repository and create the recommended conda/mamba environment:

```bash
git clone https://github.com/Balays/minitax.git
cd minitax

mamba env create -f envs/minitax-mm2fast.yml
mamba activate minitax-mm2fast

chmod +x minitax.sh minitax.complete.R minitax_validate.sh scripts/*.sh build_*.sh
```

The environment includes the R/Bioconductor dependencies, minimap2, samtools, and `taxonkit`. To update an existing environment after pulling changes:

```bash
mamba env update -f envs/minitax-mm2fast.yml --prune
mamba activate minitax-mm2fast
```

### Install mm2-fast into the active environment

```bash
scripts/install_mm2_fast_into_env.sh --env-prefix "$CONDA_PREFIX"
mm2-fast --version
```

`mm2-fast` is used as a minimap2-compatible CPU mapper. Vanilla `minimap2` can also be used by setting `mm2_path` accordingly.

## Making minitax available from PATH

The recommended approach is to create symlinks inside the active conda/mamba environment. Run this from the cloned repository after activating the environment:

```bash
mamba activate minitax-mm2fast

chmod +x minitax.sh minitax.complete.R minitax_validate.sh scripts/*.sh build_*.sh

ln -sfn "$(pwd)/minitax.sh" "$CONDA_PREFIX/bin/minitax"
ln -sfn "$(pwd)/minitax.complete.R" "$CONDA_PREFIX/bin/minitax-complete"
ln -sfn "$(pwd)/minitax_validate.sh" "$CONDA_PREFIX/bin/minitax-validate"

hash -r
```

Check:

```bash
which minitax
which minitax-complete
which minitax-validate
```

After this, activating the environment is enough to run `minitax`, `minitax-complete`, and `minitax-validate` from any working directory.

## Configuration file

The configuration file is tab-separated. The first two columns must be `argument` and `value`; optional `step` and `description` columns are allowed.

Minimal example for ONT full-length 16S with NCBI RefSeq TargetedLoci 16S:

```text
argument	value	step	description
platform	ONT	both	Sequencing platform
Vregion	V1V9	both	Full-length 16S region label
db	ncbi_refseq_16S	both	NCBI RefSeq TargetedLoci 16S database
db.dir	/mnt/d/data/databases/RefSeq_16S/ncbi_targeted_loci_16s	both	Database directory
project	MAGOR16S	both	Project label
outdir	minitax_MAGOR16S_ONT_V1V9_ncbi_refseq_16S	both	Output directory
indir	/mnt/d/data/HUNOR/MAGOR_16S/fastq	map	FASTQ input directory
fastq_suffix	.fastq.gz	map	FASTQ suffix
nproc	20	both	Worker/thread count
dt.threads	1	classify	data.table threads per worker
mm2_path	mm2-fast	map	Mapper command
mm2_ref	ncbi_refseq_16s.fna.gz	map	Reference FASTA inside db.dir
mm2_index	ncbi_refseq_16s.idx	map	minimap2 index inside db.dir
mm2_index_batch	128G	map	Avoid multi-part minimap2 indexes
mapper_backend	auto	map	Use supported backend automatically
Nsec	0	map	Seconds to sleep between samples
mapq.filt	NA	classify	Do not remove MAPQ values by default
best.mapq	TRUE	classify	Keep highest MAPQ alignments per read
keep.max.cigar	TRUE	classify	Refine tied hits by CIGAR score
methods	BestAln;SpeciesEstimate;LCA	classify	Classification methods
outputs	bam.sum;best_alignments_w_taxa;sum_taxa	classify	Output groups
```

## Building the NCBI RefSeq TargetedLoci 16S database

Use this database for full-length or near-full-length 16S rRNA reads when you want an NCBI/RefSeq targeted-loci reference rather than GTDB SSU.

The top-level wrapper is:

```text
build_NCBI_refseq_16S.sh
```

The actual implementation lives in:

```text
scripts/build_ncbi_refseq_16s_db.sh
```

This follows the same repository convention as `build_NCBI_genome_collection.sh`, which wraps `scripts/build_ncbi_minitax_db.sh`.

### Build command

```bash
mamba activate minitax-mm2fast
cd /mnt/c/GitHub/minitax

git pull
chmod +x build_NCBI_refseq_16S.sh scripts/build_ncbi_refseq_16s_db.sh

./build_NCBI_refseq_16S.sh \
  --outdir /mnt/d/data/databases/RefSeq_16S/ncbi_targeted_loci_16s \
  --min-length 1200 \
  --threads 20 \
  --minimap2 mm2-fast \
  --force
```

`--min-length 1200` is recommended for full-length ONT V1-V9 analyses to remove short partial 16S records. Use `--min-length 0` to retain all TargetedLoci records.

The builder downloads/uses RefSeq TargetedLoci Bacteria and Archaea FASTA/GBFF files plus the NCBI taxdump. It parses TaxIDs from the GBFF records, expands them to standard ranks with `taxonkit`, normalizes FASTA headers to accession-only identifiers, writes minitax-compatible taxonomy tables, and builds the minimap2 index.

### Output files

```text
ncbi_refseq_16s.fna.gz
ncbi_refseq_16s.idx
NCBI.db.tsv
NCBI.db.uni.tsv
NCBI.db.uni.spec.tsv
db_data.tsv
NCBI_16S_seqlengths.txt
metadata/ncbi_refseq_16s_accession_taxid.tsv
metadata/ncbi_refseq_16s.taxid_lineages.tsv
build_manifest.tsv
```

### Config values

```text
db	ncbi_refseq_16S
db.dir	/mnt/d/data/databases/RefSeq_16S/ncbi_targeted_loci_16s
mm2_ref	ncbi_refseq_16s.fna.gz
mm2_index	ncbi_refseq_16s.idx
mm2_index_batch	128G
```

### Sanity checks

```bash
zcat /mnt/d/data/databases/RefSeq_16S/ncbi_targeted_loci_16s/ncbi_refseq_16s.fna.gz | grep '^>' | head
head /mnt/d/data/databases/RefSeq_16S/ncbi_targeted_loci_16s/NCBI.db.uni.tsv
head /mnt/d/data/databases/RefSeq_16S/ncbi_targeted_loci_16s/metadata/ncbi_refseq_16s_accession_taxid.tsv
```

The FASTA and taxonomy table identifiers should match, for example:

```text
>NR_201932.1
NR_201932.1	...
```

## Building the NCBI genome collection database

Use this for metagenomic WGS reads mapped against RefSeq/GenBank genome assemblies.

Check selected assemblies without downloading genomes:

```bash
./build_NCBI_genome_collection.sh \
  --outdir /mnt/d/data/databases/all_NCBI_genomes \
  --dry-run \
  --max-genomes 10
```

Build the default RefSeq representative/reference collection:

```bash
./build_NCBI_genome_collection.sh \
  --outdir /mnt/d/data/databases/all_NCBI_genomes \
  --threads 32 \
  --index-batch-size 128G
```

Config values:

```text
db	all_NCBI_genomes
db.dir	/mnt/d/data/databases/all_NCBI_genomes
mm2_ref	all_NCBI_genomes.fna.gz
mm2_index	all_NCBI_genomes.idx
mm2_index_batch	128G
```

## Building and using a GTDB SSU database

`GTDB_SSU` is intended for full-length or near-full-length 16S rRNA reads when taxonomy should be consistent with GTDB-Tk and genome/MAG analyses.

A runnable GTDB SSU directory should contain:

```text
gtdb_ssu_reps.fna
gtdb_ssu_reps.idx
gtdb_ssu_taxonomy.tsv
```

`gtdb_ssu_taxonomy.tsv` is a two-column, no-header table:

```text
RS_GCF_031457235.1	d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas_E;s__Pseudomonas_E example
```

For full-length V1-V9 reads, remove very short SSU references before indexing:

```bash
seqkit seq -m 1200 gtdb_ssu_reps.raw.fna > gtdb_ssu_reps.fna
minimap2 -x map-ont -I 128G -d gtdb_ssu_reps.idx gtdb_ssu_reps.fna
```

Config values:

```text
db	GTDB_SSU
db.dir	/mnt/d/data/databases/GTDB_SSU
mm2_ref	gtdb_ssu_reps.fna
mm2_index	gtdb_ssu_reps.idx
mm2_index_batch	128G
```

## Reliable MAPQ with large minimap2 indexes

For large references, minimap2 can split the reference into multiple index parts. Mapping still works, but MAPQ is not reliable across multi-part indexes because each part is scored independently. Minitax uses MAPQ when selecting best alignments, so database indexes should be built with an `-I` value larger than the reference collection.

The bundled builders default to:

```text
mm2_index_batch	128G
```

The index builder writes a temporary index first and only replaces the final index after a successful build. It rejects builds whose minimap2 log shows more than one `loaded/built the index` block.

## Mapper backend: mm2-fast or Parabricks

By default, `mapper_backend=auto` tries NVIDIA Parabricks minimap2 only when the local CUDA/Docker/Parabricks test passes. Otherwise it uses `mm2_path`, usually `mm2-fast` or vanilla `minimap2`.

Relevant config keys:

```text
mm2_path	mm2-fast
mapper_backend	auto
parabricks_image	nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1
parabricks_num_gpus	1
parabricks_extra_flags	NA
```

Test GPU support:

```bash
scripts/test_CUDA.sh --pull --self-test
```

## Running minitax

Validate everything:

```bash
minitax-validate minitax_config.txt
```

Validate only one stage:

```bash
minitax-validate minitax_config.txt --step database
minitax-validate minitax_config.txt --step map
minitax-validate minitax_config.txt --step classify
```

Run mapping:

```bash
minitax minitax_config.txt
```

Useful mapper options:

```bash
minitax minitax_config.txt --dry-run
minitax minitax_config.txt --force
minitax minitax_config.txt --no-validate
```

Run classification/summarisation:

```bash
minitax-complete minitax_config.txt
```

## Classification methods

| Method | Meaning |
|---|---|
| `BestAln` | Resolves tied alignments by taxonomic support/probability and produces one final taxon assignment per read. |
| `LCA` | Assigns the lowest common ancestor shared by retained alignments. |
| `SpeciesEstimate` | Keeps all retained alignments and splits each read's weight across them; this is not a one-read-one-taxon method. |
| `RandAln` | Selects one retained alignment randomly for unresolved tied reads. |

Use `SpeciesEstimate`, not `SpecEst`, in the config.

## Important outputs

Typical output directories/files include:

```text
outdir/bam/
outdir/bam.sum/
outdir/best_alignments_w_taxa/
outdir/sum_taxa/
outdir/classification_timing.tsv
outdir/metadata.tsv
```

For `BestAln` and `LCA`, `sum_taxa/*_tax.identity.tsv` is the read-level final taxonomy table. It has one row per read/sample and includes the full lineage/rank columns plus `tax.identity` and `tax.identity.level`.

For `SpeciesEstimate`, the key abundance column is `norm_count`, because reads are fractionally assigned across retained alignments.

## Prebuilt database bundles

Do not commit database binaries to git. Publish full runnable bundles containing the FASTA, minimap2 index, taxonomy tables, sequence lengths, manifests, and checksums.

Create split archive parts:

```bash
scripts/package_minitax_db.sh \
  --db-dir /mnt/d/data/databases/minitax_bacteria_refseq \
  --name minitax_bacteria_refseq_YYYYMMDD.minimap2-I128G \
  --outdir /mnt/d/data/databases/minitax_bacteria_refseq/bundles \
  --part-size 45G
```

Recommended hosting:

- Zenodo for DOI, metadata, README, checksums, and small files.
- Hugging Face datasets for oversized archive parts when the full bundle exceeds Zenodo's practical quota.

## R/Bioconductor dependencies

The recommended environment file installs the needed R and Bioconductor packages. Core packages include:

```text
Rsamtools
GenomicAlignments
data.table
dplyr
tidyr
ggplot2
readr
stringi
stringr
future.apply
phyloseq
```
