#!/bin/bash

# Check if a first argument is provided
if [ -z "$1" ]; then
    echo "Please provide the path to the configuration file as the first argument."
    exit 1
fi

# If provided, use the first argument as the config file path
CONFIG_FILE="$1"


# Load the configuration file into an associative array in bash
declare -A config
while IFS=$'\t' read -r argument value step description; do
    config[$argument]=$value
done < $CONFIG_FILE

## Parameters
Nsec="${config[Nsec]}"
mm2="${config[mm2_path]}"
NPROC="${config[nproc]}"

## Reference index path
index="${config[db.dir]}/${config[mm2_index]}"
## TO-DO put in a step that will make an index if the reference is a fasta using the mm2_ref argument


# Create the output directory based on config values
outdir=minitax_"${config[project]}_${config[platform]}_${config[Vregion]}_${config[db]}"
if [ "${config[outdir]}" == "NA" ]; then
    mkdir -p "${outdir}"
fi
bamoutdir="${outdir}"/bam
mkdir "$bamoutdir"

## Input directory
indir="${config[indir]}" #
## fastq file suffix
suffix="${config[fastq_suffix]}"

echo 'The following files will be mapped: '
ls -a "$indir"/*"$suffix"

## Run minimap2
if [ "${config[platform]}" == "Illumina" ]; then
	# Determine if the reads are paired or merged
	if [ "${config[reads]}" == "paired" ]; then
		
		suffix="${config[fastq_suffix]}"
		pattern="${config[fastq_pair_pattern]}"
		if [ $pattern == "NA" ]; then
			echo pattern will not be used
			pattern=""
		fi
		
		if [ "${config[debug]}" == "TRUE" ]; then
			# List all the unique base names of the samples
			echo "This was the file pattern: " "$pattern"
			echo "Constructed grep pattern: ${pattern}_R[1-2]${suffix}"
			echo "Files matching pattern:"
			ls "$indir" | grep -E "${pattern}_R[1-2]${suffix}"
			
			echo "Base names after sed operation:"
			ls "$indir" | grep -E "${pattern}_R[1-2]${suffix}" | sed "s|${pattern}_R[1-2]${suffix}||"
			echo "Unique base names:"
			ls "$indir" | grep -E "${pattern}_R[1-2]${suffix}" | sed "s|${pattern}_R[1-2]${suffix}||" | uniq
		fi
	
		echo 'starting minimap2 using ' "$NPROC" 'cores...'
		#for filepath in "$indir"/"${pattern}_R"[1-2]"$suffix"; do
			#base=$(basename "$filepath" "_R"[1-2]"$suffix")
			#echo "$filepath"
			#echo "$base"
		for base in $(ls $indir | grep -E "${pattern}_R[1-2]${suffix}" | sed "s|${pattern}_R[1-2]${suffix}||" | uniq);do
			R1="${base}${pattern}_R1${suffix}"
			R2="${base}${pattern}_R2${suffix}"
			
			if [ "${config[debug]}" == "TRUE" ]; then
				echo "Constructed R1 path: ${indir}/${R1}"
				echo "Constructed R2 path: ${indir}/${R2}"
			fi
			
			echo 'started mapping ' "$R1" 'and ' "$R2" ' to ' "$index" ' at:'
			date
			$mm2 -ax sr -t "$NPROC" -Y -C5 --split-prefix -un -N $Nsec ${index} ${indir}/${R1} ${indir}/${R2} \
			| samtools view -b -@ ${NPROC} - \
			| samtools sort - -o ${bamoutdir}/${base}.bam -@ ${NPROC}
			samtools index -@ ${NPROC} ${bamoutdir}/${base}.bam
			echo 'mapping done, on:'
			date
			
		done
	elif [ "${config[reads]}" == "merged" ]; then
		suffix="${config[fastq_suffix]}"
		
		echo 'starting minimap2 step for merged reads...'
		for mergedRead in $(ls $indir/*${suffix});do
			base=$(basename $mergedRead $suffix)
			
			if [ "${config[debug]}" == "TRUE" ]; then
				echo "Constructed merged read path: ${mergedRead}"
			fi
			echo 'started mapping ' "$mergedRead" ' to ' "$index" ' at:'
			date
			$mm2 -ax sr -t "$NPROC" -Y -C5 --split-prefix -un -N $Nsec ${index} ${mergedRead} \
			| samtools view -b -@ ${NPROC} - \
			| samtools sort - -o ${bamoutdir}/${base}.bam -@ ${NPROC}
			samtools index -@ ${NPROC} ${bamoutdir}/${base}.bam
			echo 'mapping done, on:'
			date
		done
	else
		echo "Unknown read type for Illumina platform in configuration. Please specify 'paired' or 'merged'."
		exit 1
	fi
fi

if [ "${config[platform]}" == "ONT" ]; then
	suffix="${config[fastq_suffix]}"
	echo 'started mapping with ONT settings...'
	for fastq in $(ls $indir/*$suffix); do
		base=$(basename $indir/$fastq $suffix)
		echo 'started mapping ' "$base" ' to ' "$index" ' at:'
		date
		$mm2 -ax map-ont -t "$NPROC" -Y -C5 --split-prefix -un -N $Nsec ${index} ${fastq} \
		| samtools view -b -@ ${NPROC} - \
		| samtools sort - -o ${bamoutdir}/${base}.bam -@ ${NPROC}
		samtools index -@ ${NPROC} ${bamoutdir}/${base}.bam 
		echo 'mapping done, on:'
		date
done
fi


if [ "${config[platform]}" == "PacBio" ]; then
	suffix="${config[fastq_suffix]}"
	echo 'started mapping with PacBio settings...'
	for fastq in $(ls $indir/*$suffix); do
		base=$(basename $indir/$fastq $suffix)
		echo 'started mapping ' "$base" ' to ' "$index" ' at:'
		date
		$mm2 -ax map-pb -t "$NPROC" -Y -C5 --split-prefix -un -N $Nsec ${index} ${fastq} \
		| samtools view -b -@ ${NPROC} - \
		| samtools sort - -o ${bamoutdir}/${base}.bam -@ ${NPROC}
		samtools index -@ ${NPROC} ${bamoutdir}/${base}.bam 
		echo 'mapping done, on:'
		date
done
fi

#Rscript minitax.R minitax_config.txt


