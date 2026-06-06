#!/usr/bin/env bash
set -Eeuo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$project_dir"

mkdir -p fastq_5

for barcode in 01 02 03 04 05; do
  sample="A_barcode${barcode}.fastq.gz"
  source_path="../fastq/${sample}"
  target_path="fastq_5/${sample}"
  if [[ ! -e "fastq/${sample}" ]]; then
    echo "Missing linked FASTQ: fastq/${sample}" >&2
    exit 1
  fi
  ln -sfn "$source_path" "$target_path"
done

echo "Prepared FASTQ subset:"
ls -lhL fastq_5/*.fastq.gz
