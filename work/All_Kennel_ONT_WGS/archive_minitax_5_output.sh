#!/usr/bin/env bash
set -Eeuo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$project_dir"

config="${CONFIG_FILE:-${project_dir}/minitax_bacteria_refseq_5.config.tsv}"
[[ -f "$config" ]] || { echo "Config not found: $config" >&2; exit 1; }

output_name="$(awk -F '\t' '$1 == "outdir" { print $2; exit }' "$config")"
if [[ -z "$output_name" || "$output_name" == "NA" ]]; then
  echo "Config must define a non-empty outdir: $config" >&2
  exit 1
fi

output_target="/mnt/d/data/Gemini/All_Kennel_ONT_WGS/${output_name}"

if [[ ! -d "$output_target" ]]; then
  mkdir -p "$output_target"
  echo "Created empty output target: $output_target"
  exit 0
fi

if ! find "$output_target" -mindepth 1 -print -quit | grep -q .; then
  echo "Output target is already empty: $output_target"
  exit 0
fi

stamp="$(date +%Y%m%d_%H%M%S)"
archive="${output_target}_regular_minimap2_aborted_${stamp}"
mv "$output_target" "$archive"
mkdir -p "$output_target"

echo "Archived previous output:"
echo "  $archive"
echo "Fresh output target:"
echo "  $output_target"
