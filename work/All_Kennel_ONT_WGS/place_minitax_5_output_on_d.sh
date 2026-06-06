#!/usr/bin/env bash
set -Eeuo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$project_dir"

config="${1:-${project_dir}/minitax_bacteria_refseq_5.config.tsv}"
[[ -f "$config" ]] || { echo "Config not found: $config" >&2; exit 1; }

output_name="$(awk -F '\t' '$1 == "outdir" { print $2; exit }' "$config")"
if [[ -z "$output_name" || "$output_name" == "NA" ]]; then
  echo "Config must define a non-empty outdir: $config" >&2
  exit 1
fi

if [[ -z "${MINITAX_OUTPUT_ROOT:-}" ]]; then
  mkdir -p "$output_name"
  echo "Output directory:"
  ls -ld "$output_name"
  df -h "$output_name" 2>/dev/null || true
  exit 0
fi

output_target="${MINITAX_OUTPUT_ROOT%/}/${output_name}"

mkdir -p "$output_target"

if [[ -L "$output_name" ]]; then
  current_target="$(readlink -f "$output_name")"
  if [[ "$current_target" != "$output_target" ]]; then
    echo "Existing output symlink points elsewhere: $output_name -> $current_target" >&2
    exit 1
  fi
elif [[ -e "$output_name" ]]; then
  unexpected="$(find "$output_name" -mindepth 1 -maxdepth 1 \
    ! -name bam ! -name logs ! -name tmp ! -name mapping_manifest.tsv -print -quit)"
  if [[ -n "$unexpected" ]]; then
    echo "Refusing to replace non-empty output directory with unexpected file: $unexpected" >&2
    exit 1
  fi
  if find "$output_name"/bam "$output_name"/logs "$output_name"/tmp -mindepth 1 -print -quit 2>/dev/null | grep -q .; then
    echo "Refusing to replace output directory because bam/log/tmp already contain files." >&2
    exit 1
  fi
  rm -rf "$output_name"
  ln -s "$output_target" "$output_name"
else
  ln -s "$output_target" "$output_name"
fi

echo "Output directory:"
ls -ld "$output_name"
df -h "$output_target"
