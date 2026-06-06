#!/usr/bin/env bash
set -Eeuo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$project_dir"

config="${CONFIG_FILE:-${project_dir}/minitax_bacteria_refseq_5.config.tsv}"
[[ -f "$config" ]] || { echo "Config not found: $config" >&2; exit 1; }

outdir="$(awk -F '\t' '$1 == "outdir" { print $2; exit }' "$config")"
if [[ -z "$outdir" || "$outdir" == "NA" ]]; then
  echo "Config must define a non-empty outdir: $config" >&2
  exit 1
fi

mode="${1:-all}"

pid_file="${outdir}/run_${mode}.pid"
status_file="${outdir}/run_${mode}.status"
log_file="${outdir}/run_${mode}.log"

echo "Project: $project_dir"
echo "Output:  $(readlink -f "$outdir" 2>/dev/null || echo "$outdir")"
echo "Status:  $(cat "$status_file" 2>/dev/null || echo missing)"
echo "PID:     $(cat "$pid_file" 2>/dev/null || echo missing)"
echo

if [[ -s "$pid_file" ]]; then
  pid="$(cat "$pid_file")"
  if [[ "$pid" =~ ^[0-9]+$ ]]; then
    ps -fp "$pid" || true
  fi
fi

echo
echo "Active mapper/classifier processes:"
pgrep -af 'mm2-fast|minimap2|samtools sort|samtools view|minitax.complete|run_minitax_5' || true

echo
echo "Output size:"
du -shL "$outdir" 2>/dev/null || true
df -h "$outdir" 2>/dev/null || true

echo
echo "Manifest tail:"
tail -n 10 "$outdir/mapping_manifest.tsv" 2>/dev/null || true

echo
echo "BAM files:"
find -L "$outdir/bam" -maxdepth 1 -type f -printf '%f\t%s bytes\n' 2>/dev/null | sort || true

echo
echo "Latest minimap2 log tail:"
latest_log="$(find -L "$outdir/logs" -maxdepth 1 -name '*.minimap2.log' -type f -printf '%T@ %p\n' 2>/dev/null | sort -n | tail -n 1 | cut -d' ' -f2-)"
if [[ -n "$latest_log" ]]; then
  echo "$latest_log"
  tail -n 40 "$latest_log"
fi

echo
echo "Run log tail:"
tail -n 40 "$log_file" 2>/dev/null || true
