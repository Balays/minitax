#!/usr/bin/env bash
set -Eeuo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$project_dir"

config="${CONFIG_FILE:-${project_dir}/minitax_bacteria_refseq_5.config.tsv}"
[[ -f "$config" ]] || { echo "Config not found: $config" >&2; exit 1; }

mode="${1:-all}"
case "$mode" in
  all) run_arg="--all" ;;
  validate) run_arg="--validate-only" ;;
  dry-run) run_arg="--dry-run" ;;
  map) run_arg="--map-only" ;;
  classify) run_arg="--classify-only" ;;
  *)
    echo "Usage: $0 [all|validate|dry-run|map|classify]" >&2
    exit 2
    ;;
esac

config_indir="$(awk -F '\t' '$1 == "indir" { print $2; exit }' "$config")"
if [[ "$config_indir" == "fastq_5" || "$config_indir" == "./fastq_5" || "$config_indir" == "${project_dir}/fastq_5" ]]; then
  "${project_dir}/prepare_minitax_5_inputs.sh" >/dev/null
fi
"${project_dir}/place_minitax_5_output_on_d.sh" "$config" >/dev/null

outdir_name="$(awk -F '\t' '$1 == "outdir" { print $2; exit }' "$config")"
if [[ -z "$outdir_name" || "$outdir_name" == "NA" ]]; then
  echo "Config must define a non-empty outdir: $config" >&2
  exit 1
fi

outdir="${project_dir}/${outdir_name}"
mkdir -p "$outdir/logs"

pid_file="${outdir}/run_${mode}.pid"
status_file="${outdir}/run_${mode}.status"
log_file="${outdir}/run_${mode}.log"

if [[ -s "$pid_file" ]]; then
  old_pid="$(cat "$pid_file")"
  if [[ "$old_pid" =~ ^[0-9]+$ ]] && kill -0 "$old_pid" 2>/dev/null; then
    echo "Run already active: $mode pid=$old_pid"
    exit 0
  fi
fi

nohup bash -c '
  set +e
  project_dir="$1"
  run_arg="$2"
  status_file="$3"
  mode="$4"
  echo "mode: $mode"
  echo "started: $(date -Is)"
  echo "project_dir: $project_dir"
  echo
  "${project_dir}/run_minitax_5.sh" "$run_arg"
  status=$?
  echo
  echo "ended: $(date -Is)"
  echo "status: $status"
  echo "$status" > "$status_file"
  exit "$status"
' bash "$project_dir" "$run_arg" "$status_file" "$mode" > "$log_file" 2>&1 &

pid=$!
echo "$pid" > "$pid_file"
echo "running" > "$status_file"
echo "Started $mode run: pid=$pid"
echo "Log: $log_file"
