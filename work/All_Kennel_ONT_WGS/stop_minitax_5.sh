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

if [[ ! -s "$pid_file" ]]; then
  echo "No pid file for mode '${mode}': $pid_file"
  exit 0
fi

root_pid="$(cat "$pid_file")"
if [[ ! "$root_pid" =~ ^[0-9]+$ ]] || ! kill -0 "$root_pid" 2>/dev/null; then
  echo "No active process for mode '${mode}' pid=${root_pid}"
  echo "stopped" > "$status_file"
  exit 0
fi

collect_descendants() {
  local parent="$1"
  local child
  while read -r child; do
    [[ -n "$child" ]] || continue
    echo "$child"
    collect_descendants "$child"
  done < <(pgrep -P "$parent" 2>/dev/null || true)
}

mapfile -t descendants < <(collect_descendants "$root_pid" | awk '!seen[$0]++')
targets=("${descendants[@]}" "$root_pid")

echo "Stopping minitax ${mode} process tree:"
for pid in "${targets[@]}"; do
  ps -p "$pid" -o pid=,ppid=,cmd= 2>/dev/null || true
done

kill "${targets[@]}" 2>/dev/null || true
sleep 3

remaining=()
for pid in "${targets[@]}"; do
  if kill -0 "$pid" 2>/dev/null; then
    remaining+=("$pid")
  fi
done

if [[ "${#remaining[@]}" -gt 0 ]]; then
  echo "Forcing remaining process(es): ${remaining[*]}"
  kill -KILL "${remaining[@]}" 2>/dev/null || true
fi

echo "stopped" > "$status_file"
echo "Stopped mode '${mode}'."
