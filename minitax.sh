#!/usr/bin/env bash
set -uo pipefail

usage() {
  cat <<'EOF'
Map reads for minitax with minimap2 and samtools.

Usage:
  minitax CONFIG [options]
  ./minitax.sh CONFIG [options]

Options:
  --force        Recreate BAM files even when BAM + BAI already exist.
  --dry-run      Validate and print planned mappings without running minimap2.
  --no-validate  Skip minitax_validate.sh --step map.
  -h, --help     Show this help.

The config format is the existing minitax tab-separated config file.
Optional mapper config keys:
  mapper_backend        auto, mm2-fast, minimap2, cpu, or parabricks. Default: auto.
  parabricks_image      Docker image for Parabricks. Default: nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1
  parabricks_num_gpus   GPU count passed to Docker/Parabricks. Default: 1
  parabricks_extra_flags
  mm2_index_batch       minimap2 -I value for emergency index builds. Default: 128G
EOF
}

fail() { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "WARNING: $*" >&2; }
info() { echo "$*" >&2; }

resolve_script_dir() {
  local source="${BASH_SOURCE[0]}"
  while [[ -L "$source" ]]; do
    local dir
    dir="$(cd -P "$(dirname "$source")" >/dev/null 2>&1 && pwd)" || fail "Could not resolve script directory."
    source="$(readlink "$source")" || fail "Could not read symlink: ${BASH_SOURCE[0]}"
    [[ "$source" != /* ]] && source="$dir/$source"
  done
  cd -P "$(dirname "$source")" >/dev/null 2>&1 && pwd
}

is_na() {
  local value="${1:-}"
  [[ -z "$value" || "$value" == "NA" || "$value" == "na" || "$value" == "NULL" || "$value" == "null" ]]
}

strip_quotes() {
  local value="${1:-}"
  value="${value#"${value%%[![:space:]]*}"}"
  value="${value%"${value##*[![:space:]]}"}"
  if [[ "$value" =~ ^\".*\"$ || "$value" =~ ^\'.*\'$ ]]; then
    value="${value:1:${#value}-2}"
  fi
  printf '%s' "$value"
}

normalize_path() {
  local path
  path="$(strip_quotes "${1:-}")"
  if is_na "$path"; then
    printf '%s' ""
    return 0
  fi
  path="${path//\\//}"
  if [[ "$path" =~ ^([A-Za-z]):/(.*)$ ]]; then
    local drive="${BASH_REMATCH[1],,}"
    printf '/mnt/%s/%s' "$drive" "${BASH_REMATCH[2]}"
  else
    printf '%s' "$path"
  fi
}

resolve_command() {
  local value="$1"
  local label="$2"
  value="$(normalize_path "$value")"
  if [[ "$value" == */* ]]; then
    [[ -x "$value" ]] || fail "$label is not executable: $value"
    printf '%s' "$value"
  else
    local found
    found="$(command -v "$value" || true)"
    [[ -n "$found" ]] || fail "$label not found in PATH: $value"
    printf '%s' "$found"
  fi
}

cfg() {
  local key="$1"
  local default="${2:-}"
  local value="${config[$key]:-$default}"
  if is_na "$value"; then
    printf '%s' "$default"
  else
    printf '%s' "$value"
  fi
}

require_cfg() {
  local key="$1"
  local value="${config[$key]:-}"
  if is_na "$value"; then
    fail "Missing required config value: $key"
  fi
}

safe_name() {
  local value="$1"
  value="${value//[^A-Za-z0-9._-]/_}"
  printf '%s' "$value"
}

abs_dir() { (cd "$1" && pwd -P); }
is_under_dir() { [[ "$1" == "$2" || "$1" == "$2"/* ]]; }

parabricks_supports_preset() {
  case "$1" in
    map-ont|map-hifi|map-pbmm2|lr:hq|splice|splice:hq|splice:sr) return 0 ;;
    *) return 1 ;;
  esac
}

CONFIG_FILE=""
FORCE=0
DRY_RUN=0
RUN_VALIDATE=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --force) FORCE=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    --no-validate) RUN_VALIDATE=0; shift ;;
    -h|--help) usage; exit 0 ;;
    --*) fail "Unknown option: $1" ;;
    *)
      [[ -z "$CONFIG_FILE" ]] || fail "Unexpected argument: $1"
      CONFIG_FILE="$1"
      shift
      ;;
  esac
done

[[ -n "$CONFIG_FILE" ]] || { usage; fail "CONFIG is required."; }
[[ -f "$CONFIG_FILE" ]] || fail "Configuration file not found: $CONFIG_FILE"

script_dir="$(resolve_script_dir)"
info "minitax mapper: ${script_dir}/minitax.sh"
info "Config: $CONFIG_FILE"

if [[ "$RUN_VALIDATE" -eq 1 ]]; then
  validator="${script_dir}/minitax_validate.sh"
  if [[ -x "$validator" ]]; then
    echo "Validating mapping config..."
    bash "$validator" "$CONFIG_FILE" --step map || fail "Config validation failed."
    info "Validation passed; starting mapping setup."
  else
    warn "minitax_validate.sh is not executable or not found at ${validator}; continuing without preflight validation."
  fi
else
  info "Skipping preflight validation because --no-validate was requested."
fi

declare -A config
while IFS=$'\t' read -r argument value step description extra || [[ -n "${argument:-}" ]]; do
  argument="${argument//$'\r'/}"
  value="${value//$'\r'/}"
  if [[ "$argument" == "argument" || -z "$argument" || "$argument" == \#* ]]; then
    continue
  fi
  config["$argument"]="$value"
done < "$CONFIG_FILE"

for key in platform db db.dir indir fastq_suffix mm2_index Nsec nproc; do require_cfg "$key"; done

platform="$(cfg platform)"
db="$(strip_quotes "$(cfg db)")"
db_dir="$(normalize_path "$(cfg db.dir)")"
indir="$(normalize_path "$(cfg indir)")"
suffix="$(cfg fastq_suffix)"
samtools_bin="$(resolve_command samtools "samtools")"
nsec="$(cfg Nsec)"
nproc="$(cfg nproc)"
project="$(cfg project project)"
vregion="$(cfg Vregion WGS)"
reads="$(cfg reads paired)"
pair_pattern="$(cfg fastq_pair_pattern "")"
mm2_index="$(cfg mm2_index)"
mm2_index="${mm2_index//\\//}"
mm2_ref="$(normalize_path "$(cfg mm2_ref "")")"
mm2_index_batch="$(cfg mm2_index_batch 128G)"
mapper_backend="$(cfg mapper_backend auto)"
parabricks_image="$(cfg parabricks_image nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1)"
parabricks_num_gpus="$(cfg parabricks_num_gpus 1)"
parabricks_extra_flags="$(cfg parabricks_extra_flags "")"
debug="$(cfg debug FALSE)"
mm2_path_value="$(cfg mm2_path "")"
mm2_bin=""

[[ "$nsec" =~ ^[0-9]+$ ]] || fail "Nsec must be a non-negative integer."
[[ "$nproc" =~ ^[0-9]+$ && "$nproc" -gt 0 ]] || fail "nproc must be a positive integer."
[[ "$mm2_index_batch" =~ ^[0-9]+[KkMmGgTt]?$ ]] || fail "mm2_index_batch must look like 128G, 64000M, or 500000000."
case "$mapper_backend" in
  auto|mm2-fast|minimap2|cpu|parabricks) ;;
  *) fail "mapper_backend must be auto, mm2-fast, minimap2, cpu, or parabricks; got '$mapper_backend'." ;;
esac
if [[ "$mapper_backend" != "parabricks" ]]; then
  require_cfg mm2_path
  mm2_bin="$(resolve_command "$mm2_path_value" "minimap2")"
elif [[ -n "$mm2_path_value" ]]; then
  mm2_bin="$(resolve_command "$mm2_path_value" "minimap2")"
fi

[[ -d "$indir" ]] || fail "FASTQ input directory not found: $indir"
[[ -d "$db_dir" ]] || fail "Database directory not found: $db_dir"

default_outdir="minitax_${project}_${platform}_${vregion}_${db}"
outdir="$(normalize_path "$(cfg outdir "")")"
if is_na "$outdir"; then
  outdir="$default_outdir"
fi

bamoutdir="${outdir}/bam"
logdir="${outdir}/logs"
tmpdir="${outdir}/tmp"
mkdir -p "$bamoutdir" "$logdir" "$tmpdir" || fail "Could not create output directories under: $outdir"

index="${db_dir}/${mm2_index}"
mm2_ref_path=""
mm2_ref_rel=""
if [[ -n "$mm2_ref" ]]; then
  if [[ -s "${db_dir}/${mm2_ref}" ]]; then
    mm2_ref_path="${db_dir}/${mm2_ref}"
    mm2_ref_rel="$mm2_ref"
  elif [[ -s "$mm2_ref" ]]; then
    mm2_ref_path="$mm2_ref"
    if is_under_dir "$mm2_ref_path" "$db_dir"; then
      mm2_ref_rel="${mm2_ref_path#"$db_dir"/}"
    fi
  fi
fi

if [[ ! -s "$index" ]]; then
  if [[ -n "$mm2_ref_path" && -s "$mm2_ref_path" ]]; then
    if [[ "$DRY_RUN" -eq 1 ]]; then
      info "Would build missing minimap2 index: $index from $mm2_ref_path with -I $mm2_index_batch"
    else
      [[ -n "$mm2_bin" ]] || fail "mm2_path is required to build a missing minimap2 index."
      info "Building missing minimap2 index: $index with -I $mm2_index_batch"
      "${script_dir}/scripts/build_minimap2_index.sh" \
        --db-dir "$db_dir" \
        --fasta "$mm2_ref_path" \
        --index "$mm2_index" \
        --threads "$nproc" \
        --index-batch-size "$mm2_index_batch" \
        --minimap2 "$mm2_bin" \
        --force || fail "Failed to build minimap2 index."
    fi
  else
    fail "minimap2 index not found: $index. Set mm2_ref to a FASTA if you want the mapper to build it."
  fi
fi

manifest="${outdir}/mapping_manifest.tsv"
printf 'sample\tplatform\tpreset\tread1\tread2\tbam\tstatus\tlog\tstarted\tended\n' > "$manifest" || fail "Could not write mapping manifest: $manifest"

shopt -s nullglob

declare -a planned_samples=()
declare -a planned_presets=()
declare -a planned_read1=()
declare -a planned_read2=()

add_plan() {
  local sample="$1" preset="$2" read1="$3" read2="${4:-}"
  planned_samples+=("$sample")
  planned_presets+=("$preset")
  planned_read1+=("$read1")
  planned_read2+=("$read2")
}

plan_single_reads() {
  local preset="$1"
  local files=("$indir"/*"$suffix")
  [[ "${#files[@]}" -gt 0 ]] || fail "No FASTQ files ending with '$suffix' found in $indir"
  for fastq in "${files[@]}"; do
    local filename sample
    filename="$(basename "$fastq")"
    sample="${filename%"$suffix"}"
    add_plan "$sample" "$preset" "$fastq" ""
  done
}

plan_illumina_paired() {
  if is_na "$pair_pattern"; then
    pair_pattern=""
  fi
  local r1_suffix="${pair_pattern}_R1${suffix}"
  local r2_suffix="${pair_pattern}_R2${suffix}"
  local r1_files=("$indir"/*"$r1_suffix")
  [[ "${#r1_files[@]}" -gt 0 ]] || fail "No Illumina R1 files ending with '$r1_suffix' found in $indir"
  for r1 in "${r1_files[@]}"; do
    local filename sample r2
    filename="$(basename "$r1")"
    sample="${filename%"$r1_suffix"}"
    r2="${indir}/${sample}${r2_suffix}"
    [[ -f "$r2" ]] || fail "Missing R2 pair for sample '$sample': $r2"
    add_plan "$sample" "sr" "$r1" "$r2"
  done
  local r2_files=("$indir"/*"$r2_suffix")
  for r2 in "${r2_files[@]}"; do
    local filename sample r1
    filename="$(basename "$r2")"
    sample="${filename%"$r2_suffix"}"
    r1="${indir}/${sample}${r1_suffix}"
    [[ -f "$r1" ]] || fail "Missing R1 pair for sample '$sample': $r1"
  done
}

case "$platform" in
  Illumina)
    case "$reads" in
      paired) plan_illumina_paired ;;
      merged) plan_single_reads "sr" ;;
      *) fail "For Illumina, reads must be paired or merged; got '$reads'." ;;
    esac
    ;;
  ONT) plan_single_reads "map-ont" ;;
  PacBio) plan_single_reads "map-pb" ;;
  *) fail "platform must be Illumina, ONT, or PacBio; got '$platform'." ;;
esac

info "Output directory: $outdir"
info "BAM directory: $bamoutdir"
info "Log directory: $logdir"
info "Reference index: $index"
info "Samples planned: ${#planned_samples[@]}"

if [[ "$debug" == "TRUE" || "$debug" == "T" || "$debug" == "true" ]]; then
  for i in "${!planned_samples[@]}"; do info "  ${planned_samples[$i]}: ${planned_read1[$i]} ${planned_read2[$i]}"; done
fi

db_dir_abs="$(abs_dir "$db_dir")" || fail "Could not resolve absolute db.dir: $db_dir"
indir_abs="$(abs_dir "$indir")" || fail "Could not resolve absolute input directory: $indir"
tmpdir_abs="$(abs_dir "$tmpdir")" || fail "Could not resolve absolute tmp directory: $tmpdir"

parabricks_available=0
if [[ "$mapper_backend" == "auto" || "$mapper_backend" == "parabricks" ]]; then
  cuda_probe="${script_dir}/scripts/test_CUDA.sh"
  if [[ -f "$cuda_probe" ]] && bash "$cuda_probe" --quiet --image "$parabricks_image" >/dev/null 2>&1; then
    parabricks_available=1
  fi
fi

choose_backend() {
  local preset="$1"
  case "$mapper_backend" in
    mm2-fast|minimap2|cpu) printf '%s' "minimap2" ;;
    parabricks)
      parabricks_supports_preset "$preset" || fail "Parabricks minimap2 does not support preset '$preset'; use mapper_backend=auto or mm2-fast."
      [[ "$parabricks_available" -eq 1 ]] || fail "mapper_backend=parabricks requested, but scripts/test_CUDA.sh did not pass."
      [[ -n "$mm2_ref_path" && -n "$mm2_ref_rel" ]] || fail "Parabricks requires mm2_ref to point to a FASTA inside db.dir."
      printf '%s' "parabricks"
      ;;
    auto)
      if [[ "$parabricks_available" -eq 1 ]] && parabricks_supports_preset "$preset" && [[ -n "$mm2_ref_path" && -n "$mm2_ref_rel" ]]; then
        printf '%s' "parabricks"
      else
        printf '%s' "minimap2"
      fi
      ;;
  esac
}

first_backend="$(choose_backend "${planned_presets[0]}")"
info "Requested mapper backend: $mapper_backend"
info "Active mapper backend for '${planned_presets[0]}': $first_backend"
if [[ "$first_backend" == "parabricks" ]]; then
  info "Parabricks image: $parabricks_image"
  info "Reference FASTA: $mm2_ref_path"
else
  info "minimap2 command: $mm2_bin"
fi

declare -a mapped=()
declare -a skipped=()
declare -a failed=()
declare -a dry_planned=()

map_sample() {
  local sample="$1" preset="$2" read1="$3" read2="${4:-}"
  local safe_sample outbam tmpbam log split_prefix started ended status sample_backend
  safe_sample="$(safe_name "$sample")"
  outbam="${bamoutdir}/${sample}.bam"
  tmpbam="${tmpdir}/${safe_sample}.$$.bam"
  log="${logdir}/${sample}.minimap2.log"
  split_prefix="${tmpdir}/${safe_sample}.split"
  sample_backend="$(choose_backend "$preset")"

  if [[ "$FORCE" -eq 0 && -s "$outbam" && -s "${outbam}.bai" ]]; then
    info "Skipping existing BAM: $outbam"
    printf '%s\t%s\t%s\t%s\t%s\t%s\tskipped\t%s\t%s\t%s\n' "$sample" "$platform" "$preset" "$read1" "$read2" "$outbam" "$log" "$(date -Is)" "$(date -Is)" >> "$manifest"
    return 2
  fi

  started="$(date -Is)"
  if [[ "$DRY_RUN" -eq 1 ]]; then
    info "Planning sample '$sample' with preset '$preset' using backend '$sample_backend'..."
    printf '%s\t%s\t%s\t%s\t%s\t%s\tplanned\t%s\t%s\t%s\n' "$sample" "$platform" "$preset" "$read1" "$read2" "$outbam" "$log" "$started" "$(date -Is)" >> "$manifest"
    return 3
  fi

  info "Mapping sample '$sample' with preset '$preset' using backend '$sample_backend'..."
  info "  log: $log"
  (
    set -uo pipefail
    echo "sample: $sample"
    echo "started: $started"
    echo "backend: $sample_backend"
    echo "preset: $preset"
    echo "index: $index"
    if [[ -n "$mm2_ref_path" ]]; then echo "ref: $mm2_ref_path"; fi
    echo "read1: $read1"
    if [[ -n "$read2" ]]; then echo "read2: $read2"; fi
    echo "bam: $outbam"
    echo

    if [[ "$sample_backend" == "parabricks" ]]; then
      declare -a pb_cmd pb_extra
      pb_cmd=(docker run --rm --gpus "$parabricks_num_gpus" --volume "${db_dir_abs}:/db:ro" --volume "${indir_abs}:/reads:ro" --volume "${tmpdir_abs}:/tmpdir" --workdir /tmpdir "$parabricks_image" pbrun minimap2 --ref "/db/${mm2_ref_rel}" --index "/db/${mm2_index}" --in-fq "/reads/$(basename "$read1")" --out-bam "/tmpdir/$(basename "$tmpbam")" --preset "$preset" --num-threads "$nproc" --num-gpus "$parabricks_num_gpus" --tmp-dir /tmpdir)
      if [[ -n "$read2" ]]; then pb_cmd+=("/reads/$(basename "$read2")"); fi
      if [[ -n "$parabricks_extra_flags" ]]; then read -r -a pb_extra <<< "$parabricks_extra_flags"; pb_cmd+=("${pb_extra[@]}"); fi
      printf 'command:'; printf ' %q' "${pb_cmd[@]}"; printf '\n\n'
      "${pb_cmd[@]}"
    else
      if [[ -n "$read2" ]]; then
        "$mm2_bin" -ax "$preset" -t "$nproc" -Y -C5 --split-prefix "$split_prefix" -N "$nsec" "$index" "$read1" "$read2" |
          "$samtools_bin" view -b -@ "$nproc" - |
          "$samtools_bin" sort -@ "$nproc" -o "$tmpbam" -
      else
        "$mm2_bin" -ax "$preset" -t "$nproc" -Y -C5 --split-prefix "$split_prefix" -N "$nsec" "$index" "$read1" |
          "$samtools_bin" view -b -@ "$nproc" - |
          "$samtools_bin" sort -@ "$nproc" -o "$tmpbam" -
      fi
    fi
    "$samtools_bin" index -@ "$nproc" "$tmpbam"
  ) > "$log" 2>&1
  status=$?
  ended="$(date -Is)"

  if [[ "$status" -eq 0 ]]; then
    mv -f "$tmpbam" "$outbam"
    mv -f "${tmpbam}.bai" "${outbam}.bai"
    printf '%s\t%s\t%s\t%s\t%s\t%s\tmapped\t%s\t%s\t%s\n' "$sample" "$platform" "$preset" "$read1" "$read2" "$outbam" "$log" "$started" "$ended" >> "$manifest"
    info "Finished sample '$sample'."
    return 0
  fi

  rm -f "$tmpbam" "${tmpbam}.bai"
  printf '%s\t%s\t%s\t%s\t%s\t%s\tfailed\t%s\t%s\t%s\n' "$sample" "$platform" "$preset" "$read1" "$read2" "$outbam" "$log" "$started" "$ended" >> "$manifest"
  warn "Mapping failed for sample '$sample'. Last log lines:"
  tail -n 20 "$log" >&2 || true
  return "$status"
}

for i in "${!planned_samples[@]}"; do
  status=0
  map_sample "${planned_samples[$i]}" "${planned_presets[$i]}" "${planned_read1[$i]}" "${planned_read2[$i]}" || status=$?
  case "$status" in
    0) mapped+=("${planned_samples[$i]}") ;;
    2) skipped+=("${planned_samples[$i]}") ;;
    3) dry_planned+=("${planned_samples[$i]}") ;;
    *) failed+=("${planned_samples[$i]}") ;;
  esac
done

info
info "Mapping summary"
info "  planned: ${#planned_samples[@]}"
info "  mapped:  ${#mapped[@]}"
info "  skipped: ${#skipped[@]}"
if [[ "$DRY_RUN" -eq 1 ]]; then info "  dry-run: ${#dry_planned[@]}"; fi
info "  failed:  ${#failed[@]}"
info "  manifest: $manifest"

if [[ "${#failed[@]}" -gt 0 ]]; then
  info "Failed samples: ${failed[*]}"
  exit 1
fi

if [[ "$DRY_RUN" -eq 1 ]]; then
  info "Dry run complete; no BAM files were written."
fi
