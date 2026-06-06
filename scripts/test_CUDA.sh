#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'EOF'
Detect whether this system should use Parabricks minimap2 or mm2-fast.

Usage:
  scripts/test_CUDA.sh [options]

Options:
  --image IMAGE    Parabricks Docker image.
                   Default: nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1
  --pull           Pull the Parabricks image before testing.
  --self-test      Run a tiny pbrun minimap2 map-ont job in a temp directory.
  --quiet          Suppress diagnostics; exit 0 only when Parabricks is usable.
  -h, --help       Show this help.

The script prints recommended_mapper=parabricks when NVIDIA GPU + Docker +
Parabricks checks pass. Otherwise it prints recommended_mapper=mm2-fast.
It never writes repo-tracked files.
EOF
}

IMAGE="nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1"
PULL=0
SELF_TEST=0
QUIET=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --image) IMAGE="${2:-}"; shift 2 ;;
    --pull) PULL=1; shift ;;
    --self-test) SELF_TEST=1; shift ;;
    --quiet) QUIET=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

log() {
  if [[ "$QUIET" -eq 0 ]]; then
    echo "$*"
  fi
}

recommend_mm2_fast() {
  local reason="$1"
  log "status=failed"
  log "reason=$reason"
  log "recommended_mapper=mm2-fast"
  exit 1
}

command -v nvidia-smi >/dev/null 2>&1 || recommend_mm2_fast "nvidia-smi not found"
nvidia-smi >/dev/null 2>&1 || recommend_mm2_fast "nvidia-smi failed"
log "nvidia_smi=ok"

command -v docker >/dev/null 2>&1 || recommend_mm2_fast "docker not found"
docker info >/dev/null 2>&1 || recommend_mm2_fast "docker daemon is not reachable"
log "docker=ok"

if [[ "$PULL" -eq 1 ]]; then
  log "pulling_image=$IMAGE"
  docker pull "$IMAGE" >/dev/null || recommend_mm2_fast "failed to pull Parabricks image: $IMAGE"
fi

docker image inspect "$IMAGE" >/dev/null 2>&1 || recommend_mm2_fast "Parabricks image is not present locally; rerun with --pull"
log "parabricks_image=ok"

docker run --rm --gpus all "$IMAGE" pbrun minimap2 --version >/dev/null 2>&1 || \
  recommend_mm2_fast "Docker GPU runtime or Parabricks minimap2 failed"
log "docker_gpu_runtime=ok"

if [[ "$SELF_TEST" -eq 1 ]]; then
  tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/minitax-cuda-test.XXXXXX")"
  cleanup() {
    rm -rf "$tmpdir"
  }
  trap cleanup EXIT

  cat > "${tmpdir}/ref.fa" <<'EOF'
>tiny_ref
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOF
  cat > "${tmpdir}/reads.fastq" <<'EOF'
@tiny_read
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

  docker run --rm --gpus all \
    --volume "${tmpdir}:/workdir" \
    --workdir /workdir \
    "$IMAGE" \
    pbrun minimap2 \
    --ref /workdir/ref.fa \
    --in-fq /workdir/reads.fastq \
    --out-bam /workdir/test.bam \
    --preset map-ont \
    --num-gpus 1 \
    --num-threads 1 \
    --tmp-dir /workdir >/dev/null 2>&1 || \
      recommend_mm2_fast "Parabricks minimap2 self-test failed"

  [[ -s "${tmpdir}/test.bam" ]] || recommend_mm2_fast "Parabricks self-test did not produce a BAM"
  log "parabricks_self_test=ok"
fi

log "status=ok"
log "recommended_mapper=parabricks"
exit 0
