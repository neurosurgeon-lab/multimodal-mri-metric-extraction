#!/usr/bin/env bash
# Environment preparation
# -----------------------
# Docker Desktop/Engine with access to nipreps/fmriprep:25.2.5.
# A valid FreeSurfer licence file is required and is mounted read-only.
# Input must be a valid BIDS dataset containing task-rest BOLD data.
# Output is the standard fMRIPrep derivative in MNI152NLin2009cAsym:res-2.

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: run_fmriprep_docker.sh --bids-dir DIR --output-dir DIR --work-dir DIR
       --participant-label LABEL --fs-license FILE
       [--templateflow-dir DIR] [--image nipreps/fmriprep:25.2.5]
       [--nprocs 4] [--omp-nthreads 2] [--memory-mb 9500]
USAGE
}

bids_dir=""; output_dir=""; work_dir=""; participant=""; fs_license=""
templateflow_dir=""; image="nipreps/fmriprep:25.2.5"
nprocs="4"; omp_threads="2"; memory_mb="9500"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bids-dir) bids_dir="$2"; shift 2 ;;
    --output-dir) output_dir="$2"; shift 2 ;;
    --work-dir) work_dir="$2"; shift 2 ;;
    --participant-label) participant="${2#sub-}"; shift 2 ;;
    --fs-license) fs_license="$2"; shift 2 ;;
    --templateflow-dir) templateflow_dir="$2"; shift 2 ;;
    --image) image="$2"; shift 2 ;;
    --nprocs) nprocs="$2"; shift 2 ;;
    --omp-nthreads) omp_threads="$2"; shift 2 ;;
    --memory-mb) memory_mb="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$bids_dir" && -n "$output_dir" && -n "$work_dir" && -n "$participant" && -n "$fs_license" ]] || { usage >&2; exit 2; }
[[ -d "$bids_dir" ]] || { echo "Missing BIDS directory: $bids_dir" >&2; exit 1; }
[[ -f "$fs_license" ]] || { echo "Missing FreeSurfer licence: $fs_license" >&2; exit 1; }
command -v docker >/dev/null || { echo "docker is not available" >&2; exit 1; }

mkdir -p "$output_dir" "$work_dir"
docker_args=(
  run --rm
  -v "$bids_dir:/data:ro"
  -v "$output_dir:/out"
  -v "$work_dir:/work"
  -v "$fs_license:/opt/freesurfer/license.txt:ro"
)
if [[ -n "$templateflow_dir" ]]; then
  mkdir -p "$templateflow_dir"
  docker_args+=(-e TEMPLATEFLOW_HOME=/templateflow -v "$templateflow_dir:/templateflow")
fi

docker "${docker_args[@]}" "$image" /data /out participant \
  --participant-label "$participant" \
  --task-id rest \
  --work-dir /work \
  --output-spaces MNI152NLin2009cAsym:res-2 \
  --fs-no-reconall \
  --use-syn-sdc warn \
  --low-mem \
  --nprocs "$nprocs" \
  --omp-nthreads "$omp_threads" \
  --mem "$memory_mb" \
  --fd-spike-threshold 0.5 \
  --dvars-spike-threshold 1.5 \
  --stop-on-first-crash \
  --notrack
