#!/usr/bin/env bash
# Environment preparation
# -----------------------
# FreeSurfer >= 7.3 with a valid licence file.
# Before running: export FREESURFER_HOME=/path/to/freesurfer
#                 source "$FREESURFER_HOME/SetUpFreeSurfer.sh"
# Optional:       export FS_LICENSE=/path/to/license.txt
#
# Input:  one T1-weighted NIfTI and a de-identified subject label.
# Output: a standard recon-all subject directory under --subjects-dir.

set -euo pipefail

usage() {
  echo "Usage: $0 --t1 T1w.nii.gz --subject SUBJECT --subjects-dir DIR [--threads N]"
}

t1=""
subject=""
subjects_dir=""
threads="4"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --t1) t1="$2"; shift 2 ;;
    --subject) subject="$2"; shift 2 ;;
    --subjects-dir) subjects_dir="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$t1" && -n "$subject" && -n "$subjects_dir" ]] || { usage >&2; exit 2; }
[[ -f "$t1" ]] || { echo "Missing T1 input: $t1" >&2; exit 1; }
command -v recon-all >/dev/null || { echo "recon-all is not available; source SetUpFreeSurfer.sh" >&2; exit 1; }
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "--threads must be a positive integer" >&2; exit 2; }

mkdir -p "$subjects_dir"
export SUBJECTS_DIR="$subjects_dir"
export OMP_NUM_THREADS="$threads"
recon-all -subject "$subject" -i "$t1" -all -parallel -openmp "$threads"
