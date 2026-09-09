#!/usr/bin/env bash
# Environment preparation
# -----------------------
# External software: ANTs >= 2.5 and FSL >= 6.0.7.
# Before running: export FSLDIR=/path/to/fsl
# Ensure antsRegistrationSyNQuick.sh, antsApplyTransforms, N4BiasFieldCorrection,
# fslmaths, and bet are on PATH.
#
# Input: quantified CBF NIfTI, subject T1, and MNI template/brain/mask.
# Output: CBF and source-FOV support mask in MNI space plus transforms and QC
# intermediates. The CBF units are not changed by this registration script.

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: register_cbf_to_mni.sh --cbf FILE --t1 FILE --mni-template FILE
       --mni-brain FILE --mni-mask FILE --output-dir DIR
       [--cbf-mask-threshold 0.5] [--threads 4]
USAGE
}

cbf=""; t1=""; mni=""; mni_brain=""; mni_mask=""; output_dir=""
cbf_threshold="0.5"; threads="4"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --cbf) cbf="$2"; shift 2 ;;
    --t1) t1="$2"; shift 2 ;;
    --mni-template) mni="$2"; shift 2 ;;
    --mni-brain) mni_brain="$2"; shift 2 ;;
    --mni-mask) mni_mask="$2"; shift 2 ;;
    --output-dir) output_dir="$2"; shift 2 ;;
    --cbf-mask-threshold) cbf_threshold="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$cbf" && -n "$t1" && -n "$mni" && -n "$mni_brain" && -n "$mni_mask" && -n "$output_dir" ]] || { usage >&2; exit 2; }
for input in "$cbf" "$t1" "$mni" "$mni_brain" "$mni_mask"; do
  [[ -f "$input" ]] || { echo "Missing input: $input" >&2; exit 1; }
done
for command_name in antsRegistrationSyNQuick.sh antsApplyTransforms N4BiasFieldCorrection fslmaths bet; do
  command -v "$command_name" >/dev/null || { echo "Required command not found: $command_name" >&2; exit 1; }
done
if [[ -d "$output_dir" ]] && [[ -n "$(find "$output_dir" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "Refusing to use non-empty output directory: $output_dir" >&2
  exit 1
fi

mkdir -p "$output_dir/anat" "$output_dir/cbf" "$output_dir/registration"
export FSLOUTPUTTYPE=NIFTI_GZ
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS="$threads"
export OMP_NUM_THREADS="$threads"

t1_n4="$output_dir/anat/T1w_desc-N4.nii.gz"
t1_brain_base="$output_dir/anat/T1w_desc-N4_brain"
t1_brain="$t1_brain_base.nii.gz"
t1_brain_mask="$output_dir/anat/T1w_desc-N4_brain_mask.nii.gz"
t1_prefix="$output_dir/registration/T1w_to_MNI_"
cbf_reg="$output_dir/cbf/cbf_registration_input.nii.gz"
cbf_mask="$output_dir/cbf/cbf_registration_mask.nii.gz"
cbf_prefix="$output_dir/registration/CBF_to_T1w_"
support_native="$output_dir/cbf/cbf_source_fov_mask.nii.gz"

N4BiasFieldCorrection -d 3 -i "$t1" -o "$t1_n4" -s 2 -c '[50x50x30x20,1e-6]'
bet "$t1_n4" "$t1_brain_base" -R -f 0.25 -g 0 -m
antsRegistrationSyNQuick.sh -d 3 -f "$mni_brain" -m "$t1_brain" \
  -o "$t1_prefix" -t s -n "$threads" -p f -x "$mni_mask,$t1_brain_mask"

fslmaths "$cbf" -thr 0 -s 1 "$cbf_reg"
fslmaths "$cbf_reg" -thr "$cbf_threshold" -bin -ero "$cbf_mask"
antsRegistrationSyNQuick.sh -d 3 -f "$t1_brain" -m "$cbf_reg" \
  -o "$cbf_prefix" -t a -n "$threads" -p f -x "$t1_brain_mask,$cbf_mask"

antsApplyTransforms -d 3 -i "$cbf" -r "$mni" \
  -o "$output_dir/cbf_space-MNI.nii.gz" -n Linear -f -1 \
  -t "${t1_prefix}1Warp.nii.gz" -t "${t1_prefix}0GenericAffine.mat" \
  -t "${cbf_prefix}0GenericAffine.mat"
fslmaths "$cbf" -mul 0 -add 1 "$support_native"
antsApplyTransforms -d 3 -i "$support_native" -r "$mni" \
  -o "$output_dir/cbf_support_space-MNI.nii.gz" -n NearestNeighbor -f 0 \
  -t "${t1_prefix}1Warp.nii.gz" -t "${t1_prefix}0GenericAffine.mat" \
  -t "${cbf_prefix}0GenericAffine.mat"

echo "Complete: $output_dir/cbf_space-MNI.nii.gz"
