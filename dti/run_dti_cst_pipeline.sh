#!/usr/bin/env bash
# Environment preparation
# -----------------------
# External software: MRtrix3 >= 3.0.4, FSL >= 6.0.7, ANTs >= 2.5.
# Python >= 3.10 is used only to read BIDS JSON and write the final CSV.
# Before running: export FSLDIR=/path/to/fsl
# Ensure MRtrix3, "$FSLDIR/bin", and ANTs executables are on PATH.
#
# Input: DWI NIfTI, bval, bvec, BIDS JSON, and optional reverse-PE EPI/JSON.
# Output: preprocessed DWI, brain mask, FA/MD/AD/RD maps, native-space JHU
# CST masks, registration transforms, QC images, and cst_metrics.csv.

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage: run_dti_cst_pipeline.sh --dwi FILE --bval FILE --bvec FILE --json FILE
       --output-dir DIR [--reverse-pe FILE --reverse-pe-json FILE]
       [--standard-fa FILE] [--jhu-atlas FILE]
       [--cst-left-label 3] [--cst-right-label 4]
       [--lesion-hemisphere L|R] [--threads 8]
USAGE
}

dwi=""; bval=""; bvec=""; metadata=""; output_dir=""
reverse_pe=""; reverse_json=""; standard_fa=""; jhu_atlas=""
left_label="3"; right_label="4"; lesion_side=""; threads="8"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dwi) dwi="$2"; shift 2 ;;
    --bval) bval="$2"; shift 2 ;;
    --bvec) bvec="$2"; shift 2 ;;
    --json) metadata="$2"; shift 2 ;;
    --output-dir) output_dir="$2"; shift 2 ;;
    --reverse-pe) reverse_pe="$2"; shift 2 ;;
    --reverse-pe-json) reverse_json="$2"; shift 2 ;;
    --standard-fa) standard_fa="$2"; shift 2 ;;
    --jhu-atlas) jhu_atlas="$2"; shift 2 ;;
    --cst-left-label) left_label="$2"; shift 2 ;;
    --cst-right-label) right_label="$2"; shift 2 ;;
    --lesion-hemisphere) lesion_side="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$dwi" && -n "$bval" && -n "$bvec" && -n "$metadata" && -n "$output_dir" ]] || { usage >&2; exit 2; }
for input in "$dwi" "$bval" "$bvec" "$metadata"; do
  [[ -f "$input" ]] || { echo "Missing input: $input" >&2; exit 1; }
done
if [[ -n "$reverse_pe" || -n "$reverse_json" ]]; then
  [[ -f "$reverse_pe" && -f "$reverse_json" ]] || { echo "Supply both --reverse-pe and --reverse-pe-json" >&2; exit 2; }
fi
[[ -z "$lesion_side" || "$lesion_side" == "L" || "$lesion_side" == "R" ]] || { echo "--lesion-hemisphere must be L or R" >&2; exit 2; }
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "--threads must be a positive integer" >&2; exit 2; }
[[ -n "${FSLDIR:-}" ]] || { echo "FSLDIR is not set" >&2; exit 1; }

standard_fa="${standard_fa:-$FSLDIR/data/standard/FMRIB58_FA_1mm.nii.gz}"
jhu_atlas="${jhu_atlas:-$FSLDIR/data/atlases/JHU/JHU-ICBM-tracts-maxprob-thr25-1mm.nii.gz}"
for input in "$standard_fa" "$jhu_atlas"; do
  [[ -f "$input" ]] || { echo "Missing template/atlas: $input" >&2; exit 1; }
done
for command_name in python3 mrconvert dwidenoise mrdegibbs dwiextract mrmath mrcat \
  dwifslpreproc dwi2mask dwibiascorrect dwi2tensor tensor2metric \
  antsRegistrationSyNQuick.sh antsApplyTransforms fslmaths fslstats; do
  command -v "$command_name" >/dev/null || { echo "Required command not found: $command_name" >&2; exit 1; }
done

if [[ -d "$output_dir" ]] && [[ -n "$(find "$output_dir" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
  echo "Refusing to use non-empty output directory: $output_dir" >&2
  exit 1
fi
mkdir -p "$output_dir/work" "$output_dir/qc" "$output_dir/registration" "$output_dir/scratch"
work="$output_dir/work"
export FSLOUTPUTTYPE=NIFTI_GZ
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS="$threads"
export OMP_NUM_THREADS="$threads"

read -r pe_dir readout_time < <(python3 - "$metadata" <<'PY'
import json
import sys
with open(sys.argv[1], encoding="utf-8") as stream:
    metadata = json.load(stream)
print(metadata["PhaseEncodingDirection"], metadata["TotalReadoutTime"])
PY
)

mrconvert "$dwi" "$work/dwi_raw.mif" -fslgrad "$bvec" "$bval" -json_import "$metadata"
dwidenoise "$work/dwi_raw.mif" "$work/dwi_denoised.mif" \
  -noise "$output_dir/noise_map.nii.gz" -nthreads "$threads"
mrdegibbs "$work/dwi_denoised.mif" "$work/dwi_degibbs.mif" -nthreads "$threads"

if [[ -n "$reverse_pe" ]]; then
  dwiextract "$work/dwi_degibbs.mif" - -bzero | mrmath - mean "$work/b0_same_pe.mif" -axis 3
  mrconvert "$reverse_pe" "$work/reverse_pe.mif" -json_import "$reverse_json"
  mrmath "$work/reverse_pe.mif" mean "$work/b0_reverse_pe.mif" -axis 3
  mrcat "$work/b0_same_pe.mif" "$work/b0_reverse_pe.mif" "$work/b0_pair.mif" -axis 3
  dwifslpreproc "$work/dwi_degibbs.mif" "$work/dwi_preproc.mif" \
    -rpe_pair -se_epi "$work/b0_pair.mif" -pe_dir "$pe_dir" \
    -readout_time "$readout_time" -align_seepi \
    -eddy_options " --repol --cnr_maps --residuals" \
    -eddyqc_text "$output_dir/qc/eddy_text" -scratch "$output_dir/scratch" -nthreads "$threads"
  distortion_method="topup+eddy"
else
  dwifslpreproc "$work/dwi_degibbs.mif" "$work/dwi_preproc.mif" \
    -rpe_none -pe_dir "$pe_dir" -readout_time "$readout_time" \
    -eddy_options " --repol --cnr_maps --residuals" \
    -eddyqc_text "$output_dir/qc/eddy_text" -scratch "$output_dir/scratch" -nthreads "$threads"
  distortion_method="eddy_only"
fi

dwi2mask "$work/dwi_preproc.mif" "$work/mask_preproc.mif" -nthreads "$threads"
dwibiascorrect ants "$work/dwi_preproc.mif" "$work/dwi_biascorr.mif" \
  -mask "$work/mask_preproc.mif" -bias "$output_dir/biasfield.nii.gz" -nthreads "$threads"
dwi2mask "$work/dwi_biascorr.mif" "$work/brain_mask.mif" -nthreads "$threads"
dwi2tensor "$work/dwi_biascorr.mif" "$work/tensor.mif" -mask "$work/brain_mask.mif" -nthreads "$threads"
tensor2metric "$work/tensor.mif" \
  -fa "$output_dir/FA.nii.gz" -adc "$output_dir/MD.nii.gz" \
  -ad "$output_dir/AD.nii.gz" -rd "$output_dir/RD.nii.gz" \
  -vector "$output_dir/V1.nii.gz" -mask "$work/brain_mask.mif" -nthreads "$threads"
mrconvert "$work/dwi_biascorr.mif" "$output_dir/dwi_preprocessed.nii.gz" \
  -export_grad_fsl "$output_dir/dwi_preprocessed.bvec" "$output_dir/dwi_preprocessed.bval"
mrconvert "$work/brain_mask.mif" "$output_dir/brain_mask.nii.gz" -datatype bit
dwiextract "$work/dwi_biascorr.mif" - -bzero | mrmath - mean "$output_dir/mean_b0.nii.gz" -axis 3

reg_prefix="$output_dir/registration/dwi_to_FMRIB58_"
antsRegistrationSyNQuick.sh -d 3 -f "$standard_fa" -m "$output_dir/FA.nii.gz" \
  -o "$reg_prefix" -t s -n "$threads"
antsApplyTransforms -d 3 -i "$jhu_atlas" -r "$output_dir/FA.nii.gz" \
  -o "$output_dir/JHU_tracts_native_dseg.nii.gz" -n NearestNeighbor \
  -t "[${reg_prefix}0GenericAffine.mat,1]" -t "${reg_prefix}1InverseWarp.nii.gz"
fslmaths "$output_dir/JHU_tracts_native_dseg.nii.gz" -thr "$left_label" -uthr "$left_label" -bin "$output_dir/CST_left_mask.nii.gz"
fslmaths "$output_dir/JHU_tracts_native_dseg.nii.gz" -thr "$right_label" -uthr "$right_label" -bin "$output_dir/CST_right_mask.nii.gz"

left_voxels=$(fslstats "$output_dir/CST_left_mask.nii.gz" -V | awk '{print $1}')
right_voxels=$(fslstats "$output_dir/CST_right_mask.nii.gz" -V | awk '{print $1}')
if [[ "$left_voxels" -lt 20 || "$right_voxels" -lt 20 ]]; then
  echo "CST mask too small after registration: left=$left_voxels right=$right_voxels" >&2
  exit 1
fi

values=("$distortion_method" "$left_voxels" "$right_voxels")
for metric in FA MD AD RD; do
  values+=("$(fslstats "$output_dir/$metric.nii.gz" -k "$output_dir/CST_left_mask.nii.gz" -M)")
  values+=("$(fslstats "$output_dir/$metric.nii.gz" -k "$output_dir/CST_right_mask.nii.gz" -M)")
done

python3 - "$output_dir/cst_metrics.csv" "$lesion_side" "${values[@]}" <<'PY'
import csv
import math
import sys

output, lesion_side, distortion, left_voxels, right_voxels, *raw = sys.argv[1:]
metrics = dict(zip(
    ("fa_left", "fa_right", "md_left", "md_right", "ad_left", "ad_right", "rd_left", "rd_right"),
    map(float, raw),
))
row = {
    "distortion_correction": distortion,
    "cst_left_voxels": int(left_voxels),
    "cst_right_voxels": int(right_voxels),
    **metrics,
}
if lesion_side:
    row["lesion_hemisphere"] = lesion_side
    for metric in ("fa", "md", "ad", "rd"):
        affected = metrics[f"{metric}_{'left' if lesion_side == 'L' else 'right'}"]
        unaffected = metrics[f"{metric}_{'right' if lesion_side == 'L' else 'left'}"]
        row[f"{metric}_affected"] = affected
        row[f"{metric}_unaffected"] = unaffected
        row[f"{metric}_affected_to_unaffected_ratio"] = affected / unaffected if unaffected else math.nan
        row[f"{metric}_laterality_index"] = (
            (affected - unaffected) / (affected + unaffected)
            if affected + unaffected else math.nan
        )
with open(output, "w", newline="", encoding="utf-8") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(row))
    writer.writeheader()
    writer.writerow(row)
print(output)
PY

echo "Complete: $output_dir"
