# Multimodal MRI metric extraction code

[![DOI](https://zenodo.org/badge/1363094921.svg)](https://doi.org/10.5281/zenodo.22681564)

This folder contains the public, data-free implementation used to derive the
T1, DTI/CST, resting-state fMRI motor-network, and ASL ROI metrics represented
in Figure A. It does not contain participant images, result tables, private
laterality mappings, software licences, model weights, or access credentials.

## Quick setup

Python 3.10 or newer is recommended:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

The synthetic smoke test was run with Python 3.12, NumPy 2.3.5, pandas 3.0.5,
Nibabel 5.4.2, and Nilearn 0.13.1. The requirements file keeps Nilearn below
0.14 because that release changes confound-standardization behaviour.

External neuroimaging software is installed separately:

- FreeSurfer 7.3 or newer for T1 reconstruction.
- MRtrix3 3.0.4 or newer, FSL 6.0.7 or newer, and ANTs 2.5 or newer for DTI.
- Docker and fMRIPrep 25.2.5 for resting-state fMRI preprocessing.
- FSL 6.0.7 or newer and ANTs 2.5 or newer for ASL-to-MNI registration.

Every executable script repeats its own environment, input, and output
requirements at the top of the file.

## 1. T1: thalamus and precentral-gyrus metrics

Run FreeSurfer, then extract one row of regional metrics:

```bash
bash t1/run_freesurfer_recon.sh \
  --t1 /path/to/sub-XXXX_T1w.nii.gz \
  --subject sub-XXXX \
  --subjects-dir /path/to/freesurfer_subjects

python t1/extract_t1_roi_metrics.py \
  --subject-dir /path/to/freesurfer_subjects/sub-XXXX \
  --subject-id sub-XXXX \
  --lesion-hemisphere L \
  --output-csv /path/to/sub-XXXX_t1_metrics.csv
```

The script reads `aseg.stats` and Desikan-Killiany `lh/rh.aparc.stats`.
Reported measures are bilateral thalamic volume, thalamic volume/eTIV,
precentral mean cortical thickness, thickness SD, surface area, and gray-matter
volume. If lesion hemisphere is supplied, affected/unaffected values, their
ratio, and the asymmetry index
`2 * (affected - unaffected) / (affected + unaffected)` are also reported.

## 2. DTI/CST: FA, MD, AD, and RD

```bash
bash dti/run_dti_cst_pipeline.sh \
  --dwi /path/to/dwi.nii.gz \
  --bval /path/to/dwi.bval \
  --bvec /path/to/dwi.bvec \
  --json /path/to/dwi.json \
  --reverse-pe /path/to/reverse_pe_epi.nii.gz \
  --reverse-pe-json /path/to/reverse_pe_epi.json \
  --output-dir /path/to/dti_output
```

The pipeline applies MP-PCA denoising, Gibbs-ringing removal, eddy correction
(and TOPUP when a reverse-phase image is supplied), ANTs bias correction, and
weighted least-squares tensor fitting. FA is registered to FMRIB58_FA using
affine plus SyN registration. The JHU ICBM tractography maximum-probability
atlas (25% threshold) is transformed to native DWI space with nearest-neighbour
interpolation. Mean FA, MD, AD, and RD are then calculated in left and right
CST masks. The stored JHU label values default to 3 and 4; verify these against
the XML distributed with the installed FSL atlas before analysis.

## 3. Resting-state fMRI: M1/SMA connectivity

Preprocess BIDS data with fMRIPrep:

```bash
bash rsfmri/run_fmriprep_docker.sh \
  --bids-dir /path/to/bids \
  --output-dir /path/to/fmriprep \
  --work-dir /path/to/work \
  --participant-label XXXX \
  --fs-license /path/to/license.txt
```

Create binary AAL1 masks on the atlas grid, transform them to the exact
fMRIPrep BOLD grid using label/nearest-neighbour interpolation, and extract
connectivity:

```bash
python common/make_aal1_masks.py \
  --atlas /path/to/ROI_MNI_V4.nii \
  --output-dir /path/to/aal1_masks

python rsfmri/extract_motor_connectivity.py \
  --bold /path/to/space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz \
  --brain-mask /path/to/space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz \
  --confounds /path/to/desc-confounds_timeseries.tsv \
  --m1-left /path/to/M1_L_in_bold_space.nii.gz \
  --m1-right /path/to/M1_R_in_bold_space.nii.gz \
  --sma-left /path/to/SMA_L_in_bold_space.nii.gz \
  --sma-right /path/to/SMA_R_in_bold_space.nii.gz \
  --lesion-hemisphere L \
  --output-dir /path/to/rsfmri_metrics
```

The primary denoising model uses Friston-24 motion regressors plus the first
five aCompCor components, linear detrending, sample standardization, and a
0.01-0.08 Hz band-pass filter. Frames with FD > 0.5 mm or standardized DVARS
> 1.5 are censored together with one preceding and two following frames; the
first five volumes are also removed. Pearson correlations and Fisher-z values
are written for M1-M1 and ipsilateral/cross-hemisphere M1-SMA pairs. A global
signal regression sensitivity analysis is included.

## 4. ASL: M1 and thalamic CBF

If the CBF map is not already in atlas space, register it first:

```bash
bash asl/register_cbf_to_mni.sh \
  --cbf /path/to/cbf.nii.gz \
  --t1 /path/to/T1w.nii.gz \
  --mni-template /path/to/MNI152_T1_2mm.nii.gz \
  --mni-brain /path/to/MNI152_T1_2mm_brain.nii.gz \
  --mni-mask /path/to/MNI152_T1_2mm_brain_mask.nii.gz \
  --output-dir /path/to/asl_registration

python asl/extract_roi_cbf.py \
  --cbf /path/to/asl_registration/cbf_space-MNI.nii.gz \
  --support-mask /path/to/asl_registration/cbf_support_space-MNI.nii.gz \
  --m1-left /path/to/Precentral_L_mask.nii.gz \
  --m1-right /path/to/Precentral_R_mask.nii.gz \
  --thalamus-left /path/to/Thalamus_L_mask.nii.gz \
  --thalamus-right /path/to/Thalamus_R_mask.nii.gz \
  --units 'mL/100g/min' \
  --lesion-hemisphere L \
  --output-dir /path/to/asl_metrics
```

The extractor reports valid voxel count, coverage, mean, median, sample SD,
5th percentile, and 95th percentile for each ROI. It does not infer CBF units;
pass the units documented by the upstream quantification software. The study
used AAL1 labels 2001/2002 (precentral gyrus), 2401/2402 (SMA), and 7101/7102
(thalamus).

## Atlas and software notes

No atlas or template is redistributed here. Obtain MNI templates and the JHU
atlas from an FSL installation, AAL1 from its official distribution, and
TemplateFlow templates through TemplateFlow/fMRIPrep. Review and cite the
licence and primary publication for every external resource used in a paper.
This repository's MIT licence covers only the code in this folder, not external
software, atlases, templates, or participant data.

## Safety before publishing

`.gitignore` blocks common imaging data, outputs, private configuration, and
FreeSurfer licence files. Still review `git status` before every commit. Example
paths use placeholders, and all study-specific paths and subject lists have
been removed from the executable code.
