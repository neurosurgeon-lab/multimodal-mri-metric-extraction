#!/usr/bin/env python3
"""Extract QC-aware resting-state M1/SMA ROI-to-ROI connectivity.

Environment preparation
-----------------------
Python >= 3.10
Install: python -m pip install numpy pandas nibabel nilearn
Upstream preprocessing: fMRIPrep 25.2.5, with BOLD, brain mask, and confounds
from the same run and MNI152NLin2009cAsym resolution.

Input
-----
One preprocessed 4D BOLD NIfTI, its fMRIPrep confounds TSV and brain mask, and
four binary AAL1 masks (left/right M1 and SMA) on the exact BOLD grid.

Output
------
``connectivity_metrics.csv``, retained cleaned ROI time series, and a JSON file
recording all denoising/censoring settings. Existing output files are not
overwritten.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
from nilearn.signal import clean


MOTION = ("trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z")


def check_grid(reference: nib.spatialimages.SpatialImage, other: nib.spatialimages.SpatialImage, label: str) -> None:
    if reference.shape[:3] != other.shape[:3]:
        raise ValueError(f"{label} shape mismatch: {other.shape[:3]} != {reference.shape[:3]}")
    if not np.allclose(reference.affine, other.affine, rtol=0, atol=1e-4):
        raise ValueError(f"{label} affine mismatch")


def expanded_censor(fd: np.ndarray, dvars: np.ndarray, initial_drop: int) -> tuple[np.ndarray, np.ndarray]:
    raw = (fd > 0.5) | (dvars > 1.5)
    expanded = raw.copy()
    expanded[:-1] |= raw[1:]
    expanded[1:] |= raw[:-1]
    expanded[2:] |= raw[:-2]
    expanded[:initial_drop] = True
    return raw, expanded


def nuisance_matrix(confounds: pd.DataFrame, include_global_signal: bool) -> np.ndarray:
    columns: list[str] = []
    for base in MOTION:
        columns.extend((base, f"{base}_derivative1", f"{base}_power2", f"{base}_derivative1_power2"))
    acompcor = sorted(column for column in confounds if column.startswith("a_comp_cor_"))[:5]
    if len(acompcor) != 5:
        raise ValueError(f"Expected at least five aCompCor columns, found {len(acompcor)}")
    columns.extend(acompcor)
    if include_global_signal:
        columns.append("global_signal")
    missing = [column for column in columns if column not in confounds]
    if missing:
        raise ValueError(f"Missing nuisance columns: {missing}")
    return confounds[columns].fillna(0).to_numpy(dtype=float)


def correlation(first: np.ndarray, second: np.ndarray) -> tuple[float, float]:
    value = float(np.corrcoef(first, second)[0, 1])
    value = float(np.clip(value, -0.999999, 0.999999))
    return value, float(np.arctanh(value))


def add_pair(row: dict[str, object], label: str, first: np.ndarray, second: np.ndarray, suffix: str = "") -> None:
    value_r, value_z = correlation(first, second)
    row[f"{label}_r{suffix}"] = value_r
    row[f"{label}_fisher_z{suffix}"] = value_z


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bold", type=Path, required=True)
    parser.add_argument("--brain-mask", type=Path, required=True)
    parser.add_argument("--confounds", type=Path, required=True)
    parser.add_argument("--m1-left", type=Path, required=True)
    parser.add_argument("--m1-right", type=Path, required=True)
    parser.add_argument("--sma-left", type=Path, required=True)
    parser.add_argument("--sma-right", type=Path, required=True)
    parser.add_argument("--subject-id", default="anonymous")
    parser.add_argument("--lesion-hemisphere", choices=("L", "R"))
    parser.add_argument("--initial-drop", type=int, default=5)
    parser.add_argument("--min-retained", type=int, default=115)
    parser.add_argument("--min-roi-coverage", type=float, default=0.80)
    parser.add_argument("--tr", type=float, help="Override TR in seconds")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(f"Refusing to use non-empty output directory: {args.output_dir}")
    if args.initial_drop < 0 or args.min_retained < 2:
        raise ValueError("--initial-drop must be >= 0 and --min-retained must be >= 2")

    bold = nib.load(str(args.bold))
    if len(bold.shape) != 4:
        raise ValueError(f"BOLD must be 4D, found {bold.shape}")
    brain_image = nib.load(str(args.brain_mask))
    check_grid(bold, brain_image, "brain mask")
    brain = np.asanyarray(brain_image.dataobj) > 0
    roi_paths = {
        "m1_l": args.m1_left,
        "m1_r": args.m1_right,
        "sma_l": args.sma_left,
        "sma_r": args.sma_right,
    }
    masks: dict[str, np.ndarray] = {}
    coverage: dict[str, float] = {}
    for name, path in roi_paths.items():
        image = nib.load(str(path))
        check_grid(bold, image, name)
        full_mask = np.asanyarray(image.dataobj) > 0
        if not full_mask.any():
            raise ValueError(f"Empty ROI: {path}")
        coverage[name] = float((full_mask & brain).sum() / full_mask.sum())
        masks[name] = full_mask & brain
    minimum_coverage = min(coverage.values())
    if minimum_coverage < args.min_roi_coverage:
        raise ValueError(
            f"Minimum motor ROI coverage {minimum_coverage:.3f} is below {args.min_roi_coverage:.3f}"
        )

    confounds = pd.read_csv(args.confounds, sep="\t")
    if len(confounds) != bold.shape[3]:
        raise ValueError(f"Confound/BOLD length mismatch: {len(confounds)} != {bold.shape[3]}")
    fd = confounds["framewise_displacement"].fillna(0).to_numpy(dtype=float)
    dvars = confounds["std_dvars"].fillna(0).to_numpy(dtype=float)
    raw_censor, censor = expanded_censor(fd, dvars, args.initial_drop)
    retained = np.flatnonzero(~censor)
    if len(retained) < args.min_retained:
        raise ValueError(f"Only {len(retained)} volumes remain; require at least {args.min_retained}")

    bold_data = np.asarray(bold.dataobj, dtype=np.float32)
    signals = np.column_stack([np.nanmean(bold_data[mask], axis=0) for mask in masks.values()])
    tr = args.tr if args.tr is not None else float(bold.header.get_zooms()[3])
    cleaned = clean(
        signals,
        confounds=nuisance_matrix(confounds, include_global_signal=False),
        sample_mask=retained,
        detrend=True,
        standardize="zscore_sample",
        low_pass=0.08,
        high_pass=0.01,
        t_r=tr,
        ensure_finite=True,
    )
    cleaned_gsr = clean(
        signals,
        confounds=nuisance_matrix(confounds, include_global_signal=True),
        sample_mask=retained,
        detrend=True,
        standardize="zscore_sample",
        low_pass=0.08,
        high_pass=0.01,
        t_r=tr,
        ensure_finite=True,
    )
    names = list(masks)
    series = {name: cleaned[:, index] for index, name in enumerate(names)}
    series_gsr = {name: cleaned_gsr[:, index] for index, name in enumerate(names)}
    pairs = {
        "m1_left_right": ("m1_l", "m1_r"),
        "m1_left_sma_left": ("m1_l", "sma_l"),
        "m1_right_sma_right": ("m1_r", "sma_r"),
        "m1_left_sma_right": ("m1_l", "sma_r"),
        "m1_right_sma_left": ("m1_r", "sma_l"),
    }
    if args.lesion_hemisphere:
        affected = args.lesion_hemisphere.lower()
        other = "r" if affected == "l" else "l"
        pairs.update({
            "affected_m1_ipsilesional_sma": (f"m1_{affected}", f"sma_{affected}"),
            "affected_m1_contralesional_sma": (f"m1_{affected}", f"sma_{other}"),
        })

    row: dict[str, object] = {
        "subject_id": args.subject_id,
        "total_volumes": int(bold.shape[3]),
        "tr_seconds": tr,
        "initial_volumes_removed": args.initial_drop,
        "mean_fd_mm": float(np.mean(fd[1:])),
        "median_fd_mm": float(np.median(fd[1:])),
        "maximum_fd_mm": float(np.max(fd[1:])),
        "raw_censor_fraction": float(np.mean(raw_censor)),
        "expanded_censor_fraction": float(np.mean(censor)),
        "retained_volumes": int(len(retained)),
        "retained_duration_minutes": float(len(retained) * tr / 60.0),
        "minimum_motor_roi_coverage_fraction": minimum_coverage,
        **{f"{name}_coverage_fraction": value for name, value in coverage.items()},
    }
    if args.lesion_hemisphere:
        row["lesion_hemisphere"] = args.lesion_hemisphere
    for label, (first, second) in pairs.items():
        add_pair(row, label, series[first], series[second])
        add_pair(row, label, series_gsr[first], series_gsr[second], "_gsr_sensitivity")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "connectivity_metrics.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)
    time_series = pd.DataFrame(cleaned, columns=names)
    for index, name in enumerate(names):
        time_series[f"{name}_gsr_sensitivity"] = cleaned_gsr[:, index]
    time_series.insert(0, "retained_volume_index", retained)
    time_series.to_csv(args.output_dir / "cleaned_roi_timeseries.tsv", sep="\t", index=False)
    settings = {
        "motion_model": "Friston-24",
        "acompcor_components": 5,
        "bandpass_hz": [0.01, 0.08],
        "fd_threshold_mm": 0.5,
        "standardized_dvars_threshold": 1.5,
        "censor_expansion": "one preceding and two following volumes",
        "initial_volumes_removed": args.initial_drop,
        "global_signal_regression_primary": False,
        "global_signal_regression_sensitivity": True,
        "pearson_r_clipped_before_fisher_z": [-0.999999, 0.999999],
    }
    (args.output_dir / "settings.json").write_text(json.dumps(settings, indent=2) + "\n", encoding="utf-8")
    print(args.output_dir / "connectivity_metrics.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
