#!/usr/bin/env python3
"""Extract bilateral M1 and thalamic statistics from a quantified CBF map.

Environment preparation
-----------------------
Python >= 3.10
Install: python -m pip install numpy nibabel
The upstream ASL tool must already have produced a quantified CBF NIfTI. If
needed, use ``register_cbf_to_mni.sh`` before this script.

Input
-----
One 3D CBF map, an optional valid-FOV support mask, and four binary M1/thalamus
ROI masks on exactly the same voxel grid. Supply CBF units explicitly.

Output
------
Long-format bilateral ROI statistics and, when lesion hemisphere is supplied,
affected/unaffected mean CBF, ratio, and asymmetry. Existing outputs are never
overwritten.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import nibabel as nib
import numpy as np


def load_3d(path: Path) -> tuple[nib.spatialimages.SpatialImage, np.ndarray]:
    image = nib.load(str(path))
    if len(image.shape) != 3:
        raise ValueError(f"Expected 3D NIfTI: {path}; found {image.shape}")
    return image, np.asarray(image.dataobj, dtype=np.float32)


def check_grid(reference: nib.spatialimages.SpatialImage, other: nib.spatialimages.SpatialImage, label: str) -> None:
    if reference.shape != other.shape:
        raise ValueError(f"{label} shape mismatch: {other.shape} != {reference.shape}")
    if not np.allclose(reference.affine, other.affine, rtol=0, atol=1e-4):
        raise ValueError(f"{label} affine mismatch")


def roi_statistics(
    cbf: np.ndarray,
    support: np.ndarray,
    mask: np.ndarray,
    region: str,
    hemisphere: str,
    units: str,
) -> dict[str, object]:
    roi = mask > 0.5
    valid = roi & support & np.isfinite(cbf) & (cbf >= 0)
    values = cbf[valid]
    total = int(roi.sum())
    if total == 0:
        raise ValueError(f"Empty ROI: {region}_{hemisphere}")
    if values.size == 0:
        raise ValueError(f"No valid CBF voxels in {region}_{hemisphere}")
    return {
        "region": region,
        "hemisphere": hemisphere,
        "units": units,
        "roi_voxels": total,
        "valid_voxels": int(values.size),
        "coverage_fraction": float(values.size / total),
        "mean_cbf": float(values.mean()),
        "median_cbf": float(np.median(values)),
        "sd_cbf": float(values.std(ddof=1)) if values.size > 1 else 0.0,
        "p05_cbf": float(np.percentile(values, 5)),
        "p95_cbf": float(np.percentile(values, 95)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cbf", type=Path, required=True)
    parser.add_argument("--support-mask", type=Path)
    parser.add_argument("--m1-left", type=Path, required=True)
    parser.add_argument("--m1-right", type=Path, required=True)
    parser.add_argument("--thalamus-left", type=Path, required=True)
    parser.add_argument("--thalamus-right", type=Path, required=True)
    parser.add_argument("--units", required=True, help="For example: mL/100g/min or unconfirmed_vendor_units")
    parser.add_argument("--subject-id", default="anonymous")
    parser.add_argument("--lesion-hemisphere", choices=("L", "R"))
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise FileExistsError(f"Refusing to use non-empty output directory: {args.output_dir}")
    cbf_image, cbf = load_3d(args.cbf)
    if args.support_mask:
        support_image, support_data = load_3d(args.support_mask)
        check_grid(cbf_image, support_image, "support mask")
        support = support_data > 0.5
    else:
        support = np.ones(cbf_image.shape, dtype=bool)

    roi_paths = {
        ("M1", "L"): args.m1_left,
        ("M1", "R"): args.m1_right,
        ("Thalamus", "L"): args.thalamus_left,
        ("Thalamus", "R"): args.thalamus_right,
    }
    rows: list[dict[str, object]] = []
    for (region, hemisphere), path in roi_paths.items():
        image, mask = load_3d(path)
        check_grid(cbf_image, image, f"{region}_{hemisphere}")
        rows.append(roi_statistics(cbf, support, mask, region, hemisphere, args.units))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    long_output = args.output_dir / "roi_cbf_metrics.csv"
    with long_output.open("w", newline="", encoding="utf-8") as stream:
        fieldnames = ["subject_id", *rows[0].keys()]
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({"subject_id": args.subject_id, **row})

    if args.lesion_hemisphere:
        by_key = {(str(row["region"]), str(row["hemisphere"])): row for row in rows}
        other = "R" if args.lesion_hemisphere == "L" else "L"
        laterality_rows = []
        for region in ("M1", "Thalamus"):
            affected = float(by_key[(region, args.lesion_hemisphere)]["mean_cbf"])
            unaffected = float(by_key[(region, other)]["mean_cbf"])
            denominator = affected + unaffected
            laterality_rows.append({
                "subject_id": args.subject_id,
                "region": region,
                "units": args.units,
                "lesion_hemisphere": args.lesion_hemisphere,
                "affected_mean_cbf": affected,
                "unaffected_mean_cbf": unaffected,
                "affected_to_unaffected_ratio": affected / unaffected if unaffected else math.nan,
                "asymmetry_index": 2.0 * (affected - unaffected) / denominator if denominator else math.nan,
            })
        with (args.output_dir / "roi_cbf_laterality.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(laterality_rows[0]))
            writer.writeheader()
            writer.writerows(laterality_rows)

    settings = {
        "valid_voxel_rule": "ROI AND support mask AND finite CBF AND CBF >= 0",
        "units": args.units,
        "lesion_hemisphere": args.lesion_hemisphere,
        "asymmetry_formula": "2 * (affected - unaffected) / (affected + unaffected)",
    }
    (args.output_dir / "settings.json").write_text(json.dumps(settings, indent=2) + "\n", encoding="utf-8")
    print(long_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
