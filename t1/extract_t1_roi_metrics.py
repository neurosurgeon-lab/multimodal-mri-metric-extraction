#!/usr/bin/env python3
"""Extract thalamic and precentral-gyrus metrics from FreeSurfer outputs.

Environment preparation
-----------------------
Python >= 3.10; only the Python standard library is required.
FreeSurfer must already have completed ``recon-all -all`` for the subject.

Input
-----
The subject's FreeSurfer directory containing ``stats/aseg.stats`` and
``stats/lh.aparc.stats`` / ``stats/rh.aparc.stats``.

Output
------
One CSV row containing bilateral thalamic volume/eTIV and Desikan-Killiany
precentral thickness, surface area, and gray-matter volume. Optional lesion
laterality adds affected/unaffected, ratio, and asymmetry metrics.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


def parse_measures(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    pattern = re.compile(r"^# Measure\s+[^,]+,\s*([^,]+),.*?,\s*([-+0-9.eE]+),")
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = pattern.match(line)
        if match:
            values[match.group(1).strip()] = float(match.group(2))
    return values


def parse_aseg(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) >= 5:
            values[fields[4]] = float(fields[3])
    return values


def parse_aparc_precentral(path: Path) -> dict[str, float]:
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if fields and fields[0] == "precentral" and len(fields) >= 10:
            return {
                "surface_area_mm2": float(fields[2]),
                "gray_matter_volume_mm3": float(fields[3]),
                "mean_thickness_mm": float(fields[4]),
                "thickness_sd_mm": float(fields[5]),
            }
    raise ValueError(f"Precentral region not found in {path}")


def ratio(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator else math.nan


def add_laterality(record: dict[str, object], prefix: str, left: float, right: float, side: str) -> None:
    affected, unaffected = (left, right) if side == "L" else (right, left)
    record[f"{prefix}_affected"] = affected
    record[f"{prefix}_unaffected"] = unaffected
    record[f"{prefix}_affected_to_unaffected_ratio"] = ratio(affected, unaffected)
    denominator = affected + unaffected
    record[f"{prefix}_asymmetry_index"] = (
        2.0 * (affected - unaffected) / denominator if denominator else math.nan
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subject-dir", type=Path, required=True)
    parser.add_argument("--subject-id", help="De-identified output label; defaults to directory name")
    parser.add_argument("--lesion-hemisphere", choices=("L", "R"))
    parser.add_argument("--output-csv", type=Path, required=True)
    args = parser.parse_args()

    if args.output_csv.exists():
        raise FileExistsError(f"Refusing to overwrite: {args.output_csv}")
    aseg_path = args.subject_dir / "stats/aseg.stats"
    lh_path = args.subject_dir / "stats/lh.aparc.stats"
    rh_path = args.subject_dir / "stats/rh.aparc.stats"
    for path in (aseg_path, lh_path, rh_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    measures = parse_measures(aseg_path)
    aseg = parse_aseg(aseg_path)
    lh = parse_aparc_precentral(lh_path)
    rh = parse_aparc_precentral(rh_path)
    etiv = measures["eTIV"]
    thalamus_left = aseg.get("Left-Thalamus", aseg.get("Left-Thalamus-Proper"))
    thalamus_right = aseg.get("Right-Thalamus", aseg.get("Right-Thalamus-Proper"))
    if thalamus_left is None or thalamus_right is None:
        raise ValueError("Bilateral thalamus labels are absent from aseg.stats")

    record: dict[str, object] = {
        "subject_id": args.subject_id or args.subject_dir.name,
        "freesurfer_subject_dir": str(args.subject_dir),
        "estimated_total_intracranial_volume_mm3": etiv,
        "thalamus_left_volume_mm3": thalamus_left,
        "thalamus_right_volume_mm3": thalamus_right,
        "thalamus_left_etiv_ratio": ratio(thalamus_left, etiv),
        "thalamus_right_etiv_ratio": ratio(thalamus_right, etiv),
    }
    for hemisphere, values in (("left", lh), ("right", rh)):
        for metric, value in values.items():
            record[f"m1_precentral_{hemisphere}_{metric}"] = value

    if args.lesion_hemisphere:
        record["lesion_hemisphere"] = args.lesion_hemisphere
        for prefix, left, right in (
            ("thalamus_volume_mm3", thalamus_left, thalamus_right),
            ("thalamus_etiv_ratio", ratio(thalamus_left, etiv), ratio(thalamus_right, etiv)),
            ("m1_precentral_mean_thickness_mm", lh["mean_thickness_mm"], rh["mean_thickness_mm"]),
            ("m1_precentral_surface_area_mm2", lh["surface_area_mm2"], rh["surface_area_mm2"]),
            ("m1_precentral_gray_matter_volume_mm3", lh["gray_matter_volume_mm3"], rh["gray_matter_volume_mm3"]),
        ):
            add_laterality(record, prefix, left, right, args.lesion_hemisphere)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(record))
        writer.writeheader()
        writer.writerow(record)
    print(args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
