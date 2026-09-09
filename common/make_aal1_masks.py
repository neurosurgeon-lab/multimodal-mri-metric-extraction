#!/usr/bin/env python3
"""Create binary motor-network and thalamic masks from an AAL1 label image.

Environment preparation
-----------------------
Python >= 3.10
Install: python -m pip install numpy nibabel

Input
-----
An AAL1 ``ROI_MNI_V4`` label NIfTI. The atlas must already be in the spatial
reference required by the downstream analysis.

Output
------
Six uint8 NIfTI masks: bilateral precentral gyrus (M1), supplementary motor
area (SMA), and thalamus. Existing outputs are never overwritten.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np


LABELS = {
    "M1_L": 2001,
    "M1_R": 2002,
    "SMA_L": 2401,
    "SMA_R": 2402,
    "Thalamus_L": 7101,
    "Thalamus_R": 7102,
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, required=True, help="AAL1 label NIfTI")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    atlas = nib.load(str(args.atlas))
    data = np.asanyarray(atlas.dataobj)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, label in LABELS.items():
        output = args.output_dir / f"AAL1_{name}_mask.nii.gz"
        if output.exists():
            raise FileExistsError(f"Refusing to overwrite: {output}")
        mask = (data == label).astype(np.uint8)
        if not mask.any():
            raise ValueError(f"AAL1 label {label} ({name}) is absent from {args.atlas}")
        image = nib.Nifti1Image(mask, atlas.affine, atlas.header)
        image.set_data_dtype(np.uint8)
        nib.save(image, str(output))
        print(f"{name}: label={label}, voxels={int(mask.sum())}, output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
