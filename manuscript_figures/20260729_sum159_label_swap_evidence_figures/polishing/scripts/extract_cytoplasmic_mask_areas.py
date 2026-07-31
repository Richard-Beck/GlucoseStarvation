#!/usr/bin/env python3
"""Extract cell pixel areas from the exact masks used for cytoplasmic red."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("package_dir", type=Path)
    parser.add_argument("cytoplasm_dir", type=Path)
    args = parser.parse_args()

    manifest = pd.read_csv(args.cytoplasm_dir / "sample_manifest.csv")
    records = []
    for row in manifest.itertuples(index=False):
        mask = tifffile.imread(row.cell_mask_path)
        labels, counts = np.unique(mask, return_counts=True)
        keep = labels > 0
        records.extend(
            {
                "image_key": row.image_key,
                "object_id": int(label),
                "cell_area_px": int(count),
            }
            for label, count in zip(labels[keep], counts[keep])
        )

    out = pd.DataFrame.from_records(records)
    out_path = (
        args.package_dir
        / "derived_data"
        / "cytoplasmic_cell_area_from_exact_masks.csv.gz"
    )
    out.to_csv(out_path, index=False, compression="gzip")
    print(
        f"Wrote {len(out):,} image-object areas from "
        f"{manifest['image_key'].nunique()} exact masks."
    )


if __name__ == "__main__":
    main()
