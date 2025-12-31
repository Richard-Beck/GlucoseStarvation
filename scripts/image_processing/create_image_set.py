import os
import argparse
import pandas as pd
import numpy as np
from tifffile import imwrite, imread

# Point this to your library
import sys
sys.path.insert(0, "/home/4473331/projects/imutils/") 
from imutils.image_utils import (
    build_raw_group_map,
    make_composite,      # NOW RETURNS FLOAT32 (Clean RGB)
    make_cpose_input
)

def center_crop(image: np.ndarray, crop_size: int = 200) -> np.ndarray:
    """Crops center region. Handles (H,W,C) and (C,H,W) automatically."""
    if image.ndim == 3:
        # Heuristic: Channels-first if dim0 is smallest (e.g. 3, 2048, 2048)
        if image.shape[0] < image.shape[1] and image.shape[0] < image.shape[2]:
            h, w = image.shape[1], image.shape[2]
            return image[:, (h-crop_size)//2:(h+crop_size)//2, (w-crop_size)//2:(w+crop_size)//2]
        else:
            # Channels-last (H, W, C)
            h, w = image.shape[0], image.shape[1]
            return image[(h-crop_size)//2:(h+crop_size)//2, (w-crop_size)//2:(w+crop_size)//2, :]
    else:
        h, w = image.shape
        return image[(h-crop_size)//2:(h+crop_size)//2, (w-crop_size)//2:(w+crop_size)//2]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw_data_dir', required=True)
    parser.add_argument('--csv_path', required=True)
    parser.add_argument('--composites_dir', required=True)
    parser.add_argument('--cpose_inputs_dir', required=True)
    args = parser.parse_args()

    os.makedirs(args.composites_dir, exist_ok=True)
    os.makedirs(args.cpose_inputs_dir, exist_ok=True)

    df = pd.read_csv(args.csv_path)
    raw_groups = build_raw_group_map(args.raw_data_dir)

    print(f"\n--- Processing {len(df)} images ---")
    
    for index, row in df.iterrows():
        key = row['fname']
        print(f"Processing {index + 1}/{len(df)}: {key}")
        
        try:
            # 1. Generate CLEAN RGB Composite (Float32)
            # This handles the 32-bit vs 8-bit mixing safely
            composite_img = make_composite(key, raw_groups)
            composite_crop = center_crop(composite_img, 200)
            
            # Save as Float32 TIFF (Preserves data integrity for Classifier)
            # Warning: Some standard image viewers might display float32 TIFFs as all black
            # or all white, but ImageJ and Python will read them correctly.
            imwrite(os.path.join(args.composites_dir, f"{key}_composite.tif"), composite_crop)

            # 2. Generate Cellpose Input (Phase + Nuclei Sum)
            cpose_img = make_cpose_input(key, raw_groups)
            cpose_crop = center_crop(cpose_img, 200)
            imwrite(os.path.join(args.cpose_inputs_dir, f"{key}_cpose_input.tif"), cpose_crop)

        except ValueError as e:
            print(f"Skipping {key}: {e}")
        except Exception as e:
            print(f"Error on {key}: {e}")

    print("\n✅ Processing complete!")

if __name__ == "__main__":
    main()
