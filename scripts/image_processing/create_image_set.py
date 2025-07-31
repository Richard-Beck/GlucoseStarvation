import os
import argparse
import pandas as pd
import numpy as np
from tifffile import imwrite

###clone this somewhere and set the path to it below: https://github.com/Richard-Beck/imutils ###
import sys
sys.path.insert(0, "/home/4473331/projects/imutils/")
from imutils.image_utils import build_raw_group_map, make_composite, make_cpose_input

def center_crop(image: np.ndarray, crop_size: int = 200) -> np.ndarray:
    """
    Crops a square region from the center of an image.
    Handles (H, W, C), (C, H, W), and (H, W) formats.
    """
    if image.ndim == 3:
        # Determine if format is (C, H, W) or (H, W, C)
        if image.shape[0] < image.shape[1] and image.shape[0] < image.shape[2]:
            h, w = image.shape[1], image.shape[2]
            is_channels_first = True
        else:
            h, w = image.shape[0], image.shape[1]
            is_channels_first = False
    else:  # 2D image
        h, w = image.shape
        is_channels_first = None

    start_y = max(0, (h - crop_size) // 2)
    start_x = max(0, (w - crop_size) // 2)
    end_y = start_y + crop_size
    end_x = start_x + crop_size

    if is_channels_first:
        return image[:, start_y:end_y, start_x:end_x]
    elif is_channels_first is False:
        return image[start_y:end_y, start_x:end_x, :]
    else:  # 2D
        return image[start_y:end_y, start_x:end_x]

def main():
    """
    Loads raw images specified in a CSV, center-crops them, and saves them
    as composite and cellpose_input files.
    """
    # --- Argument Parsing ---
    parser = argparse.ArgumentParser(description="Crop raw images for Cellpose training.")
    parser.add_argument('--raw_data_dir', required=True, help='Path to the directory with raw image files.')
    parser.add_argument('--csv_path', required=True, help='Path to the CSV file with image IDs.')
    parser.add_argument('--composites_dir', required=True, help='Output directory for composite images.')
    parser.add_argument('--cpose_inputs_dir', required=True, help='Output directory for Cellpose input images.')
    args = parser.parse_args()

    # --- Configuration & Setup ---
    CROP_SIZE = 200
    os.makedirs(args.composites_dir, exist_ok=True)
    os.makedirs(args.cpose_inputs_dir, exist_ok=True)

    # --- Initialization ---
    df = pd.read_csv(args.csv_path)
    raw_groups = build_raw_group_map(args.raw_data_dir)

    # --- Image Processing ---
    print(f"\n--- Starting Image Processing for {len(df)} images ---")
    for index, row in df.iterrows():
        key = row['fname']
        print(f"Processing {index + 1}/{len(df)}: {key}")

        # Define output paths
        composite_path = os.path.join(args.composites_dir, f"{key}_composite.tif")
        cpose_input_path = os.path.join(args.cpose_inputs_dir, f"{key}_cpose_input.tif")

        # Load, create, and crop images
        composite_crop = center_crop(make_composite(key, raw_groups), CROP_SIZE)
        cpose_input_crop = center_crop(make_cpose_input(key, raw_groups), CROP_SIZE)

        # Save cropped images to disk
        imwrite(composite_path, composite_crop)
        imwrite(cpose_input_path, cpose_input_crop)

    print("\n✅ Processing complete!")

if __name__ == "__main__":
    main()
