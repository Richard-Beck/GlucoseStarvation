# main_script.py  <-- replace your current main script with this (only changed file)

import os
import argparse
import pandas as pd
import numpy as np
from tifffile import imwrite, imread

### clone this somewhere and set the path to it below: https://github.com/Richard-Beck/imutils ###
import sys
sys.path.insert(0, "/home/4473331/projects/imutils/")
from imutils.image_utils import (
    build_raw_group_map,
    make_composite,
    make_cpose_input,
    robust_normalize,     # used to normalize each fluorescence channel
)

def center_crop(image: np.ndarray, crop_size: int = 200) -> np.ndarray:
    """
    Crops a square region from the center of an image.
    Handles (H, W, C), (C, H, W), and (H, W) formats.
    """
    if image.ndim == 3:
        # Determine if format is (C, H, W) or (H, W, C)
        # crude but effective test: if first axis smaller than the others -> channels first
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

# image_utils.py

# ... (keep all existing functions like build_raw_group_map, make_composite, etc.) ...

def make_multichannel(root, raw_groups):
    """
    Creates a 5-channel image from raw TIFFs, preserving raw pixel values.
    
    The function ensures a consistent channel order and bit depth (uint16).
    If a specific channel is missing for a group, it's represented as a plane of zeros.
    
    Channel Order:
    1. phase
    2. alive
    3. aliveC1
    4. aliveC2
    5. dead
    """
    # Define the canonical order and names for the 5 output channels
    CHANNEL_ORDER = ['phase', 'alive', 'aliveC1', 'aliveC2', 'dead']
    
    # Use a dictionary to group images by their identified channel type
    raw_imgs_by_channel = {name: [] for name in CHANNEL_ORDER}
    
    # Sort file paths to find the correct channel type for each image
    for path in raw_groups.get(root, []):
        filename = os.path.basename(path)
        channel_type = filename.split('_')[-4].lower()
        
        # Handle the general 'alive' case vs specific 'aliveC1'/'aliveC2'
        if channel_type in raw_imgs_by_channel:
             raw_imgs_by_channel[channel_type].append(imread(path))
        elif 'alive' in channel_type and 'alive' in raw_imgs_by_channel:
             raw_imgs_by_channel['alive'].append(imread(path))
             
    # --- Create the output image stack ---
    
    # 1. Determine the image dimensions (H, W) from the first available image
    h, w = -1, -1
    for channel_list in raw_imgs_by_channel.values():
        if channel_list:
            h, w = channel_list[0].shape
            break
    if h == -1: # If no images were found at all for this key
        raise ValueError(f"No image files found for key: {root}")

    # 2. Create an empty 5-channel canvas. Shape is (C, H, W).
    # We use uint16 as a standard, robust bit depth for microscopy.
    multichannel_stack = np.zeros((len(CHANNEL_ORDER), h, w), dtype=np.uint16)
    
    # 3. Fill the canvas with data from the channels we found
    for i, channel_name in enumerate(CHANNEL_ORDER):
        images = raw_imgs_by_channel[channel_name]
        if images:
            # Project the maximum intensity if there's a z-stack
            max_projected_img = np.maximum.reduce(images)
            # Place it in the correct channel, ensuring consistent bit depth
            multichannel_stack[i, :, :] = max_projected_img.astype(np.uint16)
            
    return multichannel_stack

def main():
    """
    Loads raw images specified in a CSV, center-crops them, and saves them
    as composite, cellpose_input, and optional multi-channel TIFF files.
    """
    # --- Argument Parsing ---
    parser = argparse.ArgumentParser(description="Crop raw images for Cellpose training.")
    parser.add_argument('--raw_data_dir', required=True, help='Path to the directory with raw image files.')
    parser.add_argument('--csv_path', required=True, help='Path to the CSV file with image IDs.')
    parser.add_argument('--composites_dir', required=True, help='Output directory for composite images.')
    parser.add_argument('--cpose_inputs_dir', required=True, help='Output directory for Cellpose input images.')
    parser.add_argument('--multichannel_dir', required=False, default=None,
                        help='(Optional) Output directory for per-group multi-channel TIFFs (fluorescent channels).')
    args = parser.parse_args()

    # --- Configuration & Setup ---
    CROP_SIZE = 200
    os.makedirs(args.composites_dir, exist_ok=True)
    os.makedirs(args.cpose_inputs_dir, exist_ok=True)
    if args.multichannel_dir:
        os.makedirs(args.multichannel_dir, exist_ok=True)

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
        multichannel_path = None
        if args.multichannel_dir:
            multichannel_path = os.path.join(args.multichannel_dir, f"{key}_multichannel.tif")

        # Load, create, and crop images
        composite_crop = center_crop(make_composite(key, raw_groups), CROP_SIZE)
        cpose_input_crop = center_crop(make_cpose_input(key, raw_groups), CROP_SIZE)

        # Save cropped composite and cellpose inputs
        imwrite(composite_path, composite_crop)
        imwrite(cpose_input_path, cpose_input_crop)

        # Optional: create and save multichannel TIFF (each fluorescence as separate channel)
        if multichannel_path:
            try:
                multichannel = make_multichannel(key, raw_groups)
                multichannel_crop = center_crop(multichannel, CROP_SIZE)
                imwrite(multichannel_path, multichannel_crop)
            except ValueError as e:
                # no fluorescence channels found for this key; report and continue
                print(f"Skipping multichannel for {key}: {e}")

    print("\n✅ Processing complete!")

if __name__ == "__main__":
    main()

