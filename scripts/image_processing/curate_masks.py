import os
import argparse
from tifffile import imread, imwrite
from cellpose.models import CellposeModel
import torch
from pathlib import Path

### Clone this repo and set the path: https://github.com/Richard-Beck/imutils ###
import sys
sys.path.insert(0, "/home/4473331/projects/imutils/")
from imutils.curation_controller import CurationController

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Generate and curate masks for a dataset using Cellpose and imutils.")
    parser.add_argument("--raw_data_dir", type=str, required=True, 
                        help="Path to the directory with raw images for Cellpose processing.")
    parser.add_argument("--output_dir", type=str, required=True, 
                        help="Path to the base directory for saving masks.")
    parser.add_argument("--display_data_dir", type=str, default=None, 
                        help="(Optional) Path to images for display in the curation tool. "
                             "If not provided, raw_data_dir images are used for display.")
    return parser.parse_args()

def get_base_key(stem: str) -> str:
    """Removes known suffixes from a filename stem to get a base key."""
    suffixes_to_strip = [
        '_cpose_input', 
        '_composite', 
        '_curated_mask', 
        '_cpose_mask'
    ]
    for suffix in suffixes_to_strip:
        if stem.endswith(suffix):
            return stem[:-len(suffix)]
    return stem

def main():
    """
    Main function to run the mask generation and curation pipeline.
    
    This script performs the following steps:
    1. Loads raw images from `raw_data_dir`.
    2. (Optional) Loads a separate set of images from `display_data_dir` for visualization.
    3. Generates initial masks using Cellpose OR loads existing masks.
    4. Launches an interactive curation tool to refine the masks.
    5. Saves the final curated masks to the `output_dir`.
    """
    args = parse_args()

    # --- Setup Paths ---
    raw_data_dir = Path(args.raw_data_dir)
    output_dir = Path(args.output_dir)
    display_data_dir = Path(args.display_data_dir) if args.display_data_dir else raw_data_dir
    
    # Define and create output subdirectories for masks
    cpose_masks_dir = output_dir / "cpose_masks"
    curated_masks_dir = output_dir / "curated_masks"
    for path in [cpose_masks_dir, curated_masks_dir]:
        path.mkdir(parents=True, exist_ok=True)

    # --- Check for Existing Masks ---
    load_existing = False
    if any(cpose_masks_dir.iterdir()) or any(curated_masks_dir.iterdir()):
        while True:
            resp = input("Found existing masks. Load them for curation? (y/n): ").lower()
            if resp in ('y', 'yes'):
                load_existing = True
                break
            elif resp in ('n', 'no'):
                break
            print("Invalid input. Please enter 'y' or 'n'.")
    
    # --- Initialize Model (if needed) ---
    model = None
    if not load_existing:
        print("\n--- Initializing Cellpose model ---")
        model = CellposeModel(gpu=torch.cuda.is_available())

    # --- Phase 1: Batch Processing ---
    print("\n--- Batch processing started ---")
    data_for_curation = []
    
    valid_extensions = {".tif", ".tiff", ".png", ".jpg", ".jpeg"}
    raw_files = sorted([f for f in raw_data_dir.iterdir() if f.suffix.lower() in valid_extensions])
    
    # Create a quick lookup map for display images
    display_files_map = {get_base_key(f.stem): f for f in display_data_dir.iterdir() if f.suffix.lower() in valid_extensions}
    for idx, raw_path in enumerate(raw_files, start=1):
        key = get_base_key(raw_path.stem)
        print(f"Processing {idx}/{len(raw_files)}: {raw_path.name}")

        # Define paths for mask files
        cpose_path = cpose_masks_dir / f"{key}.tif"
        curated_path = curated_masks_dir / f"{key}.tif"

        # Determine the image to display in the curation tool
        display_path = display_files_map.get(key)
        if not display_path:
            print(f"  → Warning: No corresponding display image found for '{key}'. Skipping this file.")
            continue
        display_image = imread(display_path)
        
        initial_mask = None
        mask_source = "newly generated"

        if load_existing:
            if curated_path.exists():
                initial_mask = imread(curated_path)
                mask_source = "existing curated mask"
            elif cpose_path.exists():
                initial_mask = imread(cpose_path)
                mask_source = "existing Cellpose mask"

        # If we aren't loading existing masks or none were found, run Cellpose on the raw image
        if initial_mask is None:
            if model is None: # Lazy initialization
                print("\n--- Initializing Cellpose model for remaining images ---")
                model = CellposeModel(gpu=torch.cuda.is_available())
            
            raw_image = imread(raw_path)
            print("  → Running Cellpose to generate initial mask...")
            masks, _, _ = model.eval(raw_image, diameter=None)
            imwrite(cpose_path, masks)
            print(f"  → Saved initial Cellpose mask to '{cpose_path.name}'")
            initial_mask = masks
        else:
             print(f"  → Loaded {mask_source}.")

        # Queue the display image and its mask for curation
        data_for_curation.append({
            "image": display_image,
            "masks": initial_mask,
            "title": key,
            "path": curated_path
        })

    # --- Phase 2: Interactive Curation ---
    if not data_for_curation:
        print("\nNo images to process or curate. All done! ✅")
        return

    print(f"\n--- Launching curation for {len(data_for_curation)} images ---")
    controller = CurationController(
        images=[d["image"] for d in data_for_curation],
        initial_masks=[d["masks"] for d in data_for_curation],
        titles=[d["title"] for d in data_for_curation],
        mask_save_paths=[d["path"] for d in data_for_curation]
    )
    controller.start()
    print("\n🎉 Curation session finished!")

if __name__ == "__main__":
    main()
