import os
import glob
import json
import argparse
import pickle
import numpy as np
from tifffile import imread

### Clone this repo and set the path: https://github.com/Richard-Beck/imutils ###
import sys
sys.path.insert(0, "/home/4473331/projects/imutils/")
from imutils.curation_controller import CurationController
from imutils.object_classification import ObjectClassifier

def parse_args() -> argparse.Namespace:
    """Parses command-line arguments for the labeling script."""
    parser = argparse.ArgumentParser(description="Interactive object labeling and classification.")
    parser.add_argument("--images_dir", required=True, help="Directory containing the images for display.")
    parser.add_argument("--masks_dir", required=True, help="Directory containing the mask files to be labeled.")
    parser.add_argument("--output_dir", required=True, help="Directory to save session labels and the classifier model.")
    return parser.parse_args()

def get_base_key(stem: str) -> str:
    """Removes known suffixes from a filename stem to get a base key."""
    # List of suffixes to check for and remove from the end of the filename
    suffixes_to_strip = [
        '_composite', 
        '_curated_mask',
        '_cpose_input_mask',
        '_cpose_mask'
    ]
    for suffix in suffixes_to_strip:
        if stem.endswith(suffix):
            # Return the stem with the suffix removed
            return stem[:-len(suffix)]
    # If no suffixes match, return the original stem
    return stem

def main_labelling():
    """
    Main script to launch the interactive object labeling session.
    """
    args = parse_args()

    # --- Setup Paths from Arguments ---
    os.makedirs(args.output_dir, exist_ok=True)
    SESSION_LABELS_PATH = os.path.join(args.output_dir, "session_labels.json")
    CLASSIFIER_PATH = os.path.join(args.output_dir, "object_classifier.pkl")

    # --- Load Data ---
    print(f"🔎 Finding images in: {args.images_dir}")
    print(f"🔎 Finding masks in:  {args.masks_dir}")

    # Create a lookup map for images using their base key
    image_map = {get_base_key(os.path.splitext(os.path.basename(p))[0]): p for p in glob.glob(os.path.join(args.images_dir, "*.tif"))}
    mask_paths = sorted(glob.glob(os.path.join(args.masks_dir, "*.tif")))

    if not mask_paths:
        print(f"❌ No masks found in '{args.masks_dir}'.")
        return

    # Load previous session labels if they exist
    loaded_labels = None
    if os.path.exists(SESSION_LABELS_PATH):
        print(f"💾 Found existing session labels at {SESSION_LABELS_PATH}. Loading...")
        with open(SESSION_LABELS_PATH, 'r') as f:
            loaded_labels = [{int(k): v for k, v in d.items()} for d in json.load(f)]

    # --- Match Images to Masks and Prepare Data for Controller ---
    images_to_load, masks_to_load, labels_to_load, titles_to_load = [], [], [], []
    for i, mask_path in enumerate(mask_paths):
        mask_stem = os.path.splitext(os.path.basename(mask_path))[0]
        base_key = get_base_key(mask_stem)
        
        # Find the corresponding image using the base key
        image_path = image_map.get(base_key)
        if not image_path:
            print(f"  → Warning: No matching image found for mask '{os.path.basename(mask_path)}'. Skipping.")
            continue

        images_to_load.append(imread(image_path))
        masks_to_load.append(imread(mask_path))
        titles_to_load.append(base_key)

        # Load existing labels or create new empty ones
        if loaded_labels and i < len(loaded_labels):
            labels_to_load.append(loaded_labels[i])
        else:
            mask = masks_to_load[-1]
            labels_to_load.append({int(obj_id): 0 for obj_id in np.unique(mask) if obj_id != 0})

    if not images_to_load:
        print("❌ No valid image/mask pairs were found. Exiting.")
        return

    # --- Initialize Classifier and Load State ---
    classifier = ObjectClassifier()
    if os.path.exists(CLASSIFIER_PATH):
        try:
            classifier.load_state(CLASSIFIER_PATH)
        except (pickle.UnpicklingError, KeyError, EOFError) as e:
            print(f"⚠️ Could not load classifier state due to error: {e}. Starting with a fresh classifier.")

    # --- Launch Controller ---
    print(f"\n🚀 Launching interactive session for {len(images_to_load)} images...")
    controller = CurationController(
        images=images_to_load,
        initial_masks=masks_to_load,
        titles=titles_to_load,
        initial_labels=labels_to_load,
        object_classifier=classifier,
        session_path=SESSION_LABELS_PATH,
        classifier_path=CLASSIFIER_PATH
    )
    controller.start()

    print("\n🎉 Labelling session finished!")

if __name__ == "__main__":
    main_labelling()
