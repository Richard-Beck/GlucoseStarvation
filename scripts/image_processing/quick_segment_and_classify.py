import os
import csv
import glob
import argparse
import threading
import torch  # Added for GPU check
import numpy as np
from concurrent.futures import ThreadPoolExecutor
from tifffile import imread, imwrite # Added imwrite to save masks

### Clone this repo and set the path: https://github.com/Richard-Beck/imutils ###
import sys
sys.path.insert(0, "/home/4473331/projects/imutils/")
from imutils.object_classification import ObjectClassifier
from cellpose.models import CellposeModel

# --- Constants ---
MAX_WORKERS = 8  # NOTE: If you encounter CUDA errors, set this to 1
LABEL_MAP = {1: 'alive', 2: 'dead', 3: 'junk'} 
csv_writer_lock = threading.Lock() 

def get_base_key(filepath: str) -> str:
    """Extracts a clean base name from a filepath by stripping common suffixes."""
    stem = os.path.splitext(os.path.basename(filepath))[0]
    suffixes_to_strip = ['_cpose_input', '_composite']
    for suffix in suffixes_to_strip:
        if stem.endswith(suffix):
            return stem[:-len(suffix)]
    return stem

def parse_args() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Batch classify objects in images using a trained classifier.")
    parser.add_argument("--classifier_path", required=True, help="Path to the trained object_classifier.pkl file.")
    # Removed --finetuned_cellpose_model as requested
    parser.add_argument("--output_csv_path", required=True, help="Path for the output CSV results file.")
    parser.add_argument("--output_masks_dir", required=True, help="Directory to save the generated segmentation masks (TIFF).")
    parser.add_argument("--composite_dir", required=True, help="Directory with 3-channel images (Float32) for the classifier.")
    parser.add_argument("--cpose_input_dir", required=True, help="Directory with 2-channel images (Float32) for Cellpose segmentation.")
    return parser.parse_args()

def process_files(cpose_path: str, composite_path: str, output_masks_dir: str, cellpose_model: CellposeModel, classifier: ObjectClassifier, csv_writer: csv.DictWriter):
    """
    Defines the processing pipeline for a matched pair of images.
    """
    key = get_base_key(cpose_path)
    try:
        print(f"Processing: {key}")
        
        # 1. Load Cellpose Input (Float32) and run segmentation
        # Explicit cast to float32 for safety
        cpose_input_img = imread(cpose_path).astype(np.float32)
        
        # Cellpose inference (uses shared model)
        masks, _, _ = cellpose_model.eval(cpose_input_img, diameter=None)

        if masks.max() == 0:
            print(f"  -> No objects found in {key}. Skipping.")
            return 0

        # --- NEW: Save the Segmentation Mask ---
        # We save as uint16 to support >255 objects. 
        # The pixel values in this TIF correspond to the 'object_id' in your CSV.
        mask_filename = f"{key}_mask.tif"
        mask_save_path = os.path.join(output_masks_dir, mask_filename)
        imwrite(mask_save_path, masks.astype(np.uint16))

        # 2. Load the corresponding composite image (Float32)
        # Explicit cast to float32 to match ObjectClassifier requirements
        composite_img = imread(composite_path).astype(np.float32)
        
        # 3. Predict labels for all found objects (masks)
        predictions = classifier.predict_with_probabilities(composite_img, masks)
        
        # 4. Structure results for CSV output
        rows_to_write = []
        for obj_id, result in predictions.items():
            row = {
                'image_key': key,
                'mask_filename': mask_filename, # Added for easier joining later
                'object_id': obj_id,
                'predicted_label_id': result['prediction'],
                'predicted_label_name': LABEL_MAP.get(result['prediction'], 'unknown')
            }
            # Add probability columns
            for i, prob in enumerate(result['probabilities']):
                class_id = i + 1 
                class_name = LABEL_MAP.get(class_id, f"class_{class_id}")
                row[f'prob_{class_name}'] = f"{prob:.4f}"
            rows_to_write.append(row)

        # 5. Write results to CSV in a thread-safe manner
        if rows_to_write:
            with csv_writer_lock:
                csv_writer.writerows(rows_to_write)
        
        return len(rows_to_write)
    except Exception as e:
        print(f"🚨 ERROR processing key {key}: {e}")
        return 0

def main():
    """
    Robust batch classification script using pre-made images.
    """
    args = parse_args()
    
    # --- 1. Setup Directories ---
    os.makedirs(args.output_masks_dir, exist_ok=True)

    # --- 2. Initialize Models ---
    print("--- Initializing Models ---")
    
    # Load Object Classifier
    classifier = ObjectClassifier()
    if not os.path.exists(args.classifier_path):
        print(f"❌ ERROR: Classifier model not found at '{args.classifier_path}'")
        return
    classifier.load_state(args.classifier_path)

    if not classifier.is_trained:
        print("❌ ERROR: The loaded classifier is not trained.")
        return
        
    # Load Default Cellpose Model (User Syntax)
    print("⏳ Initializing default Cellpose model...")
    cellpose_model = CellposeModel(gpu=torch.cuda.is_available())

    # --- 3. Discover and Match Image Files ---
    print("\n--- Scanning and matching image files ---")
    composite_map = {get_base_key(p): p for p in glob.glob(os.path.join(args.composite_dir, "*.tif"))}
    cpose_files = glob.glob(os.path.join(args.cpose_input_dir, "*.tif"))

    tasks = []
    for cpose_path in cpose_files:
        base_key = get_base_key(cpose_path)
        composite_path = composite_map.get(base_key)
        if composite_path:
            tasks.append((cpose_path, composite_path))
        else:
            print(f"  -> Warning: No matching composite image found for '{os.path.basename(cpose_path)}'. Skipping.")

    if not tasks:
        print("❌ No matched image pairs found between the two directories.")
        return

    # --- 4. Set up CSV Writer and Process in Parallel ---
    print(f"\n--- Starting parallel classification for {len(tasks)} images using {MAX_WORKERS} workers ---")
    
    # Sort label names for consistent header order
    prob_headers = [f'prob_{name}' for name in sorted(LABEL_MAP.values())]
    headers = ['image_key', 'mask_filename', 'object_id', 'predicted_label_id', 'predicted_label_name'] + prob_headers
    
    with open(args.output_csv_path, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=headers)
        writer.writeheader()

        total_objects_processed = 0
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            # We now pass args.output_masks_dir to the worker function
            futures = [executor.submit(process_files, cp, comp_p, args.output_masks_dir, cellpose_model, classifier, writer) for cp, comp_p in tasks]
            
            for future in futures:
                total_objects_processed += future.result()

    print("\n--- Batch classification complete ---")
    print(f"✅ Successfully processed {total_objects_processed} objects.")
    print(f"   Results saved to '{args.output_csv_path}'")
    print(f"   Masks saved to '{args.output_masks_dir}'")

if __name__ == "__main__":
    main()
