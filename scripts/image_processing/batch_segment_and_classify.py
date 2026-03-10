import os
import csv
import argparse
import threading
import torch
import numpy as np
import tifffile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from skimage.segmentation import find_boundaries

import sys
sys.path.insert(0, "/home/4473331/projects/imutils/")
from imutils.object_classification import ObjectClassifier
from imutils.image_utils import build_raw_group_map, make_cpose_input, make_composite
from cellpose.models import CellposeModel

# --- Constants ---
MAX_WORKERS = 6 
LABEL_MAP = {1: 'alive', 2: 'dead', 3: 'junk'}
csv_writer_lock = threading.Lock()

def get_imagej_luts():
    """Generates a list of (3, 256) uint8 arrays for a 5-channel ImageJ composite."""
    luts = []
    
    # Ch1: Red (Fluorescence 1 - Swapped from Green)
    lut_red = np.zeros((3, 256), dtype=np.uint8)
    lut_red[0, :] = np.arange(256)
    luts.append(lut_red)
    
    # Ch2: Green (Fluorescence 2 - Swapped from Red)
    lut_green = np.zeros((3, 256), dtype=np.uint8)
    lut_green[1, :] = np.arange(256)
    luts.append(lut_green)
    
    # Ch3: Grayscale (Phase)
    lut_gray = np.zeros((3, 256), dtype=np.uint8)
    lut_gray[0, :] = np.arange(256)
    lut_gray[1, :] = np.arange(256)
    lut_gray[2, :] = np.arange(256)
    luts.append(lut_gray)
    
    # Ch4: Green (Alive Class Segmentations)
    luts.append(lut_green)
    
    # Ch5: Red (Dead Class Segmentations)
    luts.append(lut_red)
    
    return luts

def parse_args() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Batch classify objects directly from raw image data.")
    parser.add_argument("--raw_data_dir", required=True, help="Directory containing the raw image files to process.")
    parser.add_argument("--classifier_path", required=True, help="Path to the trained object_classifier.pkl file.")
    parser.add_argument("--output_csv_path", required=True, help="Path for the output CSV results file.")
    parser.add_argument("--well_ids_path", required=False, help="A file containing well ids, 1 well id on each newline.")
    parser.add_argument("--output_tiff_dir", required=False, help="Directory to save the multi-timepoint TIFFs. Required if well_ids_path is used.")
    return parser.parse_args()

def process_image_key(key: str, raw_groups: dict, cellpose_model: CellposeModel, classifier: ObjectClassifier, csv_writer: csv.DictWriter):
    """
    Standard processing pipeline for a single image key.
    """
    try:
        print(f"Processing: {key}")
        
        cpose_input_img = make_cpose_input(key, raw_groups)
        masks, _, _ = cellpose_model.eval(cpose_input_img, diameter=None)

        if masks.max() == 0:
            print(f"  -> No objects found in {key}. Skipping.")
            return 0

        composite_img_for_classifier = make_composite(key, raw_groups)
        predictions = classifier.predict_with_probabilities(composite_img_for_classifier, masks)
        
        rows_to_write = []
        for obj_id, result in predictions.items():
            row = {
                'image_key': key,
                'object_id': obj_id,
                'predicted_label_id': result['prediction'],
                'predicted_label_name': LABEL_MAP.get(result['prediction'], 'unknown')
            }
            for i, prob in enumerate(result['probabilities']):
                class_id = i + 1
                class_name = LABEL_MAP.get(class_id, f"class_{class_id}")
                row[f'prob_{class_name}'] = f"{prob:.4f}"
            rows_to_write.append(row)

        if rows_to_write:
            with csv_writer_lock:
                csv_writer.writerows(rows_to_write)
        
        return len(rows_to_write)

    except Exception as e:
        print(f"🚨 ERROR processing key {key}: {e}")
        return 0

def process_well_id(well_id: str, keys: list, raw_groups: dict, cellpose_model: CellposeModel, classifier: ObjectClassifier, csv_writer: csv.DictWriter, tiff_output_dir: str):
    print(f"[{well_id}] Starting processing ({len(keys)} timepoints)", flush=True)
    
    all_rows_to_write = []
    tiff_frames = []
    
    try:
        for i, key in enumerate(keys):
            print(f"[{well_id}] Processing frame {i+1}/{len(keys)}: {key}", flush=True)
            
            cpose_input_img = make_cpose_input(key, raw_groups)
            masks, _, _ = cellpose_model.eval(cpose_input_img, diameter=None)
            composite_img_for_classifier = make_composite(key, raw_groups)
            
            # Initialize masks only for Alive (1) and Dead (2)
            mask_alive = np.zeros(composite_img_for_classifier.shape[:2], dtype=np.float32)
            mask_dead = np.zeros(composite_img_for_classifier.shape[:2], dtype=np.float32)

            if masks.max() > 0:
                predictions = classifier.predict_with_probabilities(composite_img_for_classifier, masks)
                boundaries = find_boundaries(masks, mode='inner')
                
                for obj_id, result in predictions.items():
                    pred_label = result['prediction']
                    
                    # Format data for CSV tracking regardless of class
                    row = {
                        'image_key': key,
                        'object_id': obj_id,
                        'predicted_label_id': pred_label,
                        'predicted_label_name': LABEL_MAP.get(pred_label, 'unknown')
                    }
                    for j, prob in enumerate(result['probabilities']):
                        class_id = j + 1
                        class_name = LABEL_MAP.get(class_id, f"class_{class_id}")
                        row[f'prob_{class_name}'] = f"{prob:.4f}"
                    all_rows_to_write.append(row)

                    # Process visual arrays only for classes 1 and 2
                    if pred_label in [1, 2]:
                        obj_mask = (masks == obj_id)
                        obj_boundary = obj_mask & boundaries
                        obj_interior = obj_mask & (~boundaries)
                        
                        # Apply solid boundaries (1.0) and translucent interiors (0.3)
                        if pred_label == 1:
                            mask_alive[obj_boundary] = 1.0
                            mask_alive[obj_interior] = 0.3
                        elif pred_label == 2:
                            mask_dead[obj_boundary] = 1.0
                            mask_dead[obj_interior] = 0.3

            # Expand dimensions to prepare for 5-channel concatenation
            mask_alive_exp = np.expand_dims(mask_alive, axis=-1)
            mask_dead_exp = np.expand_dims(mask_dead, axis=-1)
            
            # Stack composite (3 channels) + Alive (1 channel) + Dead (1 channel) = 5 channels
            frame_stack = np.concatenate([composite_img_for_classifier, mask_alive_exp, mask_dead_exp], axis=-1)
            tiff_frames.append(frame_stack)

        if all_rows_to_write:
            with csv_writer_lock:
                csv_writer.writerows(all_rows_to_write)

        if tiff_frames:
            print(f"[{well_id}] Assembling TIFF from {len(tiff_frames)} frames...", flush=True)
            stack = np.array(tiff_frames, dtype=np.float32) 
            
            # Transpose to (T, C, H, W)
            stack = np.transpose(stack, (0, 3, 1, 2)) 
            
            os.makedirs(tiff_output_dir, exist_ok=True)
            out_tiff = os.path.join(tiff_output_dir, f"{well_id}.tif")
            
            print(f"[{well_id}] Writing 5-channel TIFF to {out_tiff}...", flush=True)
            
            tifffile.imwrite(
                out_tiff, 
                stack, 
                imagej=True, 
                metadata={
                    'axes': 'TCYX',
                    'mode': 'composite',
                    'LUTs': get_imagej_luts()
                }
            )
            
            if os.path.exists(out_tiff):
                size_bytes = os.path.getsize(out_tiff)
                print(f"[{well_id}] ✅ TIFF saved. File size: {size_bytes} bytes", flush=True)
            else:
                print(f"[{well_id}] 🚨 ERROR: File {out_tiff} was not created.", flush=True)

        return len(all_rows_to_write)

    except Exception as e:
        print(f"🚨 ERROR processing well {well_id}: {e}", flush=True)
        import traceback
        traceback.print_exc()
        return 0

def main():
    args = parse_args()
    
    # --- 1. Initialize Models ---
    print("--- Initializing Models ---")
    if not os.path.exists(args.classifier_path):
        print(f"❌ ERROR: Classifier model not found at '{args.classifier_path}'")
        return
        
    classifier = ObjectClassifier()
    classifier.load_state(args.classifier_path)

    if not classifier.is_trained:
        print("❌ ERROR: The loaded classifier is not trained.")
        return
        
    print("⏳ Initializing default Cellpose model...")
    cellpose_model = CellposeModel(gpu=torch.cuda.is_available())

    # --- 2. Discover and Sort Images ---
    print("\n--- Scanning for raw image groups ---")
    raw_groups = build_raw_group_map(args.raw_data_dir)
    image_keys = sorted(raw_groups.keys())
    
    if not image_keys:
        print(f"❌ No raw images found in '{args.raw_data_dir}'.")
        return

    # --- 3. Optional Grouping Logic ---
    target_wells = set()
    run_by_well = False

    if args.well_ids_path:
        if not args.output_tiff_dir:
            print("❌ ERROR: --output_tiff_dir is required when using --well_ids_path.")
            return
            
        run_by_well = True
        with open(args.well_ids_path, 'r') as f:
            target_wells = {line.strip() for line in f if line.strip()}
        
        # Group keys by well_id
        well_to_keys = defaultdict(list)
        for key in image_keys:
            # Split off the timestamp (e.g., "_00d00h00m.tif") to get well_id
            well_id = key.rsplit('_', 1)[0]
            if well_id in target_wells:
                well_to_keys[well_id].append(key)
                
        # Ensure chronological order within each well
        for well in well_to_keys:
            well_to_keys[well].sort()

    # --- 4. Set up CSV Writer and Process in Parallel ---
    print(f"\n--- Starting parallel processing using {MAX_WORKERS} workers ---")
    
    prob_headers = [f'prob_{name}' for name in sorted(LABEL_MAP.values())]
    headers = ['image_key', 'object_id', 'predicted_label_id', 'predicted_label_name'] + prob_headers
    
    total_objects_processed = 0
    
    with open(args.output_csv_path, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=headers)
        writer.writeheader()

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            if run_by_well:
                print(f"Targeting {len(well_to_keys)} specific wells for time-series extraction.")
                futures = [executor.submit(process_well_id, well, keys, raw_groups, cellpose_model, classifier, writer, args.output_tiff_dir) for well, keys in well_to_keys.items()]
            else:
                print(f"Running standard batch classification on {len(image_keys)} total images.")
                futures = [executor.submit(process_image_key, key, raw_groups, cellpose_model, classifier, writer) for key in image_keys]
            
            for future in futures:
                total_objects_processed += future.result()

    print("\n--- Processing complete ---")
    print(f"✅ Successfully processed {total_objects_processed} objects.")
    print(f"   Results saved to '{args.output_csv_path}'")

if __name__ == "__main__":
    main()
