import os
import csv
import argparse
import threading
from concurrent.futures import ThreadPoolExecutor

### Clone this repo and set the path: https://github.com/Richard-Beck/imutils ###
import sys
sys.path.insert(0, "/home/4473331/projects/imutils/")
from imutils.object_classification import ObjectClassifier
from imutils.image_utils import build_raw_group_map, make_cpose_input, make_composite
from cellpose.models import CellposeModel

# --- Constants ---
MAX_WORKERS = 6
LABEL_MAP = {1: 'alive', 2: 'dead', 3: 'junk'}
csv_writer_lock = threading.Lock()

def parse_args() -> argparse.Namespace:
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Batch classify objects from raw image data.")
    parser.add_argument("--raw_data_dir", required=True, help="Directory containing the raw image files to process.")
    parser.add_argument("--classifier_path", required=True, help="Path to the trained object_classifier.pkl file.")
    parser.add_argument("--finetuned_cellpose_model", required=True, help="Path to the fine-tuned Cellpose model file.")
    parser.add_argument("--output_csv_path", required=True, help="Path for the output CSV results file.")
    return parser.parse_args()

def process_image_key(key: str, raw_groups: dict, cellpose_model: CellposeModel, classifier: ObjectClassifier, csv_writer: csv.DictWriter):
    """
    Defines the complete processing pipeline for a single image key,
    generating images from raw data.
    """
    try:
        print(f"Processing: {key}")
        
        # 1. Generate 2-channel image FOR CELLPOSE and get masks
        cpose_input_img = make_cpose_input(key, raw_groups)
        masks, _, _ = cellpose_model.eval(cpose_input_img, diameter=None)

        if masks.max() == 0:
            print(f"  -> No objects found in {key}. Skipping.")
            return 0

        # 2. Generate 3-channel composite image FOR THE CLASSIFIER
        composite_img_for_classifier = make_composite(key, raw_groups)
        
        # 3. Use the classifier to predict labels for all found objects
        predictions = classifier.predict_with_probabilities(composite_img_for_classifier, masks)
        
        # 4. Structure the results for CSV output
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
    Robust, parallelized batch classification script that generates images from raw data.
    """
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
        
    if not os.path.exists(args.finetuned_cellpose_model):
        print(f"❌ ERROR: Fine-tuned Cellpose model not found at '{args.finetuned_cellpose_model}'")
        return
    cellpose_model = CellposeModel(gpu=True, pretrained_model=args.finetuned_cellpose_model)

    # --- 2. Discover and Sort Images ---
    print("\n--- Scanning for raw image groups ---")
    raw_groups = build_raw_group_map(args.raw_data_dir)
    image_keys = sorted(raw_groups.keys())
    
    if not image_keys:
        print(f"❌ No raw images found in '{args.raw_data_dir}'.")
        return

    # --- 3. Set up CSV Writer and Process in Parallel ---
    print(f"\n--- Starting parallel classification for {len(image_keys)} images using {MAX_WORKERS} workers ---")
    
    prob_headers = [f'prob_{name}' for name in sorted(LABEL_MAP.values())]
    headers = ['image_key', 'object_id', 'predicted_label_id', 'predicted_label_name'] + prob_headers
    
    with open(args.output_csv_path, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=headers)
        writer.writeheader()

        total_objects_processed = 0
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = [executor.submit(process_image_key, key, raw_groups, cellpose_model, classifier, writer) for key in image_keys]
            
            for future in futures:
                total_objects_processed += future.result()

    print("\n--- Batch classification complete ---")
    print(f"✅ Successfully processed {total_objects_processed} objects.")
    print(f"   Results saved to '{args.output_csv_path}'")

if __name__ == "__main__":
    main()
