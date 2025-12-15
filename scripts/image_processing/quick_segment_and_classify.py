import os
import csv
import glob
import argparse
import threading
from concurrent.futures import ThreadPoolExecutor
from tifffile import imread

### Clone this repo and set the path: https://github.com/Richard-Beck/imutils ###
import sys
sys.path.insert(0, "/home/4473331/projects/imutils/")
from imutils.object_classification import ObjectClassifier
from cellpose.models import CellposeModel

# --- Constants ---
MAX_WORKERS = 8 # Number of parallel worker threads
LABEL_MAP = {1: 'alive', 2: 'dead', 3: 'junk'} # Mapping of class IDs to names
csv_writer_lock = threading.Lock() # Lock for thread-safe CSV writing

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
    parser.add_argument("--finetuned_cellpose_model", default=None, help="Path to the fine-tuned Cellpose model file. If not provided, the default model will be used.")
    parser.add_argument("--output_csv_path", required=True, help="Path for the output CSV results file.")
    parser.add_argument("--composite_dir", required=True, help="Directory with 3-channel images for the classifier.")
    parser.add_argument("--cpose_input_dir", required=True, help="Directory with 2-channel images for Cellpose segmentation.")
    return parser.parse_args()

def process_files(cpose_path: str, composite_path: str, cellpose_model: CellposeModel, classifier: ObjectClassifier, csv_writer: csv.DictWriter):
    """
    Defines the processing pipeline for a matched pair of images.
    """
    key = get_base_key(cpose_path)
    try:
        print(f"Processing: {key}")
        
        # 1. Load pre-made images and run segmentation
        cpose_input_img = imread(cpose_path)
        masks, _, _ = cellpose_model.eval(cpose_input_img, diameter=None)

        if masks.max() == 0:
            print(f"  -> No objects found in {key}. Skipping.")
            return 0

        # 2. Load the corresponding composite image for the classifier
        composite_img = imread(composite_path)
        
        # 3. Predict labels for all found objects (masks)
        predictions = classifier.predict_with_probabilities(composite_img, masks)
        
        # 4. Structure results for CSV output
        rows_to_write = []
        for obj_id, result in predictions.items():
            row = {
                'image_key': key,
                'object_id': obj_id,
                'predicted_label_id': result['prediction'],
                'predicted_label_name': LABEL_MAP.get(result['prediction'], 'unknown')
            }
            # Add probability columns
            for i, prob in enumerate(result['probabilities']):
                class_id = i + 1 # Classifier is 0-indexed, labels are 1-indexed
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
    Robust, parallelized batch classification script using pre-made images.
    """
    args = parse_args()
    
    # --- 1. Initialize Models ---
    print("--- Initializing Models ---")
    classifier = ObjectClassifier()
    if not os.path.exists(args.classifier_path):
        print(f"❌ ERROR: Classifier model not found at '{args.classifier_path}'")
        return
    classifier.load_state(args.classifier_path)

    if not classifier.is_trained:
        print("❌ ERROR: The loaded classifier is not trained.")
        return
        
    # Load the Cellpose model
    if args.finetuned_cellpose_model:
        # If a path is provided, check if it exists and load the fine-tuned model
        if not os.path.exists(args.finetuned_cellpose_model):
            print(f"❌ ERROR: Fine-tuned Cellpose model not found at '{args.finetuned_cellpose_model}'")
            return
        print(f"🔬 Loading fine-tuned model from: {args.finetuned_cellpose_model}")
        cellpose_model = CellposeModel(gpu=True, pretrained_model=args.finetuned_cellpose_model)
    else:
        # If no path is provided, load the default model
        print("🔬 No path provided. Loading the default Cellpose model.")
        cellpose_model = CellposeModel(gpu=True)



    # --- 2. Discover and Match Image Files ---
    print("\n--- Scanning and matching image files ---")
    # Create a lookup map of {base_key: full_path} for composite images
    composite_map = {get_base_key(p): p for p in glob.glob(os.path.join(args.composite_dir, "*.tif"))}
    cpose_files = glob.glob(os.path.join(args.cpose_input_dir, "*.tif"))

    # Create a list of tasks (matched pairs of file paths)
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

    # --- 3. Set up CSV Writer and Process in Parallel ---
    print(f"\n--- Starting parallel classification for {len(tasks)} images using {MAX_WORKERS} workers ---")
    
    prob_headers = [f'prob_{name}' for name in sorted(LABEL_MAP.values())]
    headers = ['image_key', 'object_id', 'predicted_label_id', 'predicted_label_name'] + prob_headers
    
    with open(args.output_csv_path, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=headers)
        writer.writeheader()

        total_objects_processed = 0
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            # Submit tasks with unpacked paths
            futures = [executor.submit(process_files, cp, comp_p, cellpose_model, classifier, writer) for cp, comp_p in tasks]
            
            for future in futures:
                total_objects_processed += future.result()

    print("\n--- Batch classification complete ---")
    print(f"✅ Successfully processed {total_objects_processed} objects.")
    print(f"   Results saved to '{args.output_csv_path}'")

if __name__ == "__main__":
    main()
