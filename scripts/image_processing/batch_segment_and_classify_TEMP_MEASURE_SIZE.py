#!/usr/bin/env python3
import os
import sys
import csv
import argparse
import threading
import random
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import torch  # GPU check

### Clone this repo and set the path: https://github.com/Richard-Beck/imutils ###
sys.path.insert(0, "/home/4473331/projects/imutils/")
from imutils.object_classification import ObjectClassifier
from imutils.image_utils import build_raw_group_map, make_cpose_input, make_composite
from cellpose.models import CellposeModel

# --- Constants ---
MAX_WORKERS = 6
LABEL_MAP = {1: 'alive', 2: 'dead', 3: 'junk'}
SAMPLE_FRACTION = 1  # ~1%

csv_writer_lock = threading.Lock()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Batch classify objects directly from raw image data (1% sample).")
    parser.add_argument("--raw_data_dir", required=True, help="Directory containing the raw image files to process.")
    parser.add_argument("--classifier_path", required=True, help="Path to the trained object_classifier.pkl file.")
    parser.add_argument("--output_csv_path", required=True, help="Path for the output CSV results file.")
    parser.add_argument("--seed", type=int, default=42, help="RNG seed for sampling image subset.")
    parser.add_argument("--max_workers", type=int, default=MAX_WORKERS, help="Thread pool size.")
    parser.add_argument("--sample_fraction", type=float, default=SAMPLE_FRACTION, help="Fraction of images to process.")
    return parser.parse_args()


def process_image_key(
    key: str,
    raw_groups: dict,
    cellpose_model: CellposeModel,
    classifier: ObjectClassifier,
    csv_writer: csv.DictWriter
) -> int:
    """
    Process a single image key:
      1) Build float32 Cellpose input (2, H, W) in [0,1]
      2) Segment with Cellpose
      3) Build float32 composite (H, W, 3) in [0,1] for classifier
      4) Classify each object; write CSV rows (thread-safe)
      5) Include cell_size_px (pixel area) per object
    Returns number of rows written for this image.
    """
    try:
        print(f"Processing: {key}")

        # 1) Cellpose input + segmentation
        cpose_input_img = make_cpose_input(key, raw_groups)  # (2, H, W), float32 [0,1]
        masks, _, _ = cellpose_model.eval(cpose_input_img, diameter=None)

        if masks is None or masks.size == 0 or masks.max() == 0:
            print(f"  -> No objects found in {key}. Skipping.")
            return 0

        # Precompute pixel counts per label id (0 = background)
        counts = np.bincount(masks.ravel())
        if counts.size == 0:
            print(f"  -> Empty mask for {key}. Skipping.")
            return 0

        # 2) Classifier input
        composite_img_for_classifier = make_composite(key, raw_groups)  # (H, W, 3), float32 [0,1]

        # 3) Predict labels + probabilities for all objects
        predictions = classifier.predict_with_probabilities(composite_img_for_classifier, masks)

        # 4) Build rows for CSV
        rows_to_write = []
        for obj_id, result in predictions.items():
            # Safety: ensure obj_id within counts
            size_px = int(counts[obj_id]) if 0 <= obj_id < len(counts) else 0

            row = {
                'image_key': key,
                'object_id': obj_id,
                'predicted_label_id': result['prediction'],
                'predicted_label_name': LABEL_MAP.get(result['prediction'], 'unknown'),
                'cell_size_px': size_px
            }
            # Probability columns: assume index i corresponds to class_id i+1
            for i, prob in enumerate(result['probabilities']):
                class_id = i + 1
                class_name = LABEL_MAP.get(class_id, f"class_{class_id}")
                row[f'prob_{class_name}'] = f"{prob:.4f}"
            rows_to_write.append(row)

        # 5) Write rows (thread-safe)
        if rows_to_write:
            with csv_writer_lock:
                csv_writer.writerows(rows_to_write)

        return len(rows_to_write)

    except Exception as e:
        print(f"🚨 ERROR processing key {key}: {e}")
        return 0


def main():
    args = parse_args()
    random.seed(args.seed)

    # --- Initialize models ---
    print("--- Initializing Models ---")

    if not os.path.exists(args.classifier_path):
        print(f"❌ ERROR: Classifier model not found at '{args.classifier_path}'")
        return

    classifier = ObjectClassifier()
    classifier.load_state(args.classifier_path)
    if not getattr(classifier, "is_trained", False):
        print("❌ ERROR: The loaded classifier is not trained.")
        return

    print(f"⏳ Initializing default Cellpose model (GPU={torch.cuda.is_available()})...")
    cellpose_model = CellposeModel(gpu=torch.cuda.is_available())

    # --- Discover & sample images ---
    print("\n--- Scanning for raw image groups ---")
    raw_groups = build_raw_group_map(args.raw_data_dir)
    all_keys = sorted(raw_groups.keys())

    if not all_keys:
        print(f"❌ No raw images found in '{args.raw_data_dir}'.")
        return

    if not (0 < args.sample_fraction <= 1.0):
        print(f"⚠️ Invalid sample_fraction={args.sample_fraction}; defaulting to 1")
        args.sample_fraction = 1

    sample_n = max(1, int(len(all_keys) * args.sample_fraction))
    image_keys = random.sample(all_keys, sample_n)
    print(f"\n--- Sampling {sample_n} / {len(all_keys)} images (~{args.sample_fraction*100:.2f}%) ---")

    # --- CSV writer ---
    print(f"\n--- Starting parallel classification using {args.max_workers} workers ---")
    prob_headers = [f'prob_{name}' for name in sorted(LABEL_MAP.values())]
    headers = ['image_key', 'object_id', 'predicted_label_id', 'predicted_label_name', 'cell_size_px'] + prob_headers

    os.makedirs(os.path.dirname(os.path.abspath(args.output_csv_path)), exist_ok=True)
    total_objects_processed = 0

    with open(args.output_csv_path, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=headers)
        writer.writeheader()

        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            futures = [
                executor.submit(process_image_key, key, raw_groups, cellpose_model, classifier, writer)
                for key in image_keys
            ]
            for fut in futures:
                total_objects_processed += fut.result()

    print("\n--- Batch classification complete ---")
    print(f"✅ Successfully processed {total_objects_processed} objects.")
    print(f"   Results saved to '{args.output_csv_path}'")


if __name__ == "__main__":
    main()

