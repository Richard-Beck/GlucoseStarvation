import os
import sys
import argparse
import numpy as np
import cv2
import tifffile
from cellpose.models import CellposeModel

# --- Setup Paths ---
# Adjust this if your imutils path is different
sys.path.insert(0, "/home/4473331/projects/imutils/") 

# --- IMUTILS IMPORTS ---
from imutils.object_classification import ObjectClassifier
from imutils.image_utils import (
    build_raw_group_map,
    get_raw_group_from_key,
    make_cpose_input,
    make_composite,
    robust_normalize
)

# --- Constants ---
LABEL_COLORS = {
    1: (0, 255, 0),    # Alive = Green
    2: (0, 0, 255),    # Dead = Red
    3: (128, 128, 128) # Junk = Gray
}

def load_raw_channels_summed(file_paths):
    """
    Loads raw files and returns the raw aggregates (Sum) without normalization.
    Preserves original bit-depth (e.g. uint16).
    """
    phase, alive, dead = None, [], []
    
    for p in file_paths:
        lower_name = p.lower()
        img = tifffile.imread(p) # Read raw data
        
        if "_phase" in lower_name:
            phase = img
        elif "_alive" in lower_name:
            alive.append(img)
        elif "_dead" in lower_name:
            dead.append(img)
            
    # Summing preserves count information better than Max for raw quantification
    alive_sum = np.sum(alive, axis=0).astype(alive[0].dtype) if alive else None
    dead_sum = np.sum(dead, axis=0).astype(dead[0].dtype) if dead else None
    
    return phase, alive_sum, dead_sum

def create_boundaries_viz(shape_ref, masks, predictions):
    """Generates a black image with just the colored object boundaries"""
    h, w = shape_ref.shape
    vis = np.zeros((h, w, 3), dtype=np.uint8)
    
    for obj_id, pred_data in predictions.items():
        label_id = pred_data['prediction']
        color = LABEL_COLORS.get(label_id, (255, 255, 255))
        
        obj_mask = (masks == obj_id).astype(np.uint8)
        contours, _ = cv2.findContours(obj_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        
        # Draw outlines
        cv2.drawContours(vis, contours, -1, color, 1)
        
    return vis

def main():
    parser = argparse.ArgumentParser(description="Generate intermediate figure panels as TIFFs.")
    parser.add_argument("--raw_data_dir", required=True, help="Path to all_raw")
    parser.add_argument("--image_key", required=True, help="The specific image Key to process")
    parser.add_argument("--classifier_path", required=True, help="Path to object_classifier.pkl")
    parser.add_argument("--output_dir", default="figure_panels", help="Where to save images")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    base_name = os.path.join(args.output_dir, args.image_key)

    # 1. Load File Maps
    print(f"🔎 Building group map from: {args.raw_data_dir}")
    raw_groups = build_raw_group_map(args.raw_data_dir)
    
    print(f"🔎 Searching for files for key: {args.image_key}")
    files = get_raw_group_from_key(args.image_key, raw_groups)
    
    if not files:
        print("❌ No files found!")
        return

    # 2. Save Raw Input Data (Un-normalized)
    # We load these manually here to ensure we get the raw values, not the normalized ones from imutils
    print("💾 Saving Raw Channel Copies...")
    raw_phase, raw_alive, raw_dead = load_raw_channels_summed(files)
    
    if raw_phase is not None: tifffile.imwrite(f"{base_name}_01_raw_phase.tif", raw_phase)
    if raw_alive is not None: tifffile.imwrite(f"{base_name}_02_raw_alive.tif", raw_alive)
    if raw_dead is not None:  tifffile.imwrite(f"{base_name}_03_raw_dead.tif", raw_dead)

    # 3. Save Segmentation Inputs (Normalized)
    print("🧠 Generating & Saving Segmentation Inputs...")
    # This comes from imutils and IS normalized (information loss)
    seg_input = make_cpose_input(args.image_key, raw_groups) # [2, H, W]
    
    # Save the 2-channel stack used by Cellpose
    tifffile.imwrite(f"{base_name}_04_seg_input_normalized.tif", seg_input)

    # 4. Run Cellpose & Save Masks
    print("🧠 Running Cellpose...")
    model = CellposeModel(gpu=True, model_type='cyto')
    masks, _, _ = model.eval(seg_input, diameter=None, channels=[2,1]) 
    
    # Save Masks (uint16)
    tifffile.imwrite(f"{base_name}_05_masks.tif", masks.astype(np.uint16))

    # 5. Save Classifier Inputs (Normalized)
    print("🤖 Generating & Saving Classifier Inputs...")
    # This comes from imutils and IS normalized (information loss)
    clf_input = make_composite(args.image_key, raw_groups) # [H, W, 3] RGB
    
    tifffile.imwrite(f"{base_name}_06_classifier_input_normalized.tif", clf_input)

    # 6. Run Classifier
    print("🤖 Running Classifier...")
    classifier = ObjectClassifier()
    classifier.load_state(args.classifier_path)
    
    # Predict
    # Check if we need predict_with_probabilities or predict_only based on your local version
    # The uploaded file had 'predict_only', but user prompt used 'predict_with_probabilities'. 
    # Sticking to 'predict_only' as it's in the uploaded file, or adapt if your local differs.
    if hasattr(classifier, 'predict_with_probabilities'):
         predictions = classifier.predict_with_probabilities(clf_input, masks)
    else:
         predictions = classifier.predict_only(clf_input, masks)
    
    # 7. Save Boundaries Visualization
    print("🎨 Creating Boundaries Visualization...")
    # Using raw_phase shape for reference
    boundaries_viz = create_boundaries_viz(raw_phase, masks, predictions)
    tifffile.imwrite(f"{base_name}_07_boundaries_viz.tif", boundaries_viz)

    print(f"✅ Done! All Tiffs saved to: {args.output_dir}")

if __name__ == "__main__":
    main()
