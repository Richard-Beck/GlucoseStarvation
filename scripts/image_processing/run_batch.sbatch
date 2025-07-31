#!/bin/bash

#----------------------------------------------------------------#
# SBATCH - SLURM Resource Directives
#----------------------------------------------------------------#
#SBATCH --job-name=cell_classification   # A descriptive name for your job
#SBATCH --output=logs/classifier-%j.out       # Standard output log file, %j is the job ID
#SBATCH --error=logs/classifier-%j.err        # Standard error log file

#SBATCH --time=2-00:00:00                # Request 24 hours of wall time (1 day, 0 hours)
#SBATCH --qos=large                      # Quality of Service for longer jobs
#SBATCH --nodes=1                        # Request 1 compute node
#SBATCH --ntasks=1                       # Request 1 task (your script is a single process)
#SBATCH --cpus-per-task=6                # Request 6 CPU cores to match the script's workers
#SBATCH --mem=128G                        # Request 32 Gigabytes of memory
#SBATCH --gres=gpu:a30:1                 # Request one A30 GPU

#----------------------------------------------------------------#
# Environment Setup & Job Execution
#----------------------------------------------------------------#

echo "=========================================================="
echo "Job Started at $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURMD_NODENAME"
echo "=========================================================="

# --- IMPORTANT ---
# Activate your dedicated Conda environment.
# This is the only environment setup needed for Python.
# Replace the path below with the correct path to YOUR conda installation.
source /app/eb/software/Anaconda3/2024.02-1/etc/profile.d/conda.sh
conda activate cpose

# --- Run the Python script ---
echo "Python environment activated. Running the script..."
python scripts/image_processing/batch_segment_and_classify.py  --raw_data_dir all_raw --classifier_path data/train/curation_outputs/object_classifier.pkl --finetuned_cellpose_model data/train/models/cpsam-tuned --output_csv_path data/counts/batch_results.csv
echo "=========================================================="
echo "Job Finished at $(date)"
echo "=========================================================="