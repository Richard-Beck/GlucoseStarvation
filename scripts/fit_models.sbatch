#!/bin/bash

#SBATCH --job-name=deoptim_fit      # A name for your job
#SBATCH --nodes=1                   # Run all processes on a single node
#SBATCH --ntasks=1                  # Run a single task on a single node
#SBATCH --cpus-per-task=30          # Number of CPU cores per task
#SBATCH --mem=32G                   # Job memory request
#SBATCH --time=24:00:00             # Time limit (DD-HH:MM:SS or HH:MM:SS)
#SBATCH --output=logs/fit_array_%A_%a.out # Standard output and error log. %A=job ID, %a=task ID
#SBATCH --array=1-14                # Create a job array with 14 tasks

echo "======================================================"
echo "Starting SLURM Job Array"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "CPUs per Task: $SLURM_CPUS_PER_TASK"
echo "======================================================"

# --- 1. Define Script and Model ---
# MODIFIED: Set the model you want to run here ("intra_model" or "model_B")
MODEL_NAME="intra_model"

# MODIFIED: Define the new R wrapper script to run
R_SCRIPT="scripts/slurm_fit_wrapper.R"

# This URI was found in your RStudio startup script.
CONTAINER_URI="docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"


# --- 2. Load the Container Runtime Module ---
echo "Loading Apptainer module..."
module purge
module load apptainer

# Create the logs directory if it doesn't exist to prevent errors
mkdir -p logs


# --- 3. Execute the R script INSIDE the container ---
# 'apptainer exec' runs a command within the specified container environment.
# '--bind $(pwd)' makes your current working directory visible inside the container.
echo "Executing R script: ${R_SCRIPT}"
echo "Passing arguments: Task_ID=${SLURM_ARRAY_TASK_ID}, Cores=${SLURM_CPUS_PER_TASK}, Model=${MODEL_NAME}"

apptainer exec \
  --bind $(pwd) \
  "$CONTAINER_URI" \
  Rscript "$R_SCRIPT" "$SLURM_ARRAY_TASK_ID" "$SLURM_CPUS_PER_TASK" "$MODEL_NAME"

EXIT_CODE=$?
echo "======================================================"
echo "Job finished with exit code $EXIT_CODE"
echo "======================================================"

exit $EXIT_CODE

