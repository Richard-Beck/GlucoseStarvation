#! /bin/bash
#SBATCH --job-name=PF_plot
#SBATCH --output=logs/plot_%A_%a.out
#SBATCH --error=logs/plot_%A_%a.err

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0:05:00
#SBATCH --qos=small

#SBATCH --array=0-4

# --- Container ---
CONTAINER_URI="docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"
BINDS="-B /home/$USER,/share"

export OMP_NUM_THREADS=1
export TZ="US/Eastern"

# --- Hardcoded inputs ---
MODEL_NAME="model_B1cond"
CELL_LINES=("MCF10A" "MDA-MB-231" "SNU668" "SUM-159-chem" "SUM-159-fuse")

CELL_LINE="${CELL_LINES[$SLURM_ARRAY_TASK_ID]}"

echo "Starting plotting job"
echo "Cell line: ${CELL_LINE}"
echo "Model: ${MODEL_NAME}"
echo "Node: $(hostname)"

apptainer exec $BINDS $CONTAINER_URI \
Rscript MCMC/assess_pathfinder_1cond.R "${CELL_LINE}" "${MODEL_NAME}"
