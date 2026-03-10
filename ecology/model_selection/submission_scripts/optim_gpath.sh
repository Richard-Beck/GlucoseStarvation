#!/bin/bash
#SBATCH --job-name=PF_hier
#SBATCH --output=logs/%x_%A.out
#SBATCH --error=logs/%x_%A.out
#SBATCH --array=1-250
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --qos=normal

MODEL_NAME=${1:-"gpath"}
R_VAL=${2:-2}
P_VAL=${3:-2}
W_VAL=${4:-0}
CONSTRAINT_FLAG=${5:-0}
WASTE_MECH_FLAG=${6:-1}
NUM_THREADS=8

CONTAINER_URI="docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"
BINDS="-B /home/$USER,/share,/etc/passwd,/etc/group"

export OMP_NUM_THREADS=1
export TZ="US/Eastern"

id $USER > /dev/null 2>&1
sleep 2

export APPTAINER_NOHTTPS=1

apptainer exec $BINDS $CONTAINER_URI \
  Rscript ecology/model_selection/models/gpath/optim.R \
  "$MODEL_NAME" "$R_VAL" "$P_VAL" "$W_VAL" "$CONSTRAINT_FLAG" "$WASTE_MECH_FLAG" "$NUM_THREADS" "$SLURM_ARRAY_TASK_ID"