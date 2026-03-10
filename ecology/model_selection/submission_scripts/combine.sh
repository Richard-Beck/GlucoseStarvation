#!/bin/bash
#SBATCH --job-name=PF_combine
#SBATCH --output=logs/combine_%j.out
#SBATCH --error=logs/combine_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --qos=normal

MODEL_NAME=${1:-"gpath"}
R_VAL=${2:-2}
P_VAL=${3:-2}
W_VAL=${4:-0}
CONSTRAINT_FLAG=${5:-0}
WASTE_MECH_FLAG=${6:-1}

CONTAINER_URI="docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"
BINDS="-B /home/$USER,/share,/etc/passwd,/etc/group"

export APPTAINER_NOHTTPS=1

apptainer exec $BINDS $CONTAINER_URI \
  Rscript ecology/model_selection/models/gpath/combine.R \
  "$MODEL_NAME" "$R_VAL" "$P_VAL" "$W_VAL" "$CONSTRAINT_FLAG" "$WASTE_MECH_FLAG"