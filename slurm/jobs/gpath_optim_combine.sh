#!/bin/bash
#SBATCH --job-name=gpath_ocmb
#SBATCH --output=slurm/logs/%x_%j.out
#SBATCH --error=slurm/logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --qos=normal

set -euo pipefail

OUTPUT_DIR=${1:?output_dir required}
DELETE_AFTER=${2:-1}

CONTAINER_URI=${CONTAINER_URI:-"docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"}
BINDS=${BINDS:-"-B /home/$USER,/share,/etc/passwd,/etc/group"}

export TZ="US/Eastern"
export APPTAINER_NOHTTPS=1

mkdir -p slurm/logs

echo "PWD: $(pwd)"
echo "Combining optimisation outputs in ${OUTPUT_DIR}"

apptainer exec $BINDS $CONTAINER_URI \
  Rscript pipelines/gpath/optim/combine.R \
  "$OUTPUT_DIR" \
  "$DELETE_AFTER"
