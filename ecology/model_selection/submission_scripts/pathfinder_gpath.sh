#!/bin/bash
#SBATCH --job-name=PF_pathfinder
#SBATCH --output=logs/pathfinder_%j.out
#SBATCH --error=logs/pathfinder_%j.err

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --qos=normal

# ---- Argument Parsing ----
MODEL_NAME=${1:-"gpath"}
R_VAL=${2:-2}
P_VAL=${3:-2}
W_VAL=${4:-0}
CONSTRAINT_FLAG=${5:-0}
WASTE_MECH_FLAG=${6:-1}

# ---- Pathfinder Config ----
SEED=${7:-1}
PF_DRAWS=${8:-1000}
PF_SINGLE_DRAWS=${9:-1000}
PF_MAX_LBFGS=${10:-1000}
PF_INIT_ALPHA=${11:-0.01}
PF_NUM_PATHS=${12:-4}

# ---- Threads ----
NUM_THREADS=${13:-16}

# ---- Container Config ----
CONTAINER_URI="docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"
BINDS="-B /home/$USER,/share,/etc/passwd,/etc/group"

export OMP_NUM_THREADS=1
export TZ="US/Eastern"
export APPTAINER_NOHTTPS=1

echo "Starting Pathfinder on $(hostname)"
echo "Slurm job: ${SLURM_JOB_ID}"
echo "Model: ${MODEL_NAME} | Dimensions: ${R_VAL}R, ${P_VAL}P, ${W_VAL}W"
echo "Constraint Flag: ${CONSTRAINT_FLAG} | Waste Mech Flag: ${WASTE_MECH_FLAG}"
echo "Seed: ${SEED}"
echo "Pathfinder: draws=${PF_DRAWS} single_path_draws=${PF_SINGLE_DRAWS} max_lbfgs_iters=${PF_MAX_LBFGS} init_alpha=${PF_INIT_ALPHA} num_paths=${PF_NUM_PATHS}"
echo "Threads: ${NUM_THREADS}"

# ---- FORCE CACHE REFRESH ----
id $USER > /dev/null 2>&1
sleep 2

# ---- Execution ----
# pathfinder.R arguments:
# 1: MODEL_NAME
# 2: R_val
# 3: P_val
# 4: W_val
# 5: CONSTRAINT_FLAG
# 6: WASTE_MECH_FLAG
# 7: NUM_THREADS
# 8: SEED
# 9: PF_DRAWS
# 10: PF_SINGLE_DRAWS
# 11: PF_MAX_LBFGS
# 12: PF_INIT_ALPHA
# 13: PF_NUM_PATHS
apptainer exec $BINDS $CONTAINER_URI \
  Rscript ecology/model_selection/models/gpath/pathfinder.R \
  "$MODEL_NAME" \
  "$R_VAL" \
  "$P_VAL" \
  "$W_VAL" \
  "$CONSTRAINT_FLAG" \
  "$WASTE_MECH_FLAG" \
  "$NUM_THREADS" \
  "$SEED" \
  "$PF_DRAWS" \
  "$PF_SINGLE_DRAWS" \
  "$PF_MAX_LBFGS" \
  "$PF_INIT_ALPHA" \
  "$PF_NUM_PATHS"