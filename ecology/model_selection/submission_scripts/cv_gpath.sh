#!/bin/bash
#SBATCH --job-name=cv_nuts
#SBATCH --output=logs/nuts_cv_%A_%a.out
#SBATCH --error=logs/nuts_cv_%A_%a.err

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --qos=normal
#SBATCH --array=1-20
# ---- CV grid: folds x chains ----
# Default: 5 folds (cell lines) x 4 chains = 20 tasks
N_FOLDS=${13:-5}
N_CHAINS=${14:-4}


# ---- Argument Parsing (same as before) ----
MODEL_NAME=${1:-"gpath"}
R_VAL=${2:-1}
P_VAL=${3:-1}
W_VAL=${4:-0}
CONSTRAINT_FLAG=${5:-0}      # 1=strict specialization, 0=unconstrained
WASTE_MECH_FLAG=${6:-1}      # 1=additive, 0=multiplicative

# ---- NUTS Config (defaults) ----
ITER_WARMUP=${7:-500}
ITER_SAMPLING=${8:-1000}
ADAPT_DELTA=${9:-0.99}
MAX_TREED=${10:-12}

# ---- Threads ----
NUM_THREADS=${11:-16}

# ---- Script path ----
CV_SCRIPT=${12:-"ecology/model_selection/models/gpath/nuts_cv.R"}

# ---- Derive FOLD_ID and CHAIN_ID from SLURM_ARRAY_TASK_ID ----
# Task IDs 1..(N_FOLDS*N_CHAINS)
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}

MAX_TASKS=$((N_FOLDS * N_CHAINS))
if [ "$TASK_ID" -lt 1 ] || [ "$TASK_ID" -gt "$MAX_TASKS" ]; then
  echo "ERROR: TASK_ID=$TASK_ID out of range 1..$MAX_TASKS (N_FOLDS=$N_FOLDS, N_CHAINS=$N_CHAINS)"
  exit 1
fi

# 1-based indexing:
# fold = 1 + floor((task-1)/N_CHAINS)
# chain = 1 + ((task-1) % N_CHAINS)
FOLD_ID=$(( 1 + ( (TASK_ID - 1) / N_CHAINS ) ))
CHAIN_ID=$(( 1 + ( (TASK_ID - 1) % N_CHAINS ) ))

# ---- Container Config ----
CONTAINER_URI="docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"
BINDS="-B /home/$USER,/share,/etc/passwd,/etc/group"

export OMP_NUM_THREADS=1
export TZ="US/Eastern"

echo "Starting CV NUTS on $(hostname)"
echo "Slurm job: ${SLURM_JOB_ID} | array task: ${TASK_ID} (fold=${FOLD_ID}, chain=${CHAIN_ID})"
echo "Model: ${MODEL_NAME} | Dimensions: ${R_VAL}R, ${P_VAL}P, ${W_VAL}W"
echo "Constraint Flag: ${CONSTRAINT_FLAG} | Waste Mech Flag: ${WASTE_MECH_FLAG}"
echo "NUTS: warmup=${ITER_WARMUP} sampling=${ITER_SAMPLING} adapt_delta=${ADAPT_DELTA} max_treedepth=${MAX_TREED}"
echo "Threads per chain: ${NUM_THREADS}"
echo "CV script: ${CV_SCRIPT}"

# ---- FORCE CACHE REFRESH ----
id $USER > /dev/null 2>&1
sleep 2

# ---- Execution ----
# run_nuts_cv_single_chain.R expected args:
# 1: MODEL_NAME
# 2: R_val
# 3: P_val
# 4: W_val
# 5: CONSTRAINT_FLAG
# 6: WASTE_MECH_FLAG
# 7: NUM_THREADS
# 8: FOLD_ID
# 9: CHAIN_ID
# 10: ITER_WARMUP
# 11: ITER_SAMPLING
# 12: ADAPT_DELTA
# 13: MAX_TREED
export APPTAINER_NOHTTPS=1
apptainer exec $BINDS $CONTAINER_URI \
  Rscript "$CV_SCRIPT" \
  "$MODEL_NAME" \
  "$R_VAL" \
  "$P_VAL" \
  "$W_VAL" \
  "$CONSTRAINT_FLAG" \
  "$WASTE_MECH_FLAG" \
  "$NUM_THREADS" \
  "$FOLD_ID" \
  "$CHAIN_ID" \
  "$ITER_WARMUP" \
  "$ITER_SAMPLING" \
  "$ADAPT_DELTA" \
  "$MAX_TREED"