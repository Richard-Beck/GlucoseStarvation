#!/bin/bash
#SBATCH --job-name=PF_nuts
#SBATCH --output=logs/nuts_%A_%a.out
#SBATCH --error=logs/nuts_%A_%a.err

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --qos=normal


# ---- Default: run 4 chains as an array job ----
#SBATCH --array=1-4

# ---- Argument Parsing ----
MODEL_NAME=${1:-"gpath"}
R_VAL=${2:-2}
P_VAL=${3:-2}
W_VAL=${4:-0}
CONSTRAINT_FLAG=${5:-0}      # 1=strict specialization, 0=unconstrained
WASTE_MECH_FLAG=${6:-1}      # 1=additive, 0=multiplicative

# ---- NUTS Config (sensible defaults) ----
CHAIN_ID=${7:-${SLURM_ARRAY_TASK_ID:-1}}   # seed + file id
ITER_WARMUP=${8:-500}
ITER_SAMPLING=${9:-1000}
ADAPT_DELTA=${10:-0.99}
MAX_TREED=${11:-12}

# ---- Threads ----
NUM_THREADS=${12:-16}  # Stan threads per chain

# ---- Container Config ----
CONTAINER_URI="docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"
BINDS="-B /home/$USER,/share,/etc/passwd,/etc/group"

export OMP_NUM_THREADS=1
export TZ="US/Eastern"

echo "Starting NUTS on $(hostname)"
echo "Slurm job: ${SLURM_JOB_ID} | array task: ${SLURM_ARRAY_TASK_ID}"
echo "Model: ${MODEL_NAME} | Dimensions: ${R_VAL}R, ${P_VAL}P, ${W_VAL}W"
echo "Constraint Flag: ${CONSTRAINT_FLAG} | Waste Mech Flag: ${WASTE_MECH_FLAG}"
echo "Chain ID (seed): ${CHAIN_ID}"
echo "NUTS: warmup=${ITER_WARMUP} sampling=${ITER_SAMPLING} adapt_delta=${ADAPT_DELTA} max_treedepth=${MAX_TREED}"
echo "Threads per chain: ${NUM_THREADS}"

# ---- FORCE CACHE REFRESH ----
id $USER > /dev/null 2>&1
sleep 2

# ---- Execution ----
# nuts.R Script Arguments expected:
# 1: MODEL_NAME
# 2: R_val
# 3: P_val
# 4: W_val
# 5: CONSTRAINT_FLAG
# 6: WASTE_MECH_FLAG
# 7: NUM_THREADS
# 8: CHAIN_ID
# 9: ITER_WARMUP
# 10: ITER_SAMPLING
# 11: ADAPT_DELTA
# 12: MAX_TREED
export APPTAINER_NOHTTPS=1
apptainer exec $BINDS $CONTAINER_URI \
  Rscript ecology/model_selection/models/gpath/nuts.R \
  "$MODEL_NAME" \
  "$R_VAL" \
  "$P_VAL" \
  "$W_VAL" \
  "$CONSTRAINT_FLAG" \
  "$WASTE_MECH_FLAG" \
  "$NUM_THREADS" \
  "$CHAIN_ID" \
  "$ITER_WARMUP" \
  "$ITER_SAMPLING" \
  "$ADAPT_DELTA" \
  "$MAX_TREED"