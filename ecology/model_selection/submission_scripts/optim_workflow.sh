#!/bin/bash
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
MODEL_NAME=${1:-"gpath"}
R_VAL=${2:-2}
P_VAL=${3:-2}
W_VAL=${4:-0}
CONSTRAINT_FLAG=${5:-0}
WASTE_MECH_FLAG=${6:-1}

# Submit array and capture ID
OPTIM_JOB_ID=$(sbatch --parsable ecology/model_selection/submission_scripts/optim_gpath.sh "$MODEL_NAME" "$R_VAL" "$P_VAL" "$W_VAL" "$CONSTRAINT_FLAG" "$WASTE_MECH_FLAG")
echo "Optimization array submitted: $OPTIM_JOB_ID"

# Submit consolidation with afterok dependency
COMBINE_JOB_ID=$(sbatch --parsable --dependency=afterany:$OPTIM_JOB_ID ecology/model_selection/submission_scripts/combine.sh "$MODEL_NAME" "$R_VAL" "$P_VAL" "$W_VAL" "$CONSTRAINT_FLAG" "$WASTE_MECH_FLAG")
echo "Consolidation job submitted: $COMBINE_JOB_ID"