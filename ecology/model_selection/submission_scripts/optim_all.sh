#!/bin/bash

# Target run IDs
RUN_IDS=(
  "1R_1P_0W_C0_M1"
  "1R_1P_1W_C0_M0"
  "1R_1P_1W_C0_M1"
  "2R_1P_0W_C0_M1"
  "2R_1P_1W_C0_M0"
  "2R_1P_1W_C0_M1"
  "2R_2P_0W_C0_M1"
  "2R_2P_0W_C1_M0"
  "2R_2P_1W_C0_M0"
  "2R_2P_1W_C0_M1"
  "2R_2P_1W_C1_M0"
  "2R_2P_1W_C1_M1"
)

MODEL_NAME="gpath"
SCRIPT_PATH="ecology/model_selection/submission_scripts/optim_workflow.sh"

echo "Submitting ${#RUN_IDS[@]} jobs to SLURM..."

for run_id in "${RUN_IDS[@]}"; do
  # Split the string by underscore
  IFS='_' read -r r_part p_part w_part c_part m_part <<< "$run_id"
  
  # Strip the alphabetic characters to isolate the integers
  R_VAL=${r_part//R/}
  P_VAL=${p_part//P/}
  W_VAL=${w_part//W/}
  C_VAL=${c_part//C/}
  M_VAL=${m_part//M/}
  
  echo "Submitting: $run_id"
  
  # Submit job with extracted arguments
  sbatch "$SCRIPT_PATH" "$MODEL_NAME" "$R_VAL" "$P_VAL" "$W_VAL" "$C_VAL" "$M_VAL"
  
  # Brief pause to avoid hammering the scheduler
  sleep 0.5
done

echo "All jobs submitted."