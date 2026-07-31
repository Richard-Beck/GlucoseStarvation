#!/usr/bin/env bash
set -euo pipefail

ACTION=${1:?usage: run.sh ACTION CONFIG [phase]}
CONFIG=${2:?release config required}
shift 2

case "$ACTION" in
  prepare)
    scripts/agentRrunner.sh scripts/parameter_estimation/build_counts.R "$CONFIG"
    scripts/agentRrunner.sh scripts/parameter_estimation/build_stan_data.R "$CONFIG"
    scripts/agentRrunner.sh scripts/parameter_estimation/materialize_datasets.R "$CONFIG"
    scripts/agentRrunner.sh scripts/parameter_estimation/build_fit_plan.R "$CONFIG"
    scripts/agentRrunner.sh scripts/parameter_estimation/validate_release.R "$CONFIG" prepare
    ;;
  submit-optim|submit-nuts)
    release_root=$(scripts/agentRrunner.sh scripts/parameter_estimation/release_path.R "$CONFIG")
    dry_run=0
    dataset_filter=""
    deferred_optim_jobs=""
    while (($#)); do
      case "$1" in
        --dry-run)
          dry_run=1
          shift
          ;;
        --datasets)
          [[ "$ACTION" == "submit-optim" ]] || { echo "--datasets is only supported for submit-optim" >&2; exit 1; }
          [[ $# -ge 2 ]] || { echo "--datasets requires a file" >&2; exit 1; }
          dataset_filter=$2
          shift 2
          ;;
        --defer-after-optim)
          [[ "$ACTION" == "submit-nuts" ]] || { echo "--defer-after-optim is only supported for submit-nuts" >&2; exit 1; }
          [[ $# -ge 2 ]] || { echo "--defer-after-optim requires a wave job_ids.tsv" >&2; exit 1; }
          deferred_optim_jobs=$2
          shift 2
          ;;
        *)
          echo "Unknown submission option: $1" >&2
          exit 1
          ;;
      esac
    done
    if [[ "$ACTION" == "submit-optim" ]]; then
      scripts/parameter_estimation/submit_optim.sh "$release_root" "$dry_run" "$dataset_filter"
    else
      scripts/parameter_estimation/submit_nuts.sh "$release_root" "$dry_run" "$deferred_optim_jobs"
    fi
    ;;
  status)
    scripts/agentRrunner.sh scripts/parameter_estimation/status.R "$CONFIG"
    ;;
  derive-plan)
    scripts/agentRrunner.sh scripts/parameter_estimation/make_derived_plan.R "$CONFIG"
    ;;
  derive-optim)
    scripts/agentRrunner.sh scripts/parameter_estimation/build_optimization_assessment.R "$CONFIG"
    ;;
  submit-derived)
    dry_run=0
    while (($#)); do
      case "$1" in
        --dry-run) dry_run=1; shift ;;
        *) echo "Unknown submission option: $1" >&2; exit 1 ;;
      esac
    done
    scripts/agentRrunner.sh scripts/parameter_estimation/make_derived_plan.R "$CONFIG"
    release_root=$(scripts/agentRrunner.sh scripts/parameter_estimation/release_path.R "$CONFIG")
    scripts/parameter_estimation/submit_derived.sh "$release_root" "$CONFIG" "$dry_run"
    ;;
  derived-status)
    scripts/agentRrunner.sh scripts/parameter_estimation/derived_status.R "$CONFIG"
    ;;
  strategy-plan)
    while (($#)); do
      case "$1" in
        --strategy-config)
          [[ $# -ge 2 ]] || { echo "--strategy-config requires a file" >&2; exit 1; }
          export GPATH_STRATEGY_CONFIG=$2
          shift 2
          ;;
        *) echo "Unknown strategy-plan option: $1" >&2; exit 1 ;;
      esac
    done
    scripts/agentRrunner.sh scripts/parameter_estimation/strategy_simulations/make_plan.R "$CONFIG"
    ;;
  submit-strategies)
    dry_run=0
    strategy_config=""
    while (($#)); do
      case "$1" in
        --dry-run) dry_run=1; shift ;;
        --strategy-config)
          [[ $# -ge 2 ]] || { echo "--strategy-config requires a file" >&2; exit 1; }
          strategy_config=$2
          export GPATH_STRATEGY_CONFIG=$strategy_config
          shift 2
          ;;
        *) echo "Unknown strategy submission option: $1" >&2; exit 1 ;;
      esac
    done
    scripts/agentRrunner.sh scripts/parameter_estimation/strategy_simulations/make_plan.R "$CONFIG"
    release_root=$(scripts/agentRrunner.sh scripts/parameter_estimation/release_path.R "$CONFIG")
    scripts/parameter_estimation/strategy_simulations/submit.sh \
      "$release_root" "$CONFIG" "$dry_run" \
      "${strategy_config:-$(dirname "$CONFIG")/strategy_simulations.json}"
    ;;
  strategy-status)
    while (($#)); do
      case "$1" in
        --strategy-config)
          [[ $# -ge 2 ]] || { echo "--strategy-config requires a file" >&2; exit 1; }
          export GPATH_STRATEGY_CONFIG=$2
          shift 2
          ;;
        *) echo "Unknown strategy-status option: $1" >&2; exit 1 ;;
      esac
    done
    scripts/agentRrunner.sh scripts/parameter_estimation/strategy_simulations/status.R "$CONFIG"
    ;;
  validate-strategies)
    while (($#)); do
      case "$1" in
        --strategy-config)
          [[ $# -ge 2 ]] || { echo "--strategy-config requires a file" >&2; exit 1; }
          export GPATH_STRATEGY_CONFIG=$2
          shift 2
          ;;
        *) echo "Unknown validate-strategies option: $1" >&2; exit 1 ;;
      esac
    done
    scripts/agentRrunner.sh scripts/parameter_estimation/strategy_simulations/validate.R "$CONFIG"
    ;;
  reduce-strategy-endpoints)
    scripts/agentRrunner.sh scripts/parameter_estimation/strategy_simulations/reduce_endpoints.R "$CONFIG" "$@"
    ;;
  validate-derived)
    scripts/agentRrunner.sh scripts/parameter_estimation/validate_derived.R "$CONFIG"
    ;;
  validate)
    scripts/agentRrunner.sh scripts/parameter_estimation/validate_release.R "$CONFIG" "${1:-prepare}"
    ;;
  *)
    echo "Unknown action: $ACTION" >&2
    exit 1
    ;;
esac
