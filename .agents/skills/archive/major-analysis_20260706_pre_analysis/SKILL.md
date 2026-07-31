---
name: major-analysis
description: "Plan, approve, execute, and validate significant glucose-starvation analysis packages that are time-consuming, compute-expensive, scientifically consequential, or operationally complex. Use when Codex needs to run or design major analyses such as new model fits, NUTS/optimization batches, large simulations, image reprocessing, classifier training, substantial dataset rebuilding, or analysis packages that may affect manuscript claims, figures, or maintained workflows."
---

# Major Analysis

## Purpose

Use this skill for substantial analysis packages that need a decision gate before execution. The goal is to keep analyses reproducible, fast enough for the HPC environment, scientifically auditable, and integrated with the existing project rather than creating one-off repo bloat.

This skill does not replace manuscript figure drafting, polishing, or integration workflows. It produces validated analysis outputs that those downstream workflows may later consume.

## Operating Rule

Do not launch expensive jobs, long local runs, broad image processing, new model fitting, or irreversible dataset writes until the user has approved a brief analysis plan. Lightweight discovery, input inspection, feasibility checks, and smoke tests are allowed before approval when they reduce ambiguity.

If the user already supplied explicit approval for a specific plan, proceed with execution and record that approval in the analysis notes.

## Workflow

1. Orient to the project.
   - Read `docs/project_map.txt` before planning.
   - Identify the maintained workflow family, canonical input objects, relevant prior outputs, and downstream manuscript or figure dependencies.
   - Prefer maintained scripts and helpers over new code.

2. Create a brief human-readable plan and stop for approval.
   - Write the plan before implementation under `agent-dev/major_analyses/<YYYYMMDD_slug>/analysis_plan.md` unless the user specifies another root.
   - Keep it short enough for real review.
   - Include the scientific question, direct inputs, output contract, reusable-code plan, compute/SLURM plan, validation plan, expected runtime class, and approval status.
   - End with explicit user options: approve, revise scope, defer, or ask for a feasibility-only pass.

3. After approval, implement narrowly.
   - Reuse existing functions from `R/`, `src/`, maintained `scripts/`, and pipeline utilities wherever possible.
   - Put broadly reusable functions in an appropriate shared file such as `R/<domain>_utils.R` or another existing utility module, with a short comment only when the intent is not obvious.
   - Keep entrypoint scripts thin: parse arguments, call shared helpers, write declared outputs, and record provenance.
   - Keep package-specific glue inside the analysis root only when it is not reusable.

4. Deploy appropriately.
   - Run R scripts through `scripts/agentRrunner.sh`.
   - Run small jobs directly from the terminal.
   - Use SLURM for parallel or batch work, GPU jobs, large-memory jobs, or expected runtime over about 5 minutes.
   - Use arrays or independent shards for embarrassingly parallel work; combine results with a deterministic reducer.
   - Check current QOS limits before full submission with `sacctmgr show qos format=Name,Priority,MaxTRESPU,MaxJobsPU,MaxSubmitJobsPU,GrpTRES -P`.
   - Choose the smallest QOS/resource request that comfortably fits the job, and scale arrays to minimize wall time without exceeding per-user or group limits.
   - Smoke test one small unit before submitting a full batch.

5. Validate and summarize.
   - Run the validation plan promised in `analysis_plan.md`.
   - Write `validation_report.md` or `validation_report.json` under the analysis root.
   - Record commands, SLURM job ids, logs, package versions when relevant, git commit, important input checksums, and output checksums or manifests.
   - State which outputs are ready for downstream figure/manuscript work and which are not.

6. Decide whether to update `docs/project_map.txt`.
   - Update it when the analysis adds or changes a maintained entrypoint, canonical output path, reusable helper family, workflow branch, or manuscript-relevant result location.
   - Do not bloat it with one-off run details; point to the stable script and output root.

## Data Lineage

Prefer a simple active dataset path. A good major-analysis package depends on a canonical source such as `data/stan_ready_data.Rds`, `data/stan_ready_data_with_area.Rds`, maintained count tables, or an explicit, reproducible pipeline output.

Avoid depending on a long chain of derived objects, especially when provenance is unclear. If a derived object is necessary:

- Trace and record its direct builder, input files, parameters, commit, and checksum.
- Prefer rebuilding it from canonical inputs inside the package or adding a small maintained builder.
- Treat dubious provenance as a blocker or user decision, not as a silent dependency.
- Keep new derived objects near their workflow family, with a manifest explaining how to regenerate them.

## Output Layout

Use one analysis root for coordination artifacts:

```text
agent-dev/major_analyses/<YYYYMMDD_slug>/
  analysis_plan.md
  contract.json
  run_manifest.json
  notes.md
  validation_report.md
  logs/
  slurm/
```

Place durable scientific outputs in established project locations when they are meant to be reused:

- Tables and audit exports: `data/report_exports/<analysis_slug>/`
- Figures or review plots: `figures/<analysis_slug>/`
- Model run outputs: the maintained `data/runs/...` family when the analysis belongs to an existing model pipeline
- Job records: `slurm/runs/<analysis_slug>/` or the relevant existing SLURM run family

Avoid scattering outputs across ad hoc folders. If a location deviates from these conventions, explain why in `notes.md`.

## Plan Contents

Use this minimal structure for `analysis_plan.md`:

```markdown
# <Analysis Title>

Status: tentative pending user approval

## Question
<One short paragraph describing what this package will decide.>

## Direct Inputs
- <Canonical input path> - <why it is appropriate>
- <Any derived input> - <provenance status and rebuild plan>

## Implementation Plan
- Reused code:
- New reusable helpers:
- Thin entrypoints:
- Output root:

## Compute Plan
- Runtime class:
- Local vs SLURM:
- Parallelization:
- QOS/resource request:
- Smoke test:

## Validation Plan
- Input checks:
- Smoke/sanity checks:
- Scientific checks:
- Reproducibility/provenance checks:

## Output Contract
- Tables:
- Figures/review artifacts:
- Manifests/logs:
- Downstream consumers:

## Approval
Approve, revise, defer, or request feasibility-only review before execution.
```

## Contract Fields

Write `contract.json` before or during implementation when the package will produce reusable outputs. Keep it compact:

```json
{
  "contract_type": "major_analysis",
  "schema_version": "0.1",
  "analysis_root": "agent-dev/major_analyses/<YYYYMMDD_slug>",
  "approval_status": "approved|pending|revised|deferred",
  "question": "",
  "direct_inputs": [],
  "derived_inputs": [],
  "entrypoints": [],
  "reusable_helpers": [],
  "output_roots": [],
  "compute_plan": {
    "uses_slurm": false,
    "qos": "",
    "parallelization": "",
    "smoke_test": ""
  },
  "validation_outputs": [],
  "project_map_decision": ""
}
```

## Validation Expectations

Select validation checks that match the risk of the analysis. At minimum, every major analysis should include:

- Input integrity: expected rows, keys, dimensions, line/ploidy/glucose/time coverage, missingness, and checksum/provenance where relevant.
- Smoke execution: one small unit or reduced dataset before full deployment.
- Output integrity: declared files exist, are nonempty, have expected schemas, and can be loaded by downstream code.
- Scientific sanity: compare against an existing baseline, null model, prior maintained output, simulated expectation, or independently recomputed summary.
- Reproducibility: record commands, script paths, parameters, git state, SLURM ids, logs, seeds, and package/session information.

If validation fails, do not promote outputs to downstream figure or manuscript workflows. Record the blocker and the smallest next action.

## Finish Criteria

A major-analysis package is complete only when:

- The approved plan, contract, notes, run manifest, logs or job ids, and validation report are saved under the analysis root.
- Durable outputs are in their declared locations and have manifests or checksums.
- Validation passed or unresolved failures are clearly marked as blockers.
- Reusable code is documented by path and does not duplicate an existing helper.
- `docs/project_map.txt` was updated or the decision not to update it is recorded.
