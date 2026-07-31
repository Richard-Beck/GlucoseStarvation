---
name: analysis
description: "Scope, execute, validate, and document analysis work using coherent analysis packages, controlled compute, provenance, and durable outputs. Use only when explicitly invoked by the user, or when an agent is assigned analysis work or analytical assessment within a multi-agent workflow. Do not invoke merely because a standalone task involves analysis, data processing, coding, debugging, or manuscript work."
---


# Analysis

## Purpose

Use this skill to perform analysis work in a structured, reproducible manner.

Inputs may include direct user instructions, a delegated analysis request, or a
feedback packet requiring analytical assessment. Infer candidate analytical
questions, scope, outputs, and constraints from the request and relevant project
context.

This skill owns analysis execution and analysis-scoped conclusions.

## Workflow

### 1. Scope the analysis

* Identify candidate analytical questions and the decisions or outputs they
  would support.
* Inspect enough relevant data, provenance, code, and project context to
  determine which candidate questions are real, answerable, and in scope.
* Identify known inputs, missing inputs, constraints, approval state, and likely
  validation requirements.
* Prefer canonical inputs and maintained implementations.
* Use subagents for broad discovery when this would preserve the coordinating
  agent’s context or materially accelerate inspection.

This step ends when candidate questions have been confirmed, combined, rejected,
deferred, or marked as blocked, and the remaining work can be represented as one
or more analysis packages.

### 2. Define analysis packages

Organize the work into the smallest number of coherent packages that can be
executed and validated without unnecessary coordination or context burden.

Split work when separation would materially enable:

* independent or parallel progress;
* context isolation for implementation-heavy work;
* different compute or approval boundaries;
* separate completion, failure, or deferral;
* independent validation.

Each package should identify:

* the analytical question or decision;
* known inputs or the required discovery scope;
* the expected output or evidence;
* dependencies on other packages;
* its completion, blocker, or approval-gate condition.

Do not create separate packages merely for inspection, implementation,
execution, validation, or provenance recording. These are stages within a
package.

Repeated execution of the same analysis across datasets, parameters, or
conditions is normally one sharded package.

### 3. Advance each package

Within each package, proceed with safe and reversible work, including:

* input and provenance checks;
* feasibility assessment;
* small summaries;
* lightweight implementation;
* smoke tests;
* reduced-data runs;
* validation of existing outputs.

Continue until the package:

* completes;
* becomes blocked;
* is shown to require no analytical change; or
* reaches a major compute or irreversible-write boundary.

Before expensive compute, broad reprocessing, major model fitting, or
irreversible writes, create a brief analysis plan and obtain approval unless
that specific work has already been approved.

Follow applicable `AGENTS.md` and local project instructions for execution,
scheduling, resources, runners, and storage.

### 4. Execute approved work

* Reuse maintained code where practical.
* Keep entrypoints thin and place reusable logic in appropriate shared modules.
* Smoke-test a small unit before broad execution when practical.
* Record commands, parameters, code state, inputs, outputs, and execution
  identifiers at a level appropriate to the analysis.

### 5. Validate

* Check input integrity, output integrity, and scientific plausibility.
* Compare against baselines, simulations, null expectations, or independently
  recomputed summaries where useful.
* Use independent validation when it materially improves confidence.
* Do not promote failed or inadequately validated outputs for reuse.

### 6. Synthesize

* Answer the original analytical questions.
* State what was completed, what outputs were produced, what validation was
  performed, and whether outputs are suitable for reuse.
* Record important assumptions, unresolved failures, deferred work, and
  blockers.
* Return analysis-scoped evidence without claiming completion outside the
  analysis scope.

## Subagent Delegation

Subagents may be used at any stage when they materially improve speed or preserve
the coordinating agent’s context.

Good candidates include:

* broad repository or input discovery;
* provenance tracing;
* bounded implementation work;
* repetitive data inspection or processing;
* debugging and log-heavy execution;
* independent validation.
* Discrete analysis packages

The coordinating agent should retain:

* the original request;
* confirmed analytical questions;
* package boundaries and dependencies;
* scientific constraints;
* approval state;
* acceptance criteria;
* final synthesis.

Do not delegate solely because more than one package exists.

When delegating, provide:

* a bounded task;
* known starting points or a defined search scope;
* relevant constraints;
* permitted write scope;
* expected return and evidence;
* required validation;
* stop conditions.

Require a concise return containing:

* findings or work completed;
* files inspected or changed;
* commands or jobs run;
* outputs produced;
* validation status;
* assumptions and blockers.

Subagents must not revert or overwrite unrelated work.

## Output Layout

For substantial analysis work, use one coordination root:

```text
major_analyses/<YYYYMMDD_slug>/
  analysis_plan.md
  contract.json
  run_manifest.json
  notes.md
  validation_report.md
  logs/
```

Create only files and optional subdirectories that are useful. Place durable
scientific outputs according to project conventions or caller instructions.

## Analysis Plan

Use this minimal structure when an approval plan is required:

```markdown
# <Analysis Title>

Status: <pending approval | approved | in progress | complete | blocked>

## Question Or Instruction
<What this analysis will decide or produce.>

## Work Packages
- <Package> - <status>

## Direct Inputs
- <Canonical input> - <why it is appropriate>

## Derived Inputs
- <Derived input or "none"> - <provenance status>

## Implementation Plan
- Reused code:
- New reusable helpers:
- Thin entrypoints:
- Output root:

## Execution Plan
- Runtime or resource class:
- Execution approach:
- Parallelization or delegation:
- Smoke test:

## Validation Plan
- Input checks:
- Execution checks:
- Scientific checks:
- Reproducibility checks:

## Output Contract
- Tables or structured outputs:
- Figures or review artifacts:
- Model or processed-data outputs:
- Manifests and logs:
- Return evidence:

## Approval
Approve, revise, defer, or request a feasibility-only pass.
```

## Contract Fields

Write `contract.json` when the analysis will produce durable or reusable outputs:

```json
{
  "contract_type": "analysis",
  "schema_version": "0.3",
  "analysis_root": "major_analyses/<YYYYMMDD_slug>",
  "approval_status": "approved|pending|revised|deferred|not_required",
  "question_or_instruction": "",
  "work_packages": [],
  "blocked_major_compute_packages": [],
  "direct_inputs": [],
  "derived_inputs": [],
  "constraints": [],
  "entrypoints": [],
  "reusable_helpers": [],
  "output_roots": [],
  "execution_plan": {
    "approach": "",
    "parallelization": "",
    "smoke_test": ""
  },
  "validation_outputs": [],
  "return_evidence": []
}
```

Omit irrelevant fields rather than inventing content.

## Validation Expectations

Match validation to the risk and intended reuse of the analysis. Substantial
analyses should normally include:

* input dimensions, keys, factors, conditions, and missingness;
* a smoke execution on a reduced unit;
* output existence, schema, and loadability checks;
* an appropriate scientific sanity check;
* recorded commands, parameters, code state, logs, and execution identifiers.

If validation fails, do not promote the outputs for reuse. Record the blocker and
the smallest useful next action.

## Finish Criteria

The analysis is complete when:

* each analytical question is answered, deferred, or clearly blocked;
* declared outputs exist or their absence is explained;
* validation is adequate for the intended use;
* important inputs, parameters, commands, code state, and outputs are recorded;
* assumptions and scientific decisions are visible;
* the caller receives a concise synthesis with paths to supporting artifacts.

A supported no-change or not-applicable conclusion is a valid analysis result.

