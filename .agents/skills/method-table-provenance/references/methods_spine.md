------------------------------------------------------------------------

# Methods Spine Workflow

Use this workflow to turn a classified provenance table into a concise, reader-facing Methods outline.

This workflow does **not** write manuscript prose.

## Inputs

Required provenance table:

``` text
id | parent | what | why | comment
```

Required classification table:

``` text
id | node_kind | include_in_methods | notes
```

`include_in_methods` must be one of:

``` text
include_main
include_supplement
include_code_availability
exclude
```

## Output

A Methods spine consisting of one outline packet per bucket, plus an overseer reconciliation report.

## Buckets

Use these fixed buckets:

``` text
study_design_and_scope
measurement_generation
data_assembly_and_qc
derived_variables_and_summaries
model_or_algorithm_specification
fitting_inference_and_model_selection
simulation_or_prediction
sensitivity_and_robustness
code_data_and_reproducibility
```

## Roles

### Overseer

The Overseer assigns one bucket to each worker, checks coverage, reconciles overlaps, and produces the final Methods spine.

### Bucket worker

Each worker handles one bucket. The worker searches the classified provenance table for relevant nodes and produces a terse, node-supported outline.

Workers do **not** write prose.

## Overseer Workflow

1.  Validate that every provenance-table `id` appears exactly once in the classification table.

2.  Assign one bucket per worker subagent. Escalate to user if subagents are unavailable.

3.  Give each worker:

    -   the provenance table;
    -   the classification table;
    -   the full bucket list;
    -   their assigned bucket.

4.  Tell workers to produce bullet outlines only. Every bullet must cite supporting node IDs.

5.  Collect bucket packets.

6.  Check that every `include_main` node is represented in at least one bucket or explicitly accounted for.

7.  Identify nodes appearing in multiple buckets.

8.  Reconcile overlaps:

    -   Accept shared nodes when they genuinely support multiple buckets.
    -   Flag conflicts when two buckets appear to explain the same method, assign different meaning to the same node, or disagree on primary versus supporting status.

9.  Produce the final Methods spine and reconciliation report.

## Bucket Worker Workflow

1.  Read the assigned bucket definition.

2.  Search the classified provenance table for relevant nodes.

3.  Prioritize nodes marked:

    -   `include_main`
    -   `include_supplement`
    -   `include_code_availability`

4.  Use `exclude` nodes only as context, not as outline support, unless the node appears to hide a consequential method detail.

5.  Inspect referenced files/scripts/configs only briefly, and only when needed to understand what the node contributes.

6.  Produce a terse outline. Each bullet should state one Methods-relevant point and cite node IDs.

7.  Return the packet using this format:

``` markdown
## Methods bucket: `<bucket_name>`

### Scope summary

<One to three sentences explaining what this bucket covers.>

### Bullet outline

- <Methods-relevant point.> Supported by: `node_id`, `node_id`.
- <Methods-relevant point.> Supported by: `node_id`.

### Node table

| node_id | role | include_in_methods | terse finding |
|---|---|---|---|
| `node_id` | primary | include_main | <Finding.> |
| `node_id` | supporting | include_supplement | <Finding.> |
| `node_id` | code/repro | include_code_availability | <Finding.> |
| `node_id` | context | exclude | <Only if needed.> |

### Open questions

- <Ambiguity, possible overlap, or item needing deeper inspection.>
```

## Bucket Definitions

### `study_design_and_scope`

Use for nodes defining what was studied, analysis populations, units, comparison structure, inclusion/exclusion logic, and primary versus secondary/sensitivity scope.

### `measurement_generation`

Use for nodes defining how raw observations were generated, acquired, annotated, segmented, classified, scored, assayed, or measured.

### `data_assembly_and_qc`

Use for nodes defining how observations and metadata were combined into analysis-ready data, including filtering, calibration application, missing-data handling, censoring treatment, and QC rules affecting retained data.

### `derived_variables_and_summaries`

Use for nodes defining variables, transformed quantities, summary metrics, endpoints, contrasts, ratios, rates, smoothing rules, aggregation rules, or descriptive features.

### `model_or_algorithm_specification`

Use for nodes defining model or algorithm structure: equations, states, parameters, likelihoods, observation models, priors, constraints, model families, classifiers, decision rules, or null/reference models.

### `fitting_inference_and_model_selection`

Use for nodes defining optimization, sampling, training, inference, parameter estimation, convergence checks, diagnostic thresholds, model comparison, or model selection.

### `simulation_or_prediction`

Use for nodes defining simulations, fitted predictions, counterfactuals, scenario tests, initial conditions, simulated schedules/inputs, objectives, endpoint scoring, or simulated output summaries.

### `sensitivity_and_robustness`

Use for nodes defining alternative scopes, sensitivity analyses, leave-one-out analyses, held-out tests, threshold sensitivity, model-family sensitivity, initial-condition sensitivity, or robustness comparisons.

### `code_data_and_reproducibility`

Use for nodes identifying scripts, configs, manifests, repository assets, saved outputs, software references, run specifications, data availability assets, or code availability assets.

## Reconciliation Report Format

``` markdown
# Methods Spine Reconciliation

## Bucket index

| bucket | status | notes |
|---|---|---|

## Include-main coverage

| node_id | bucket | status | notes |
|---|---|---|---|

## Duplicate node assignments

| node_id | buckets | decision |
|---|---|---|

## Unresolved issues

- <Issue requiring human review or deeper inspection.>

## Final recommended spine

1. <Recommended section/bucket order>
2. <Recommended section/bucket order>
```

## Rules

-   Do not draft manuscript prose.
-   Do not narrate provenance mechanics.
-   Do not use graph atoms as section headings.
-   Do not bundle mechanically by `node_kind`.
-   Do not delete excluded nodes from provenance; just keep them out of the outline unless they matter as context.
-   Every outline bullet must cite node IDs.
-   If a required Methods point is missing from the graph, list it under open questions.
