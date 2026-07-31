# Node Classification

Classify each row in a structurally validated canonical provenance table by its
graph role and its local Methods obligation. Do not create a second inclusion
table.

## Contents

- [Inputs and output](#inputs-and-output)
- [Node kinds](#node-kinds)
- [Methods inclusion](#methods-inclusion)
- [Reason codes](#reason-codes)
- [Workflow](#workflow)

## Inputs and Output

Require the canonical table:

```text
id | parent | what | why | comment | lock_target | lock_selector | lock_kind | sha256 | hash_status
```

Require passing panel-leaf and lock verification before reusing an existing
classification.

Return exactly:

```text
id | node_kind | include_in_methods | reason_code | notes
```

Represent every input id exactly once. Use only the controlled values below.
Keep `notes` concise; populate it for supplementary/code-availability rows,
terminal gaps, deprecated inputs, and rows needing source inspection.

## Node Kinds

```text
display_node
source_node
data_artifact
computation_step
code_or_config_reference
provenance_marker
```

Use `display_node` for the current manuscript panel leaves. Rendering, layout,
composition, and package-assembly nodes should already have been removed from a
valid canonical graph.

Use `source_node` for raw, manual, received, or external inputs;
`data_artifact` for generated data, results, or model outputs;
`computation_step` for Methods-relevant transformations;
`code_or_config_reference` for fixed code, model definitions, configs, mappings,
or specifications; and `provenance_marker` only for explicit terminal gaps or
boundaries.

## Methods Inclusion

```text
include_main
include_supplement
include_code_availability
exclude
```

Use `include_main` when the row locally defines or materially specifies:

- study design, population, units, measurements, or inclusion rules;
- preprocessing, calibration, variables, endpoints, contrasts, or summaries;
- models, algorithms, estimation, inference, or model selection;
- QC, validation, uncertainty, simulation, prediction, sensitivity, or a
  consequential threshold or assumption.

Use `include_supplement` for legitimate but granular, exhaustive, secondary, or
diagnostic methodological detail.

Use `include_code_availability` when the row mainly identifies reproducibility
assets such as scripts, configs, manifests, software environments, run specs,
or saved output locations. A code/config row may instead be `include_main` when
it is the actual method specification.

Use `exclude` for current panel endpoints, implementation plumbing, cache or
serialization mechanics, workflow markers, terminal gap markers, and unused
legacy inputs. A structurally valid input should not contain rendering or
historical-overlay rows.

Do not use a `manual_review` limbo category. Assign the best local value and
record uncertainty in `notes`.

## Reason Codes

Use exactly one:

```text
defines_design_or_population
defines_measurement
defines_preprocessing_or_calibration
defines_inclusion_exclusion
defines_variable_or_endpoint
defines_model_or_algorithm
defines_estimation_or_inference
defines_model_selection_or_comparison
defines_qc_or_validation
defines_uncertainty
defines_simulation_or_prediction
defines_sensitivity_or_robustness
defines_threshold_or_assumption
needed_for_reproducibility
display_only
implementation_only
workflow_provenance_only
deprecated_or_legacy
terminal_boundary_or_gap
```

Keep the reason local and functional. Duplicate collapse and section placement
belong to construction and spine workflows, not classification. Do not classify
rows as interpretive or rhetorical.

## Workflow

1. Read every row and inspect referenced sources only when local meaning is
   unclear.
2. Assign `node_kind`, `include_in_methods`, and `reason_code` independently.
3. Add a note only when it changes downstream handling.
4. Validate controlled values and exact one-to-one id coverage.

For more than 40 rows, workers may classify disjoint batches of at most 40 rows.
Merge batches without changing ids or creating multiple classifications per row,
then validate the complete table.
