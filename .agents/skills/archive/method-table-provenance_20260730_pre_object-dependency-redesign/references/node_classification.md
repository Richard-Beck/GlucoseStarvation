Node Classification Workflow.

Use this optional workflow to classify each provenance-table node by both:

1) its structural role in the provenance graph; and
2) whether it creates a Methods obligation.

This is one integrated node-classification workflow. Do not produce a separate
secondary inclusion table.

Required inputs: A single provenance table with structure
id | parent | what | why | comment

Required outputs: A single output table with structure
id | node_kind | include_in_methods | reason_code | notes

Output validation:
1) All input ids are represented exactly once in the output table.
2) `node_kind` exactly matches one entry in the Node Kind table below.
3) `include_in_methods` exactly matches one entry in the Methods Inclusion table below.
4) `reason_code` exactly matches one entry in the Reason Code table below.
5) `notes` may be empty, but should be nonempty when:
   - `include_in_methods = include_supplement`
   - `include_in_methods = include_code_availability`
   - `reason_code = deprecated_or_legacy`
   - `reason_code = terminal_boundary_or_gap`
   - the row may need source/code inspection before drafting.

## Node Kind table

node_kind                   Meaning
━━━━━━━━━━━━━━━━━━━━━━━━━━  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
display_node                Figure/panel/render/layout/visual encoding nodes.
──────────────────────────  ─────────────────────────────────────────────────────────────────────
source_node                 Raw, manual, received, or external inputs.
──────────────────────────  ─────────────────────────────────────────────────────────────────────
data_artifact               Any generated data/result/model-output artifact.
──────────────────────────  ─────────────────────────────────────────────────────────────────────
computation_step            Action/run/function step that transforms inputs into outputs.
──────────────────────────  ─────────────────────────────────────────────────────────────────────
code_or_config_reference    Fixed code, model definition, config, model list, mapping, or spec.
──────────────────────────  ─────────────────────────────────────────────────────────────────────
provenance_marker           Synthetic terminal/gap/boundary marker; not a real analysis object.

## Methods Inclusion table

include_in_methods          Meaning
━━━━━━━━━━━━━━━━━━━━━━━━━━  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
include_main                The node locally defines something a reader needs to understand the
                            primary methods.
──────────────────────────  ─────────────────────────────────────────────────────────────────────
include_supplement          The node is legitimate methodological detail, but too granular,
                            exhaustive, diagnostic, or secondary for the main Methods.
──────────────────────────  ─────────────────────────────────────────────────────────────────────
include_code_availability   The node matters for reproducibility, but belongs in repository,
                            data, code, or workflow availability language rather than Methods prose.
──────────────────────────  ─────────────────────────────────────────────────────────────────────
exclude                     The node is not Methods-facing.

Do not use `manual_review` or any equivalent limbo category. Handle ambiguity in
`notes` while still assigning the best local classification.

## Reason Code table

reason_code
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
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

Keep `reason_code` local and functional. Do not use relation-dependent reasons
such as "duplicate of another included node"; duplicate collapse, module
reduction, and section placement belong to later graph-reduction or drafting
workflows.

Do not classify a node as "interpretive" or "rhetorical." A row may contain
interpretive language, but this workflow classifies what the node is or does in
the graph and whether it creates a Methods obligation.

## Inclusion guidance

Use `include_main` when the node locally defines or materially specifies:

- study design or analysis population
- experimental, sample, or computational units
- data acquisition or measurement procedure
- preprocessing that changes the meaning of observations
- calibration or standardization
- inclusion/exclusion criteria
- variable, feature, endpoint, contrast, or summary metric
- statistical, mathematical, mechanistic, or algorithmic model
- estimation, fitting, optimization, or inference procedure
- model comparison or model selection rule
- validation, QC, or diagnostic criterion that affects interpretation
- uncertainty quantification
- simulation or prediction protocol
- sensitivity/robustness analysis design
- consequential threshold, filter, or assumption

Use `include_supplement` for legitimate methodological detail that is too
granular for the main Methods, such as:

- detailed parameter/configuration tables
- extensive QC inventories
- non-primary sensitivity branches
- detailed implementation settings
- exhaustive model lists
- long diagnostic summaries
- secondary validation details
- reproducibility-level workflow detail that is more than code availability but
  less than main Methods

Use `include_code_availability` when the node mainly identifies reproducibility
assets, such as:

- script paths
- configuration files
- manifest files
- package inventories
- run specifications
- exact saved output locations
- software environment files
- cached reproducibility artifacts

This does not mean the node is scientifically irrelevant. A model-definition
file can be `include_main` if it is the place where the model is specified. A
script that merely executes an already described method is usually
`include_code_availability`.

Use `exclude` when the node is not Methods-facing, such as:

- figure rendering
- layout/composition
- visual styling
- panel integration
- cache materialization
- local file conversion
- serialization/deserialization
- mechanical reshaping
- mechanical joining/indexing that does not define a variable or population
- workflow status markers
- round/version markers
- audit markers
- provenance-only terminal markers
- deprecated/legacy artifacts not used to define reported methods

## Classification rule tree

1) Does this node locally define what was studied, measured, included, excluded,
   modeled, estimated, simulated, validated, or summarized?
   -> `include_main`, unless the detail is clearly too granular.

2) Does this node locally define a legitimate methodological detail that is
   secondary, exhaustive, diagnostic, or too detailed for the main Methods?
   -> `include_supplement`.

3) Does this node mainly identify code, configuration, run manifests, saved
   outputs, or reproducibility assets?
   -> `include_code_availability`.

4) Is this node only display, rendering, implementation plumbing, workflow
   provenance, status/versioning, or an unused legacy artifact?
   -> `exclude`.

Assign `node_kind` independently from `include_in_methods`. For example, a
`code_or_config_reference` can be `include_main` when it defines a model, and a
`computation_step` can be `exclude` when it only performs implementation
plumbing.

## Notes

Use `notes` for concise local nuance, not as an uncontrolled category system.
Good notes include:

- Read source/config if this node is used for drafting; it appears to contain
  the model specification rather than merely executing it.
- Exclude from main Methods; this is a rendering artifact, but retain in
  code/provenance documentation.
- Include only if this threshold affects retained observations or reported
  endpoints.
- This row appears to be a terminal provenance boundary; do not convert into
  Methods prose.
- Local classification only. Whether this should be collapsed with related nodes
  is a later graph-reduction decision.

## Workflow

1) Assess rowcount of input table.
2) If rowcount <= 40, you are a classification worker. If rowcount > 40, you are
   a classification overseer.
3) Proceed to role-specific workflow.

## Classification Worker Workflow

1) Read each row thoroughly.
2) Assign `node_kind`, `include_in_methods`, and `reason_code` using only the
   local row and any necessary source inspection.
3) Add `notes` when validation requires it or when a local caveat matters for
   drafting.
4) Validate your output table.
5) Follow any instructions about where to save your output.

## Classification Overseer Workflow

1) Split input table rows into batches of <= 40 rows.
2) Determine output location. Default to a tmp folder near the input for
   intermediate batch tables, and a single file near the input for final merged
   classification table.
3) Assign each batch to a subagent. If subagents are unavailable, escalate to the
   user for permission or additional guidance.
4) Once agents are finished, combine each batch into a single table.
5) Validate final table.
6) Delete tmp batch folder.
