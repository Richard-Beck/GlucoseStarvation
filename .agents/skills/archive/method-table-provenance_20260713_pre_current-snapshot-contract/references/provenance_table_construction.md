# Provenance Table Construction Workflow

Use this optional workflow when the user asks to build, revise, or merge
method-table provenance rows for manuscript figures or panels.

The output is a single construction table, not polished Methods prose. The
table must include the five semantic provenance columns and the lock columns
below:

```text
id | parent | what | why | comment | lock_target | lock_selector | lock_kind | sha256 | hash_status
```

`lock_target` is the canonical repo-relative path whose bytes were hashed for
that row. For durable object rows, this should usually be the object itself. For
in-memory/action rows, this may be the most informative directly adjacent
durable upstream or downstream artifact. Use `NA` only when no useful hashable
target exists.

`lock_selector` records any row, panel, object, or fragment selector such as
`#panel_b` or `::case_e`. The selector is metadata only; the SHA256 should cover
the full `lock_target` file. Use `NA` when no selector is needed.

`lock_kind` should be one of `file`, `code`, `proxy_file`,
`representative_file`, `terminal`, `external`, or `unresolved`.

`hash_status` should state how the digest relates to the row, using values such
as `computed_self`, `computed_downstream_proxy`, `computed_upstream_proxy`,
`computed_code_proxy`, `computed_representative`, `not_applicable`, `missing`,
`ambiguous`, or `unresolved`.

## Inputs

Required:

- A set of manuscript figures or panels with enough metadata to trace each panel
  back through its generating artifacts and transformations.

This metadata may appear in any equivalent figure package. Common concrete
sources include:

- figure or panel manifest rows;
- package-level provenance tables;
- final panel or figure image paths;
- legends or semantic notes;
- package inventories;
- generating scripts, configs, commands, saved tables, model outputs, or other
  upstream artifacts.

Optional:

- User instructions, feedback, corrections, or scope constraints.
- Existing per-panel or merged provenance tables.

Treat named files such as `figure_set_manifest.csv`, `polishing/provenance.csv`,
`integrated_figure_legends.md`, package `legend.md`, and semantic interpretation
markdown as implementation examples, not universal requirements.

## Roles

- Overseer: assigned multiple panels, or assigned reconciliation over multiple
  per-panel chains.
- Tracer: assigned one manuscript panel to trace into a method-table provenance
  chain.
- Merger: assigned exactly two existing provenance chains to reconcile into an
  additional merged chain.

## Shared Guidelines

Use `id` for the generated object or action node represented by the row. Durable
objects should use their stable filename or a concise normalized object id.
Prefer canonical repo-relative paths for durable files when practical. If a
package-local shorthand path is used, make the package root explicit in
`comment`. In-memory transformations may use a concise synthetic id such as
`F2b_plot_table`.

Use `parent` for the immediate upstream object or action used to generate `id`.
Do not use `parent` as a grouping label or section name.

Allow multiple parents only when `id` is produced by a single merge, model,
plot, render, or integration step that directly consumes multiple immediate
upstream streams. When using multiple parents, make `what` name the combining
action and use `comment` to state each parent's role. Do not use multiple
parents to skip an intermediate action or in-memory object.

Keep `parent` parseable: use backticked parent ids or paths separated by
semicolons, and do not mix prose into the parent cell.

A row should represent one durable object, named in-memory object, or named
action boundary that Methods writing may need to distinguish. Do not combine
data construction, plot-table derivation, plot rendering, and figure integration
in one row when those steps have different inputs, scripts, assumptions, or
Methods relevance.

Prefer row boundaries at saved files, named functions, named plot/data objects,
model fits, joins/merges across streams, population-defining filters, and
terminal provenance boundaries.

Use `what` to state the generation or transformation performed. Keep it
operational and specific.

Use `why` to state the role of that generation step in the analysis chain. Keep
panel-specific rhetorical framing near the panel/plot rows; upstream data-object
rows should have broader generation-facing reasons.

Use `comment` for concrete references and constraints:

- generated object path(s);
- source object path(s);
- script/function/command references;
- row counts, filters, model scopes, or thresholds when verified;
- validation status or accepted exceptions if relevant to the method object;
- terminal reason when the trace reaches a raw, external, or unresolved boundary.

## Hash and Lock Guidelines

Generate SHA256 locks during provenance construction for every row where a
useful hashable target can be identified. Do not emit a separate lock table as
the primary output; record lock information in the provenance table columns.

For ordinary durable file rows, compute SHA256 over the file bytes and set:

- `lock_target` to the canonical repo-relative file path;
- `lock_selector` to `NA` unless the row represents a named fragment, panel, or
  selected table view;
- `lock_kind` to `file` or `code`;
- `sha256` to the computed digest;
- `hash_status` to `computed_self`.

For rows representing a selected table view, panel, or file fragment, hash the
whole base file and record the selector in `lock_selector`. Do not compute
selector-specific hashes unless the user explicitly asks for that extra audit
detail.

For in-memory objects or action nodes, hash a nearby durable artifact that would
be expected to change if the in-memory object or action changed. Prefer, in
order:

1. a directly generated downstream file or table;
2. the most specific durable upstream input consumed by the action;
3. the generating script, function file, or config, when no better data artifact
   is available.

Set `lock_kind` to `proxy_file` or `code`, and set `hash_status` to
`computed_downstream_proxy`, `computed_upstream_proxy`, or
`computed_code_proxy`. Explain the proxy choice in `comment` when it is not
obvious from the row.

For rows representing many run files, do not hash a manifest or spec as a
substitute for file contents unless the row is specifically about that manifest
or spec. Prefer a representative content file that is most likely to reveal
relevant drift for the row, set `lock_kind = representative_file`, set
`hash_status = computed_representative`, and state in `comment` why that file was
chosen and what larger file set or run scope it represents.

If no single representative file is honest enough for a many-file row, either
split the provenance into more specific durable-file rows or set `sha256 = NA`
and `hash_status = unresolved` or `ambiguous`, with a short explanation in
`comment`.

For directory bundles, do not compute recursive directory hashes by default.
Choose a representative durable content file, split the row into more specific
file rows, or mark the lock unresolved. Directory paths may still appear in
`comment` for orientation.

For terminal, external, missing, or unresolved boundaries, use `sha256 = NA` and
set `hash_status` to `not_applicable`, `external`, `missing`, `ambiguous`, or
`unresolved` as appropriate. A terminal row must still explain why the trace
stops.

If a known checksum already appears in upstream metadata, verify it against the
local artifact when practical. If the checksum is accepted from metadata rather
than recomputed, set `hash_status = metadata_checksum` and cite the metadata
source in `comment`.

Prefer structured metadata and scripts over prose summaries when deciding
immediate parent relationships.

For a single panel, return a compact Markdown table with the required semantic
and lock columns. For multiple panels, panel-scoped worker outputs may be used
as temporary intermediate artifacts, but the final deliverable must be one
integrated provenance table with the semantic columns and lock columns. Add
short notes only for scope issues that affect Methods writing, such as
broader-than-used manifest inputs, unresolved reconciliation ambiguity, or
semantic caveats that prevent overclaiming.

Do not write manuscript prose unless the user explicitly asks for prose after
reviewing the table.

## Overseer Workflow

For multiple requested panels, delegate exactly one panel per tracer subagent. If
a callable subagent or multi-agent facility is unavailable, stop before producing
the multi-panel table and ask the user whether to enable subagents or narrow the
request to a single panel.

Each tracer must return rows only for its assigned panel and should not rewrite
rows for another panel.

For multi-panel tracing, preserve each tracer's returned row set as an
intermediate panel-scoped section or file until reconciliation. The overseer may
add a short index or status summary, but must not compact, coalesce, or rewrite
panel row sets before the merger workflow.

When multi-panel tracing is complete, run the merger workflow:

1. Treat the per-panel tables as fixed inputs; preserve them unmodified and write
   merged chains as additional outputs.
2. Use `scripts/score_chain_overlap.py` to select high-overlap chain pairs to
   merge.
3. Run pairwise merge tasks in parallel by subagent delegation, one selected pair
   per merger subagent. Default to 8 parallel disjoint merges per batch unless
   the user specifies otherwise.
4. Give each merger exactly two chain tables, plus optional one-line overlap
   hints.
5. After each merge round, recompute overlap using merged tables plus any
   retained tables not modified in the previous round.
6. Repeat until a single integrated table exists, preserving and reconciling the
   lock columns as part of the table.

If merger subagents are unavailable, stop before producing reconciled
multi-chain output and ask the user whether to enable subagents or narrow the
request to one pair.

## Tracer Workflow

For a single requested panel, perform the trace directly.

Read sources in this order when they exist in the supplied figure package:

1. Panel identity metadata: exact figure/panel id, displayed content, source
   package, and row role.
2. Package or panel provenance metadata: panel image path, generating script or
   command, direct data inputs, context inputs, and final composite output.
3. Package inventory, if available.
4. Legend or semantic notes for manuscript-facing terminology, visible object
   identity, source clarification, caveats, and provenance hints.
5. Generating scripts and upstream scripts only as needed to resolve exact
   transformations, inputs, or terminal provenance.

Trace in this order:

1. Establish the leaf.
   - Identify the exact manuscript-visible panel and its displayed content.
   - Create a leaf row only after identifying the immediate plotted
     table/object/action that generated the displayed panel.
2. Find the package-level panel row or equivalent trace anchor.
   - Confirm subpanel image path, generating script or command, direct data
     inputs, context inputs, and final composite output.
   - If figure-level inputs are broader than the panel actually uses, narrow the
     trace using provenance metadata, plotting scripts, and semantic notes.
3. Trace upstream one edge at a time.
   - For each row, ask: "What immediate object or action generated this id?"
   - Read scripts or structured metadata to verify the edge before writing it
     when the parent is not obvious.
   - Do not skip material action nodes. If a function or named in-memory object
     filters, summarizes, joins, models, maps aesthetics, computes displayed
     quantities, or renders a file, represent that action/object as a row before
     tracing to its input table(s).
   - Continue from saved tables to their generating scripts, then to input
     tables/runs, until the trace reaches a raw, external, or unresolved
     boundary.
4. Terminal rows.
   - Do not stop at a vague "reporting boundary."
   - A terminal row must say why it is terminal, such as raw image groups,
     experimental/acquisition protocol outside current file provenance, an
     external fixed input, or an unresolved provenance gap.
   - If provenance is missing or ambiguous, write a terminal gap row rather than
     inventing a parent.
5. Write the table.
   - Order rows from leaf to upstream root unless the user asks otherwise.
   - Return one table for the assigned panel.
   - Keep paths concise but specific enough to locate the object. Prefer
     canonical repo-relative paths for lock targets.
   - Use backticks for ids, paths, functions, and field names.
   - Compute SHA256 values for directly hashable lock targets before returning
     the table.
   - For in-memory/action rows, choose a useful upstream or downstream proxy
     artifact to hash whenever possible.
   - Use explicit `hash_status` values for rows that cannot be hashed, rather
     than leaving lock columns blank.

## Merger Workflow

When assigned exactly two existing provenance chains, return one additional
merged table with the semantic provenance columns and lock columns.

The merged table must preserve meaningful provenance from both inputs, reconcile
duplicate nodes where justified, and keep distinct nodes distinct when they
represent different objects or actions.

Guidelines:

- Preserve immediate-parent semantics.
- Merge rows only when they describe the same generated object or action node.
- Prefer the most specific verified path/object id and the most complete
  upstream trace.
- When duplicate rows are already determined to describe the same durable
  artifact, prefer a canonical resolvable path over package-local shorthand.
  Preserve useful aliases in `comment`.
- For in-memory or action nodes, keep the clearest normalized object id, but
  preserve the most useful lock proxy available. Prefer downstream proxy hashes
  over upstream proxy hashes when both are equally specific.
- Preserve `lock_target`, `lock_selector`, `lock_kind`, `sha256`, and
  `hash_status` values during reconciliation. If two duplicate durable-artifact
  rows disagree on a checksum for the same lock target, keep the conflict
  visible and report it. If two duplicate in-memory/action rows use different
  proxy targets, choose the more specific proxy and mention the alternative in
  `comment` when it remains informative.
- Preserve panel-specific uses in `why` or `comment` when they matter for
  Methods writing.
- If reconciliation would hide a real difference, keep separate rows.
- If one input traces deeper than the other, keep the deeper trace.
