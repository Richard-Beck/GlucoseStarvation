---
name: method-table-provenance
description: Build method-table provenance for manuscript figure subpanels in the GlucoseStarvation project. Use when Codex is asked to trace one or more manuscript-visible figure subpanels, plotted quantities, or panel-level data objects backward through figure-set manifests, package provenance, legends, semantic interpretations, scripts, and intermediate data objects into an id/parent/what/why/comment table for Methods writing.
---

# Method Table Provenance

## Purpose

Use this skill to generate Methods-facing provenance rows for manuscript subpanels. The output is a construction table, not polished Methods prose.

The table schema is:

``` text
id | parent | what | why | comment
```

## Roles

-   Overseer: assigned multiple subpanels, or assigned reconciliation over multiple per-subpanel chains.
-   Tracer: assigned one manuscript subpanel to trace into a method-table provenance chain.
-   Merger: assigned exactly two existing provenance chains to reconcile into an additional merged chain.

## Shared Guidelines

Start from the requested integration root when supplied. If it is not supplied, use the current integration root recorded in `docs/project_map.txt`.

Use `id` for the generated object or action node represented by the row. Durable objects should use their stable filename or a concise normalized object id. In-memory transformations may use a concise synthetic id such as `F2b_plot_table`.

Use `parent` for the immediate upstream object or action used to generate `id`. Do not use `parent` as a grouping label or section name.

Allow multiple parents only when `id` is produced by a single merge, model, plot, render, or integration step that directly consumes multiple immediate upstream streams. When using multiple parents, make `what` name the combining action and use `comment` to state each parent's role. Do not use multiple parents to skip an intermediate action or in-memory object.

A row should represent one durable object, named in-memory object, or named action boundary that Methods writing may need to distinguish. Do not combine data construction, plot-table derivation, plot rendering, and figure integration in one row when those steps have different inputs, scripts, assumptions, or Methods relevance. Prefer row boundaries at saved files, named functions, named plot/data objects, model fits, joins/merges across streams, population-defining filters, and terminal provenance boundaries.

Use `what` to state the generation or transformation performed. Keep it operational and specific.

Use `why` to state the role of that generation step in the analysis chain. Keep panel-specific rhetorical framing near the panel/plot rows; upstream data-object rows should have broader generation-facing reasons.

Use `comment` for concrete references and constraints:

-   generated object path(s);
-   source object path(s);
-   script/function/command references;
-   row counts, filters, model scopes, or thresholds when verified;
-   validation status or accepted exceptions if relevant to the method object;
-   terminal reason when the trace stops.

Prefer structured CSV/JSON/RDS metadata and scripts over prose summaries when deciding immediate parent relationships.

For a single subpanel, return a compact Markdown table with the five required columns. For multiple subpanels, return or write distinct subpanel-scoped tables with the five required columns, plus an optional short index. After each table or after the indexed set, add short notes only for important scope issues that affect Methods writing, such as a manifest input list being broader than the subpanel's actual inputs, unresolved reconciliation ambiguity, or a semantic caveat that prevents overclaiming.

Do not write manuscript prose unless the user explicitly asks for prose after reviewing the table.

## Overseer Workflow

For multiple requested subpanels, delegate exactly one subpanel per tracer subagent. If a callable subagent or multi-agent facility is unavailable, stop before producing the multi-subpanel table and ask the user whether to enable/permit subagents or narrow the request to a single subpanel.

Each tracer must return rows only for its assigned subpanel and should not rewrite rows for another subpanel.

For multi-subpanel requests, preserve each tracer's returned row set as a distinct subpanel-scoped section or file. The overseer may add a short index or status summary, but must not merge, compact, coalesce, or rewrite subpanel row sets unless the user explicitly asks for a consolidated table.

When multi-subpanel tracing is complete run the merger workflow as follows:

1.  Treat the per-subpanel tables as fixed inputs; preserve them unmodified and write merged chains as additional outputs.
2.  Use `scripts/score_chain_overlap.py` to select high-overlap chain pairs to merge.
3.  Run pairwise merge tasks in parallel by subagent delegation, one selected pair per merger subagent. Default to 8 parallel disjoint merges per batch unless the user specifies otherwise. If merger subagents are unavailable, stop before producing reconciled multi-chain output and ask the user whether to enable subagents or narrow the request to one pair.
4.  Give each merger exactly two chain tables, plus optional one-line overlap hints.
5.  After each merge round, recompute overlap using merged tables plus any retained tables not modified in the previous round.
6.  Repeat until a single integrated table exists.

## Tracer Workflow

For a single requested subpanel, perform the trace directly.

Read sources in this order:

1.  `figure_set_manifest.csv`: find the exact `figure_id` and `panel_id` row.
2.  Package-level provenance path from the manifest, usually `polishing/provenance.csv`.
3.  Relevant package inventory under `package_inventories/`, if available.
4.  `integrated_figure_legends.md` or the package-local `legend.md` for manuscript-facing terminology.
5.  Linked semantic interpretation markdown for visible object identity, source clarification, caveats, and provenance hints.
6.  Generating scripts and upstream scripts only as needed to resolve exact transformations, inputs, or terminal provenance.

Trace in this order:

1.  Establish the leaf.
    -   Identify the exact manuscript-visible subpanel, its normalized figure id, panel id, source package, short content, and row role from `figure_set_manifest.csv`.
    -   Create a leaf row only after identifying the immediate plotted table/object/action that generated the displayed subpanel.
2.  Find the package-level panel row.
    -   Use the package `provenance.csv` to confirm subpanel image path, generating script or command, direct data inputs, context inputs, and final composite output.
    -   If the integration manifest lists figure-level inputs that are broader than the subpanel actually uses, narrow the trace using the package provenance, plotting script, and semantic interpretation.
3.  Trace upstream one edge at a time.
    -   For each row, ask: "What immediate object or action generated this id?"
    -   Read scripts or structured metadata to verify the edge before writing it when the parent is not obvious.
    -   Do not skip material action nodes. If a function or named in-memory object filters, summarizes, joins, models, maps aesthetics, computes displayed quantities, or renders a file, represent that action/object as a row before tracing to its input table(s).
    -   Continue from saved tables to their generating scripts, then to their input tables/runs, until the trace reaches an explicit terminal.
4.  Terminal rows.
    -   Do not stop at a vague "reporting boundary."
    -   A terminal row must say why it is terminal, such as raw image groups, experimental/acquisition protocol outside current file provenance, an external fixed input, or an unresolved provenance gap.
    -   If provenance is missing or ambiguous, write a terminal gap row rather than inventing a parent.
5.  Write the table.
    -   Order rows from leaf to upstream root unless the user asks otherwise.
    -   Return one table for the assigned subpanel.
    -   Keep paths concise but specific enough to locate the object.
    -   Use backticks for ids, paths, functions, and field names.

## Merger Workflow

When assigned exactly two existing provenance chains, return one additional merged `id | parent | what | why | comment` table.

The merged table must preserve meaningful provenance from both inputs, reconcile duplicate nodes where justified, and keep distinct nodes distinct when they represent different objects or actions.

The purpose of the merger workflow is to remove duplicate rows and refine/improve existing rows. Rough guidelines:

-   Preserve immediate-parent semantics from this skill.

-   Merge rows only when they describe the same generated object or action node.

-   Prefer the most specific verified path/object id and the most complete upstream trace.

-   Preserve panel-specific uses in `why` or `comment` when they matter for Methods writing.

-   If reconciliation would hide a real difference, keep separate rows.

-   If one input traces deeper than the other, keep the deeper trace.

-   Use your own judgement to emit a merged chain that best captures the intent of this skill.

## Optional workflows

Additional optional workflows are documented in references/ folder of this skill.
Optional workflows should not be read or opened by default. Use only when user requests.

name: Node classification
description: Use to classify nodes according to a pre-existing taxonomy.
source: references/node_classification.md

name: Methods Spine 
description: Build a bucket-scoped, node-supported Methods outline from a classified method-table provenance graph. Use after node classification, before Methods prose drafting.
source: references/methods_spine.md

name: Graph Analysis (DEPRECATED)
description: Use to deterministically analyze and reduce a method table graph.
source: references/graph_analysis.md

name: Methods Drafting 
description: Use to draft methods text from a pre-existing methods spine.
source: references/methods_drafting.md
