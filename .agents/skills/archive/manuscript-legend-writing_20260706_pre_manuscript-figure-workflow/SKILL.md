---
name: archived-manuscript-legend-writing
description: "Archived historical copy of manuscript-legend-writing, superseded by manuscript-figure-workflow on 2026-07-06. Use only when explicitly auditing or comparing the pre-promotion legend-writing workflow."
---

# Manuscript Legend Writing

## Purpose

Use this skill to produce or validate journal-facing figure legends after figure-set integration and before manuscript rendering. The skill owns legend prose, legend assembly, and legend validation; it should not draft Results text, regenerate figures, remap evidence, edit the claim graph, or rerun analyses.

## Source Order

1. Read `docs/project_map.txt`.
2. Treat the requested integration root as the source of truth.
3. Prefer these local sources when present:
   - `figure_set_manifest.csv`
   - `figure_numbering_crosswalk.csv`
   - `final_images/`
   - `semantic_interpretation_index.csv`
   - `claim_graph_integrated.json`
   - `section_change_assessment/section_change_assessment.csv`
   - `section_change_assessment/section_assessments/`
   - `unresolved_claim_decisions.md`
   - `section_change_assessment/unresolved_assessment_blockers.md`
   - legacy `manuscript_spine_planning/` contract files, only when working from an older integration root that already contains them
   - upstream polishing legends referenced by the manifest's `legend_path` column
   - user feedback, notes, visual QC, and package inventories relevant to the figure
4. If existing `integrated_figure_legends.md` or `legend_revisions/` files are present, revise them in place unless the user asks for a clean replacement.

## Output Contract

Default outputs live in the integration root unless the user specifies another root:

```text
<integration_root>/
  integrated_figure_legends.md
  legend_revisions/              # optional per-figure working files
    Figure_1.md
    Figure_S9_continued.md
  legend_validation_report.md
```

`integrated_figure_legends.md` is the renderer-facing contract:

```markdown
# Figure Legends

## Figure 1. <short title>
<journal-facing legend body>

## Figure S1. <short title>
<journal-facing legend body>

## Figure S9 continued. <short title>
<journal-facing legend body>
```

Contract rules:

- Include exactly one `##` block for every rendered figure image in `figure_set_manifest.csv`.
- Labels must map reversibly to manifest/final-image IDs: `Figure_1` -> `Figure 1`, `Figure_S8A` -> `Figure S8A`, and `Figure_S9_continued` -> `Figure S9 continued`.
- Do not include no-output disposition rows as figure legend blocks.
- Use stable main-then-supplement order from the manifest or crosswalk.
- Give every block a nonempty short title and nonempty body.
- Keep audit/provenance notes out of the integrated legend body.
- If per-figure `legend_revisions/*.md` files contain audit notes, merge only the journal-facing legend block into `integrated_figure_legends.md`.

## Writing Workflow

1. Build the figure list from the manifest and final images; record expected labels before writing.
2. For each figure, gather visible content and intent from the final image, manifest rows, semantic interpretations, upstream polishing legend, claim graph constraints, section-change assessment notes when relevant, unresolved decisions, and figure-specific feedback.
3. If more than one figure legend needs writing or heavy revision, use subagents in waves. Prefer one rendered figure legend block per subagent with a disjoint output such as `legend_revisions/Figure_3.md`.
   Group multiple legend blocks in one subagent only for tightly coupled continued/split figures or when subagent/tooling limits make one-per-legend delegation impractical; record the reason in the audit notes or `legend_validation_report.md`.
4. Require each subagent to write a complete figure legend plus separate audit notes if needed. The central agent merges only journal-facing legend text.
5. Assemble `integrated_figure_legends.md` in the contract format.
6. Write `legend_validation_report.md` with pass/block status, expected figure count, present legend count, missing/extra legend labels, source files inspected, accepted exceptions, and unresolved legend decisions.

## Legend Standard

For each figure, make the legend understandable without the main text while avoiding a mini-Methods or duplicate Results section.

Include:

- Figure number and short title or summary.
- Panel-by-panel descriptions in panel order.
- Definitions of symbols, colors, lines, abbreviations, units, scale bars, and error bars where relevant.
- Statistical essentials where relevant: what intervals/error bars show, sample size or `n`, test used, and significance notation.
- Only enough method detail to interpret the figure.
- Necessary scope labels, such as SUM-159-fuse-excluded versus all-cell-line, where the figure uses them.
- Explicit caveats needed to read the figure correctly, such as model-implied latent states versus direct measurements.

Exclude:

- Source roots, generation commands, source legend paths, checksums, revision notes, and workflow labels.
- Long procedural detail that belongs in Methods.
- Repeated Results interpretation that is not needed to read the figure.
- Claims stronger than the integrated claim graph, current section-change assessment context, or visible figure content supports.

## Validation

Before finishing, check:

- `integrated_figure_legends.md` exists and parses into one block per expected figure.
- No expected figure is missing and no extra figure labels are present.
- Figure labels map to manifest/final-image IDs, including continued and letter-suffixed supplemental figures.
- Every manifest panel for each figure is described or an omission is documented in `legend_validation_report.md`.
- Titles and bodies are nonempty and journal-facing.
- Internal project jargon, source paths, commands, and provenance are absent from legend bodies.
- Unresolved-decision files were checked and relevant caveats were handled or recorded as unresolved.
- Legend output does not contain Results prose sidecar material, claim-graph edits, or new analysis claims.

## Handoff To Rendering

The manuscript renderer should consume `integrated_figure_legends.md` only after `legend_validation_report.md` reports pass or accepted exceptions. If the renderer finds missing, extra, malformed, or stale legend blocks, it should launch this skill or a legend-writing subagent rather than writing legends inline in renderer code.
