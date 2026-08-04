# Redraft Workflow DAG Survey

Generated: 2026-07-07 15:23:49 EDT

## Scope

This report records a survey of
`.agents/skills/plan-redraft/references/redraft_workflow_dag.json` against the
actual direct artifact dependencies implied by the skill definitions for the
redraft workflow phases.

Abstract, Introduction, and Discussion were intentionally ignored for this
survey.

The DAG declares that each edge represents one direct artifact consumption, not
a stage-order rule or transitive dependency. The survey therefore treated an
edge as correct only when the consumer skill's input/source contract actually
expects the producer's artifact.

## Sources Inspected

- `.agents/skills/plan-redraft/references/redraft_workflow_dag.json`
- `.agents/skills/plan-redraft/SKILL.md`
- `.agents/skills/analysis/SKILL.md`
- `.agents/skills/manuscript-figure-workflow/SKILL.md`
- `.agents/skills/manuscript-figure-workflow/references/polishing.md`
- `.agents/skills/manuscript-figure-workflow/references/figure_set_integration.md`
- `.agents/skills/manuscript-figure-workflow/references/semantic_interpretation.md`
- `.agents/skills/claim-graph-integration/SKILL.md`
- `.agents/skills/manuscript-legend-writing/SKILL.md`
- `.agents/skills/results-text/SKILL.md`
- `.agents/skills/method-table-provenance/SKILL.md`
- `.agents/skills/method-table-provenance/references/node_classification.md`
- `.agents/skills/method-table-provenance/references/methods_spine.md`
- `.agents/skills/method-table-provenance/references/methods_drafting.md`
- `.agents/skills/manuscript-section-change-assessment/SKILL.md`
- `.agents/skills/manuscript-feedback-preflight/SKILL.md`
- `.agents/skills/manuscript-draft-rendering/SKILL.md`

## Survey Findings

The DAG is broadly consistent with the redraft narrative order, but it is not a
complete direct artifact-consumption graph.

The strongest supported direct dependencies are:

- `analysis -> figures`: figure generation and polishing normally consume
  current analysis outputs or no-change analysis decisions.
- `figures -> claim_graph_integration`: claim graph integration requires a
  current evidence package, which the integrated figure set provides.
- `figures -> legends`: integrated legend writing consumes the integrated figure
  package, including manifest, final images, semantic interpretation index, and
  package-local legends.
- `figures -> results`: Results requires a figure-description source, preferably
  semantic interpretations produced or indexed as part of figure integration.
- `figures -> methods`: Methods-facing provenance starts from integrated figure
  manifests, package provenance, semantic interpretations, legends or package
  legends, and generating scripts.
- `legends -> final_assembly_validation`: final rendering consumes integrated
  legend text, though the desired placement and validation ownership for legends
  remains unsettled.
- `results -> final_assembly_validation`: final rendering consumes Results
  sidecars.
- `methods -> final_assembly_validation`: final rendering consumes Methods prose
  or equivalent section sidecars.

## Discrepancies And Feedback

| Survey item | Survey assessment | User feedback / disposition |
|---|---|---|
| `manuscript-feedback-preflight` is absent from the DAG. | It is a hard gate in the current `manuscript-draft-rendering` and `manuscript-section-change-assessment` skill definitions, but absent from the DAG. | This absence is intentional. The workflow is expected to deprecate `manuscript-feedback-preflight`, so do not add it to the DAG on this basis. |
| `manuscript-section-change-assessment` is absent from the DAG. | This is a real omission relative to the current skill contracts. Section-change assessment consumes prior draft text, current figure set, semantic interpretations, claim graph, claim/evidence audit, feedback preflight, and legends where applicable; its outputs may inform Results, legends, and rendering. | This is a key unresolved design issue. It needs separate consideration rather than an immediate mechanical DAG patch. |
| Missing `legends -> results`. | `results-text` can consume figure legends as substitute or augmenting figure-description sources. | Leave as optional only. The current workflow should not prefer or encode this path. |
| Missing `claim_graph_integration -> final_assembly_validation`. | `manuscript-draft-rendering` consumes or validates `claim_graph_integrated.json` and requires user-fixed claims to appear in forceful manuscript prose. | Confirmed as a true discrepancy. |
| Missing `legends -> methods`. | `method-table-provenance` reads `integrated_figure_legends.md` or package-local `legend.md` for manuscript-facing terminology, in addition to figure manifest, package provenance, semantic interpretations, and scripts. | Mark for further examination. |
| `figures -> methods` labeled only as `figure provenance` is incomplete. | The Methods workflow is more than direct figure provenance: method-table provenance, optional node classification, Methods spine, and Methods drafting together form the Methods stage. | Accepted with caveat. Treat those substeps collectively as "Methods" for the workflow phase, but the artifact label may need refinement. |
| `legends -> final_assembly_validation` labeled only as `legend text` is incomplete. | Rendering currently expects `integrated_figure_legends.md` and a `legend_validation_report.md` with pass or accepted exceptions. | Keep the general dependency under review. User is unsure where legends should belong. Preference is to remove `legend_validation_report.md` as a standalone artifact and make legend validation either part of legend writing or part of final report validation, as appropriate. |
| `claim_graph_integration -> results` may overstate a hard dependency. | `results-text` requires figure-description sources; claim graph excerpts, especially locked claims, are optional. | No explicit follow-up decision recorded in this turn. |
| `section_change_assessment -> results` is absent. | `results-text` may receive section-change assessment annotations as optional manuscript-context inputs. | No explicit follow-up decision recorded beyond the broader section-change assessment design issue. |
| `section_change_assessment -> final_assembly_validation` is absent. | `manuscript-draft-rendering` prefers section-change assessment outputs when present. | No explicit follow-up decision recorded beyond the broader section-change assessment design issue. |

## Open Design Questions

1. How should section-change assessment fit into the redraft workflow DAG, if at
   all? The current skill makes it a substantial dependency node, but the user
   has flagged this as unresolved design work rather than an immediate DAG edit.
2. Should `claim_graph_integration -> final_assembly_validation` be added as a
   direct edge for the final rendered manuscript's locked-claim and graph
   validation checks?
3. Should Methods consume legends directly, or should legend-derived terminology
   remain folded into the broader figure/integration package consumed by
   Methods?
4. Where should integrated legends live in the workflow: inside figures/legends,
   as a separate phase, or as part of final assembly? Related preference:
   eliminate `legend_validation_report.md` as an independent durable artifact
   unless a future design explicitly preserves it.
5. Should the `figures -> methods` artifact label be broadened from `figure
   provenance` to something like `figure package and methods provenance inputs`,
   while still treating method-table provenance, node classification, Methods
   spine, and Methods drafting as one Methods phase?

## Maintenance Note

This report is a reference asset only. It does not change
`redraft_workflow_dag.json`, `SKILL.md`, or `docs/project_map.txt`.
