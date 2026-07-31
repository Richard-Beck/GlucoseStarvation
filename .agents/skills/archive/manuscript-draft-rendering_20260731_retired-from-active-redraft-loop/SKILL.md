---
name: archived-manuscript-draft-rendering
description: "Archived historical copy of manuscript-draft-rendering, retired from the active redraft loop on 2026-07-31. Use only when explicitly auditing or comparing the former standalone manuscript-rendering workflow."
---

# Manuscript Draft Rendering

## Purpose

Use this skill to create or revise a manuscript draft from an already integrated evidence package. This is downstream of manuscript integration and figure polishing. It should render, rewrite, embed, and validate; it should not remap evidence, edit the claim graph, regenerate figures, rerun analyses, or treat the task as integration unless the user explicitly asks for that.

## Source Order

1. Read `docs/project_map.txt`.
2. Treat the requested integration/draft root as the source of truth.
3. Prefer these source files when present:
   - `figure_set_manifest.csv`
   - `final_images/`
   - `claim_graph_integrated.json`
   - `integration_report.md`
   - existing prose draft, prior section sidecars, or legacy `manuscript_bullet_draft.md`
   - `section_change_assessment/section_change_assessment.csv` and `section_change_assessment/section_assessments/` when available
   - `feedback_archive/manuscript-feedback-preflight.md`
   - `manuscript_sections/results/*.md` or another explicit Results-text sidecar
   - `integrated_figure_legends.md`
   - `legend_validation_report.md`
4. Preserve claim graph contents and evidence mappings unless explicitly asked to change them.

## Feedback Preflight Gate

- Before drafting or revising any non-Results manuscript prose or figure legends from this skill, verify that `<integration_root>/feedback_archive/manuscript-feedback-preflight.md` exists for the target integration/draft root.
- If the preflight audit is absent, stale for a different integration/draft root, or missing the retained residues/unresolved drafting checks needed for the requested text pass, run the `manuscript-feedback-preflight` skill first.
- Do not treat ad hoc searches through scattered feedback files as a substitute for the preflight audit.
- Read the preflight audit before launching legend or non-Results prose work. Pass the audit path plus section-relevant retained residues and unresolved drafting checks to those text-writing subagents.
- If preflight cannot be completed to standard because required subagent tooling is unavailable, stop before non-Results manuscript prose or legend drafting/revision and report the blocker. Do not proceed with a main-agent solo feedback-discovery fallback.
- When only mechanically re-rendering an already drafted manuscript with no prose changes, still consume the existing preflight audit if present and record any missing-audit exception in the audit/provenance material.

## Drafting Rules

- Make every `user_fixed: true` claim explicit in the manuscript text.
- State locked claims directly and forcefully first; attach scope limits or caveats immediately after.
- For model-supported selection claims, state the selection claim clearly while saying it is model/simulation-supported rather than direct competition evidence.
- For working hypotheses, state the hypothesis explicitly and label it as a working hypothesis.
- Do not let caution erase the claim. Avoid replacing a locked claim with vague language such as "is associated with" when the graph says the claim is supported.
- Keep caveats visible: direct uptake versus apparent uptake, latent model states versus measured variables, context dependence, exceptions, transfer limits, and sensitivity analyses.

## Results Text Boundary

- Do not draft, revise, or delegate Results prose inside this skill.
- Treat Results sidecars, usually under `manuscript_sections/results/*.md`, as renderer inputs produced by the `results-text` skill.
- If Results text is missing, stale, malformed, or requested for revision, use `results-text` as the owning skill for drafting/delegation before continuing the render.
- The renderer should inspect and inject Results sidecars; it should not silently rewrite Results prose in code.
- Keep Results sidecars available as review artifacts. Inject only the manuscript-facing Results body into the HTML while retaining sidecar notes in audit/review material.

## Text Injection

- Keep Results text outside the renderer source. The renderer should read Results sidecar files and inject their manuscript-facing body into the HTML during generation.
- Keep figure legends outside the renderer source. The renderer should read `integrated_figure_legends.md` and inject its legend blocks during generation.
- Prefer the same sidecar-and-injection pattern for other manuscript sections when prose is substantial or expected to be revised: abstract, introduction, discussion, methods, or a whole-document markdown source can all be external inputs.
- Small rendering labels, anchors, style text, and generated audit boilerplate may remain in the renderer, but scientific prose should be source data rather than hard-coded Python/R strings.
- Validate that every required section sidecar exists before rendering. Fail loudly if an expected Results sidecar is missing or malformed.

## Journal-Facing Language

- Remove or translate internal project jargon from manuscript-facing prose and captions.
- Replace version/workflow labels with content descriptions:
  - `v3`, `v4`, `redraft`, `integration package` -> the combined evidence, the manuscript draft, or the figure set.
  - `WP1`/`WP2`/`WP3`/`WP4`/`WP5` -> measurement, model-free feature, mechanistic model, posterior parameter, or selection-simulation analysis.
  - `FG1`/`FG2`/`FG3`/`FG4` -> omit or translate to the scientific figure topic.
  - `no-SUM` -> `SUM-159-fuse-excluded`, with a definition on first use if the figure label retains the shorthand.
  - `MA1`/`MA2` -> leave-one-cell-line-out posterior robustness or initial-cell-density sensitivity analysis.
- Keep file paths, source roots, script names, commands, and generation bookkeeping out of journal-facing captions.
- Put provenance in an `Audit Trail` section or sidecar note, not in the main figure legend.

## Legend Contract And Delegation

- Do not draft or substantially revise figure legends inside this skill. Use the `manuscript-legend-writing` skill when legends are missing, stale, malformed, or requested for revision.
- Before rendering, check that `integrated_figure_legends.md` exists and follows the renderer contract:
  - top-level `# Figure Legends`;
  - exactly one `## Figure ...` block per rendered figure image in `figure_set_manifest.csv`;
  - labels map reversibly to manifest IDs, such as `Figure_1` -> `Figure 1`, `Figure_S8A` -> `Figure S8A`, and `Figure_S9_continued` -> `Figure S9 continued`;
  - every block has a nonempty short title and body;
  - no no-output disposition rows appear as figure legend blocks.
- Prefer a `legend_validation_report.md` from `manuscript-legend-writing` that reports pass or accepted exceptions. If it is absent, the renderer may run a deterministic parse/check, but should not treat content quality as reviewed.
- If legend writing is required, launch `manuscript-legend-writing` with a preference for one legend-writing subagent per rendered figure legend block, each writing a disjoint file such as `legend_revisions/Figure_3.md`.
- Give each figure-legend subagent the integration root, its figure-specific manifest rows, final image, upstream polishing legend path, semantic interpretation paths, relevant claim-graph constraints, section-change assessment notes when relevant, unresolved decisions, and relevant manuscript-feedback-preflight residues/checks.
- Group multiple legend blocks in one subagent only for tightly coupled continued/split figures or when subagent/tooling limits make one-per-legend delegation impractical; record the reason in the legend audit material.
- Do not ask legend writers to rediscover raw feedback unless the preflight audit identifies a targeted unresolved source check.
- The renderer may consume `legend_revisions/` only after `manuscript-legend-writing` has merged those revisions into `integrated_figure_legends.md`.
- Fail loudly rather than falling back to manifest `short_content` for final manuscript captions, unless the user explicitly asks for a diagnostic preview only.

## HTML Rendering

- Generate self-contained HTML for review unless the user asks for linked assets.
- Embed final figure PNGs directly as `data:image/png;base64,...`; do not rely on local `final_images/*.png` links for portable drafts.
- Keep the renderer reproducible from local source files. Prefer updating a renderer script over hand-editing generated HTML.
- Render main figures near their Results callouts and include supplemental figures in a figure index or supplement section.
- Use audit/provenance text only outside the journal-facing manuscript body.

## Validation

Before finishing, check:

- Renderer compiles and regenerates the HTML.
- Every expected figure appears exactly once unless intentional.
- Every `<img>` source is embedded or otherwise matches the user's requested delivery format.
- HTML anchors are unique.
- `feedback_archive/manuscript-feedback-preflight.md` exists for the target integration/draft root or a no-prose-change rendering exception is recorded.
- Non-Results manuscript prose addresses or explicitly defers the preflight audit's relevant retained residues and unresolved drafting checks.
- All `user_fixed` claims appear in forceful manuscript prose.
- Results sections are injected from sidecar files generated with the `results-text` skill.
- Renderer source does not hard-code Results body paragraphs.
- Results sidecars retain manuscript-facing text and `Revision Notes`.
- Obvious internal jargon is absent from journal-facing prose and captions.
- `integrated_figure_legends.md` parses cleanly, has exactly one legend block for every rendered figure in the manifest, and has no extra figure blocks.
- `legend_validation_report.md` is present with pass status or accepted exceptions, unless the user explicitly requested a diagnostic render before legend review.
- The existing claim graph still validates if the validator is available.

Update `docs/project_map.txt` only when the produced artifact or workflow is durable enough to affect project organization.
