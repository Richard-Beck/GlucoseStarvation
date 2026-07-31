---
name: archived-manuscript-figure-draft-reviewer
description: "Archived historical copy of manuscript-figure-draft-reviewer, retired from the active redraft loop on 2026-07-31. Use only when explicitly auditing or comparing the former independent figure-draft review workflow."
---

# Manuscript Figure Draft Reviewer

## Purpose

Use this skill to perform an independent first-principles audit of a manuscript figure drafting package.

The reviewer protects the user's actual intent from being silently transformed by intermediate workflow artifacts. Treat drafter-authored files such as `feedback_intake.md`, `drafting_panels.md`, directive tables, and `review_report.html` as evidence to inspect, not as authoritative statements of the user's intent.

For visual-QC expectations, read `.agents/references/manuscript_figure_style.md`.

## Review Standard

Answer four questions:

1. Did the draft deliver what the user actually asked for?
2. Did the draft preserve what the user did not ask to change?
3. Can both claims be verified from raw feedback, prior figure/code, current local code, and rendered outputs?
4. Is the package reviewable enough for a human to approve or request revision?

The two central failure modes are:

- failure to deliver something the user asked for
- modification of something the user did not ask to change

Demand revision when either failure is present or when the required comparison cannot be performed.

## Inputs

Read the artifacts needed to reconstruct user intent and compare versions from first principles:

- all relevant raw user feedback bearing on the target figure, including prior-version, sibling-package, cross-package, polishing, integration, or redraft feedback when it affects figure scope
- the last user-reviewed version of each affected figure or subpanel, usually from the prior polishing step or prior manuscript iteration
- the code that generated that last user-reviewed version
- the current local drafting package, including local scripts, copied prior-code baselines, code diffs, generated PNGs, legends, report manifest, `prior_panel_disposition.csv`, `prior_code_fidelity.csv`, and `review_report.html`
- drafter-authored summaries and status tables only after raw feedback has been read or while actively auditing whether those summaries preserve the raw feedback

If the correct prior user-reviewed figure or code cannot be identified, mark the affected comparison `blocked` and demand revision unless the user explicitly asked for a novel figure with no prior-version fidelity requirement.

## Core Workflow

1. Identify the figure package scope and the current recommended outputs.
2. Read raw user feedback relevant to that scope. Do not rely on the drafting package's directive IDs or report table as the source of truth.
3. Write a short independent reconstruction of user intent before evaluating the delivered draft.
4. Identify the last user-reviewed version for each affected figure/subpanel and the code that generated it.
5. Audit `prior_code_fidelity.csv`, copied prior-code baselines, active local scripts, and `code_diffs/` before accepting any inherited-panel claim.
6. Compare prior code to copied local baseline code and copied local baseline code to current local drafting code. Check for changes to data source, filtering, model inputs, statistical summaries, panel role, visual encoding, labels, ordering, layout, annotations, or interpretation.
7. Compare prior rendered outputs to current rendered outputs. Use visual inspection to catch unrequested changes that are not obvious from code.
8. Check whether every requested change is delivered visibly and materially in the draft outputs, not merely mentioned in prose or hidden in an appendix.
9. Audit the review surface: the human-facing report should make the recommended drafts, important alternatives, prior-code fidelity status, blockers, and caveats easy to inspect.
10. Inspect final PNGs directly under the shared manuscript figure style expectations.
11. Produce a reviewer report with a proceed/hold decision and concrete required revisions.

## Delivery Check

For each independent raw-feedback requirement, classify the draft as:

- `delivered`: materially present in the relevant output
- `partially_delivered`: present but weakened, incomplete, hard to inspect, or moved to a less appropriate output
- `missing`: absent from the delivered draft
- `deferred`: explicitly assigned to another workflow with a defensible reason
- `blocked`: not evaluable because required data, prior version, code, or output is absent
- `not_a_drafting_request`: feedback is text-only, integration-only, or otherwise outside figure drafting scope

Treat missing, diluted, or appendix-only handling of a user-requested figure change as a substantive finding.

## Unrequested-Change Check

Compare the current local subpanel generation code against the code for the last user-reviewed version. Flag any current-code deviation that is not requested or clearly implied by raw user feedback.

Common unrequested deviations include:

- changing data sources, inclusion/exclusion rules, model runs, filters, normalization, or aggregation
- replacing raw data with summaries when raw display was requested, or replacing summaries with raw displays without instruction
- changing axis definitions, units, scales, color encodings, symbols, facet variables, panel order, labels, annotations, or legend meanings
- adding, dropping, merging, or repurposing panels outside the requested scope
- changing figure numbering, subpanel identity, or the scientific role of a panel
- carrying forward an exploratory redesign as the recommended draft without explicit user support

If code comparison cannot be performed because prior code is missing, current code is not local, or provenance is ambiguous, mark the relevant item `blocked` and demand revision.

## Prior Code Fidelity Gate Audit

For every panel that descends from a prior user-reviewed artifact, require local evidence that the drafter started from prior reviewed code when that code exists.

Check that:

- `prior_code_fidelity.csv` exists and has one row for every affected prior, replaced, novel, or blocked panel/subpanel
- each inherited row identifies `panel_id`, `inheritance_class`, `prior_png_path`, `prior_code_path`, `copied_baseline_code_path`, `active_local_script_path`, `diff_path`, `allowed_change_directive_ids`, `fidelity_status`, and any `blocker`
- inherited panels classified as `inherited_preserve`, `inherited_targeted_fix`, or `inherited_move` have copied baseline code inside the local drafting package, usually under `prior_code/<panel_id>/`
- active inherited-panel scripts live in the local drafting package, usually under `scripts/`, and are copied/adapted from the local baseline rather than newly written from scratch
- each inherited panel has a plain-text diff under `code_diffs/` comparing copied baseline code to the active local script
- every substantive code diff is mapped to a raw user-feedback directive id, a purely mechanical local-path/output-path relocation, or a documented blocker
- panels classified as `inherited_replace`, `novel_no_prior`, or `blocked_prior_missing` have a defensible rationale grounded in raw user feedback, the approved plan, or an explicit missing-provenance blocker

Return `BLOCKED` when prior reviewed code should exist but cannot be identified, copied, or compared and the panel is not explicitly novel. Return `NEEDS_REVISION` when prior code exists but was not copied locally, when the draft sources old manuscript/polishing/ad-hoc scripts directly, when an inherited panel lacks a diff, when substantive diffs are not tied to raw user feedback or mechanical relocation, or when a copied PNG substitutes for local generation code.

## Report And Manifest Audit

Review `review_report.html`, manifests, legends, and package notes as downstream presentation artifacts.

Check that:

- recommended drafts appear where a reviewer can find them quickly
- blockers and partial delivery are visible near the top of the report
- prior-code fidelity status is visible near the top of the report, not only in a provenance appendix
- every reviewable final PNG has a legend or review guide
- the report does not claim feedback coverage that is unsupported by raw feedback and rendered outputs
- important requested changes are not only described in prose
- exploratory alternatives are clearly separated from targeted-fix recommendations

Do not accept a drafter-authored directive table as proof of coverage. Use it only to identify possible omissions or contradictions.

## Visual QC

Inspect rendered PNGs directly. Check for:

- figure-level titles, captions, headers, or manuscript prose embedded in figure images
- clipped labels, data, legends, colorbars, or scale bars
- unreadable text or symbols at intended print size
- overlapping plot elements, labels, annotations, or panels
- ambiguous panel order or inconsistent panel labels
- inconsistent color, symbol, line, or annotation meanings across related figures
- excessive spacing, density, or in-image prose
- reused raster material that preserves known defects

If a serious review-visible defect remains in a recommended output, demand revision unless it is explicitly documented as a blocker or user-approved exception.

## Decision Status

Use one of these statuses:

- `PASS`: ready for user review or the next workflow stage
- `PASS_WITH_NOTES`: minor issues only; no blocker to user review
- `NEEDS_REVISION`: important delivery, fidelity, report, provenance, or visual-QC issues remain
- `BLOCKED`: the reviewer cannot evaluate required claims because raw feedback, prior user-reviewed outputs, prior code, current local code, or rendered outputs are missing or ambiguous

Prefer `NEEDS_REVISION` or `BLOCKED` when the draft cannot prove both delivery and preservation.

## Reviewer Report

Write a concise report, usually `draft_review.md` or `reviewer_report.md` in the drafting root unless the user requests another location.

Include:

- decision status
- scope reviewed
- raw feedback sources read
- independent summary of user intent
- prior user-reviewed figure/code sources used for comparison
- prior-code fidelity findings, including missing copies, missing diffs, direct dependencies on old scripts, or unjustified code changes
- missing-delivery findings
- unrequested-change findings
- report/manifest/legend findings
- visual-QC findings
- required revisions and any residual risks

Lead with findings ordered by severity. For each finding, name the affected figure/subpanel, cite the relevant file paths, explain the evidence, and state the required action.
