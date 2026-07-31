---
name: archived-manuscript-feedback-preflight
description: "Archived historical copy of manuscript-feedback-preflight, retired from the active redraft loop on 2026-07-31. Use only when explicitly auditing or comparing the former pre-drafting feedback audit."
---

# Manuscript Feedback Preflight

## Purpose

Run a final raw-feedback memory pass before drafting or revising manuscript prose. The goal is not to summarize feedback or plan figure work; it is to preserve manuscript-facing signals that may have fallen through workflow boundaries.

Examples of in-scope signals:

- liked prose, titles, phrases, or framing to preserve;
- preferred terminology or wording corrections;
- biological interpretations, asides, uncertainties, caveats, or future-work ideas;
- comments explicitly described as Results-text, prose, narrative, or outside the active figure/planning workflow;
- warnings that a figure/artifact could imply the wrong scientific claim.

## Required Inputs

- A target integration/draft root for the manuscript text pass.
- A `feedback-sources` list containing direct/raw feedback source paths.

If no `feedback-sources` list is present, create one before interpretation. Prefer:

```text
<integration_root>/feedback_archive/feedback-sources.txt
```

If the integration root cannot be inferred, ask the user for it before proceeding.

## Feedback-Sources Creation

Create the list from direct/raw user-feedback files only. Search the relevant integration/redraft roots for filenames such as:

```text
*user-feedback*.txt
*user_feedback*.txt
*raw*transcript*
plan_user_feedback.txt
plan_feedback_transcript.md
```

Exclude parsed or agent-authored feedback artifacts unless no raw source exists:

```text
feedback.md
feedback_intake.md
feedback_sources.md
feedback_understanding/
subagent_outputs/
user_review_mapping.md
user_review_lint.md
review_report.html
critic_review*.md
```

Write paths relative to the repository root. Verify every listed path exists. If duplicate snapshots are retained, keep them as source paths; do not interpret duplicates during source-list creation.

## Mandatory Subagents

This skill requires subagent tooling. Before reading feedback sources, verify that a callable subagent/delegation tool is available. If it is not visible, use tool discovery if available. If no callable subagent tool is available, abort the workflow and tell the user that manuscript feedback preflight cannot be done to standard without subagents.

When multiple feedback files are listed, launch exactly one subagent per feedback file. Do not batch multiple files into one subagent. Do not let the main agent perform the per-file interpretation itself and call that equivalent.

Each subagent must inspect its assigned raw feedback file plus the current manuscript artifacts needed to decide whether a feedback nugget still matters. It must not edit files.

## Per-File Subagent Task

Give each subagent a prompt with this shape:

```text
You are reviewing one raw user-feedback file immediately before manuscript text drafting/revision for the GlucoseStarvation manuscript.

Source file: <path>
Current integration/draft root: <path>

Extract only manuscript-facing signals likely to be missed by the standard redrafting workflow.
Do not summarize implementation feedback. For each candidate nugget, trace the artifact the user was discussing and check whether that artifact still exists or has been superseded in the current manuscript integration/draft root.

Return Markdown with one record per surviving candidate:
- raw_source
- user_comment_anchor: short quote/paraphrase plus line number if available
- artifact_discussed
- current_artifact_status: present | changed | removed | superseded | unclear
- still_relevant: yes | no | unclear
- manuscript_residue: what should survive into prose revision
- confidence: high | medium | low
- follow_up_needed

Also list rejected candidates briefly with the reason they were rejected.
```

## Main-Agent Synthesis

After all subagents return, synthesize only surviving manuscript-facing residues. Do not produce a general feedback summary.

For each retained item, include:

```text
source:
artifact_discussed:
current_artifact_status:
manuscript_residue:
confidence:
follow_up_needed:
```

Drop candidates when:

- the artifact no longer exists and no manuscript-facing implication remains;
- the feedback was local figure QA already resolved by the final artifact;
- later feedback clearly superseded it;
- it only says a figure/package was approved or rejected;
- it is an implementation instruction with no prose, claim, narrative, caveat, or interpretation residue.

Keep candidates when:

- the final artifact still carries an interpretive risk raised by the user;
- the user liked wording or framing that should be preserved;
- the user supplied terminology that should guide prose;
- the user proposed or hinted at a biological interpretation, caveat, uncertainty, or future direction;
- the user identified narrative order or cross-section placement intent;
- the comment was explicitly for Results text or outside the figure-generation workflow.

## Output

Write a concise Markdown audit at:

```text
<integration_root>/feedback_archive/manuscript-feedback-preflight.md
```

The audit must include:

- the `feedback-sources` list used or created;
- one line per subagent assignment;
- retained manuscript-facing residues;
- rejected or obsolete feedback categories;
- unresolved artifact-status checks that need user or manuscript-drafting attention.

Do not update manuscript prose in this skill unless the user separately requests drafting/revision after the preflight audit is complete.
