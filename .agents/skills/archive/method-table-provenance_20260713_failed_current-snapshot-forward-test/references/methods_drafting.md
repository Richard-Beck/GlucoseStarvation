# Methods Drafting Workflow

Use this workflow after `methods-spine`.

The goal is to draft a coherent manuscript Methods section while keeping
Methods-relevant judgment with the Methods writer and keeping raw provenance,
file dumps, and execution traces out of manuscript prose.

This workflow is overseer-led. The Overseer manages the multi-turn process, but
the Methods writer owns every draft and revision of manuscript prose.

## Contents

- [Overseer contract](#overseer-contract)
- [Shared Methods guidelines](#shared-methods-guidelines)
- [Role turn contracts](#role-turn-contracts)
- [Operating principle](#operating-principle)

## Overseer Contract

The Overseer manages process quality, information flow, and role boundaries.
All Overseer-relevant instructions are collected in this section.

### Required inputs

* Methods spine.
* Structurally validated current canonical provenance table, if needed for node lookup.
* Classification table derived from that table, if needed for coverage checks.
* Unresolved-issues list from the Methods spine.

### Required outputs

1. Writer-authored Methods draft for each writer turn.
2. Writer-authored source-inspection requests, when the writer needs more facts.
3. Source-inspection fact cards.
4. Reviewer points and reviewer follow-up decisions for each review round.
5. Writer-authored revised draft(s) after reviewer points.
6. Final overseer process and coverage audit.

### Overseer responsibilities

* Deploy subagents with fresh contexts to fulfill writer, reviewer, and inspector
  roles, and escalate to the user if subagents are unavailable.
* Create role-specific assignment packets.
* Preserve writer-authored drafts, revisions, and information requests.
* Batch source-inspection requests.
* Merge duplicate or overlapping requests.
* Split requests that are too broad.
* Route requests to File inspectors.
* Review returned fact cards for quality.
* Pass approved fact cards to the Methods writer.
* Pass reviewer points back to the Methods writer for revision.
* Pass the full review history back to the Methods reviewer on follow-up review
  turns.
* Check final coverage against the Methods spine.
* Audit role compliance and stop the workflow if ownership boundaries are
  violated.

The Overseer does **not** decide whether a detail is Methods-relevant at the
initial request stage. That is the Methods writer's responsibility. The Overseer
may challenge requests only when they are unclear, duplicative, too broad, or
obviously impossible to answer from the requested sources.

The Overseer must not write or revise Methods prose. The Overseer may assemble
files, add audit metadata outside the manuscript text, and request another
writer turn. If the Overseer edits manuscript prose after review, the run is
process-invalid and must be reported as failed.

### Subagent prompt construction

When invoking a writer, reviewer, or file-inspector subagent, the Overseer should
use a minimal role packet rather than restating this workflow in the prompt. The
packet should identify the subagent's role, the current workflow stage, the
allowed input artifacts, the expected output artifact, and the sections of this
file the subagent should read before acting.

Every role packet should include the path to this file:

`/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/GlucoseStarvation/.agents/skills/method-table-provenance/references/methods_drafting.md`

Use this reading map:

- Methods writer: read the opening workflow description, `Overseer Contract`
  paragraphs on writer ownership and role boundaries, `Overseer orchestration
  workflow` steps 1-6, and `Shared Methods Guidelines`.
- File inspector: read the opening workflow description, `Overseer Contract`
  paragraphs on role boundaries, `Overseer orchestration workflow` steps 3-5,
  and the fact-card requirements in the inspector assignment packet.
- Methods reviewer: read the opening workflow description, `Overseer Contract`
  paragraphs on role boundaries, `Overseer orchestration workflow` steps 7-8,
  and `Shared Methods Guidelines`.

Do not replace the reading map with paraphrased workflow guidance. Short
orientation notes are acceptable only when they point to the relevant sections
and clarify the local assignment.

### Overseer orchestration workflow

1. Prepare or update the writer packet.
   * Include the path to this file and the writer entries from the reading map
     above.
   * Include the Methods spine, unresolved-issues list, and any constraints from
     the user.
   * Add approved fact cards, prior writer draft(s), reviewer points, review
     history, and project-owner TODOs when they become available.
   * Exclude raw source-file contents unless the writer is explicitly authorized
     to inspect them.

2. Invoke the Methods writer.
   * The writer returns a draft plus any source-inspection requests and open
     issues.
   * Preserve the writer output as a writer-authored artifact.

3. If the writer requests source inspection, batch requests procedurally.
   * Merge duplicate requests.
   * Split overly broad requests.
   * Group related requests into batches.
   * Reject only malformed or unanswerable requests.
   * Do not reject a request because the Overseer thinks the detail is
     unimportant.

4. Invoke File inspectors.
   * Assign each inspector one request or coherent batch.
   * Include the path to this file and the file-inspector entries from the
     reading map above.
   * Require fact cards in the standard format.

5. Review fact cards procedurally.
   * Approve a fact card if it answers the request, is terse, separates facts
     from uncertainty, avoids implementation dumps, and flags limits honestly.
   * Return or revise a fact card if it narrates code, includes irrelevant
     plotting or file-layout detail, introduces unsupported interpretation, fails
     to answer the request, or hides uncertainty.
   * Do not convert fact-card limitations into prose yourself.
   * Add approved fact cards to the writer packet and return to step 2.

6. Apply the draft-readiness gate.
   * If the latest writer output says source inspection is not needed and any
     remaining factual gaps are classified as project-owner TODOs or open issues,
     go to step 7.
   * If the latest writer output includes source-inspection requests, do not go
     to review. Use steps 3-5 to update the writer packet, then return to
     step 2.
   * If the latest writer output is ambiguous about whether more information is
     needed, return to step 2 with a procedural clarification request.

7. Invoke the Methods reviewer.
   * Include the path to this file and the reviewer entries from the reading map
     above.
   * For the first review, provide only the draft and Methods guidelines.
   * For follow-up reviews, provide the updated draft, the Methods guidelines,
     and the full review history. The review history may include prior reviewer
     points, prior writer responses or open issues, and the prior reviewed
     draft(s). It must not include provenance tables, source-inspection requests,
     fact cards, node classifications, or file contents.
   * Preserve the reviewer points and reviewer follow-up decisions unmodified.

8. Apply the reviewer gate.
   * The reviewer, not the Overseer, decides whether prior reviewer points have
     been satisfied.
   * If the reviewer returns no remaining points to address, go to step 9.
   * If the reviewer returns unresolved prior points or allowed new points raised
     directly by the modified Methods text, add those points to the writer packet
     and return to step 2.
   * The Overseer may ask the reviewer to restate an unclear point, but must not
     decide that a reviewer concern has been satisfied.

9. Terminate or continue.
   * Terminate only when the writer has no further source-inspection requests for
     the current draft, high-priority details are resolved or explicitly marked
     TODO/project-owner input, the reviewer has no remaining points to address,
     all major spine areas are represented, no known unresolved provenance issue
     is overclaimed, and remaining open issues require human/project-owner input
     rather than further file inspection.
   * Do not loop indefinitely to polish prose. Iteration is for resolving
     Methods-relevant factual gaps and reviewer-identified methodological
     problems.

10. Produce the final overseer audit.
    * Audit coverage and process compliance.
    * If the writer did not perform required writer steps, report the run as
      process-invalid rather than filling the gap as Overseer.

### Final audit checks

* All major spine areas are represented.
* All `include_main` obligations are covered, collapsed into broader Methods
  descriptions, or explicitly deferred.
* Primary and sensitivity scopes are consistent.
* Known provenance gaps are not overclaimed.
* Fact-card details were not distorted.
* Provenance language did not leak into prose.
* Reviewer points were cleared by the reviewer or remain explicitly deferred to
  project-owner input.
* Full review history was provided to the reviewer on follow-up review turns.
* Reviewer follow-up points were limited to prior reviewer concerns and issues
  raised directly by modified Methods text.
* Writer-authored drafts and revisions exist for every prose stage.
* The Overseer did not write or revise manuscript prose.

Final audit format:

```markdown
# Final Methods Draft Audit

## Process Compliance

<Whether each role stayed within its boundary. Mark process-invalid if not.>

## Coverage

<Coverage against the Methods spine.>

## Scope Consistency

<Primary/sensitivity/excluded scope handling.>

## Fact-Card Use

<Whether inspected details were used accurately.>

## Reviewer Response

<Whether reviewer points were cleared by the reviewer, remain open for writer
revision, or are deferred to project-owner input.>

## Remaining TODOs

- <TODO or project-owner question.>

## Provenance Leakage

- <Any remaining internal/provenance terms.>

## Recommendation

<Ready for human review / needs targeted writer revision / needs more inspection /
process-invalid.>
```

## Shared Methods Guidelines

These guidelines are visible to both the Methods writer and the Methods reviewer.

### What the Methods section should do

The Methods section should explain, in a reader-facing way, how the study was
performed and how the reported quantities were generated.

A good Methods draft should let a reader understand:

1. **What was studied.**
   Describe the experimental system, analysis units, comparison groups, and the
   main analysis populations.

2. **What was measured.**
   Explain how the raw observations were generated, including assays, imaging,
   annotations, classifications, calibrations, or other measurement steps.

3. **How the analysis data were assembled.**
   Explain how raw observations, metadata, calibration data, and quality-control
   rules were combined into analysis-ready data.

4. **How reported variables and summaries were defined.**
   Define the quantities that appear in the Results: counts, rates, areas,
   ratios, contrasts, endpoints, scores, uncertainty summaries, or other derived
   features.

5. **What models or algorithms were used.**
   Describe the model structure, algorithm, classifier, likelihood,
   parameterization, assumptions, or decision rules at the level needed to
   understand the analysis.

6. **How models or algorithms were fit, selected, or evaluated.**
   Explain fitting, optimization, sampling, training, model comparison, model
   selection, convergence checks, validation, and diagnostic criteria when
   relevant.

7. **How simulations or predictions were performed.**
   Describe the simulated system, initial conditions, inputs, scenarios,
   endpoints, and scoring rules.

8. **How robustness was assessed.**
   Explain sensitivity analyses, alternative analysis populations, leave-one-out
   analyses, held-out tests, threshold checks, or other stress tests.

9. **How the work can be reproduced.**
   Summarize where code, data, configuration files, and important outputs are
   available, without turning the Methods into an execution log.

### How to use the Methods spine

Use the Methods spine as an evidence map, not as prose.

The spine tells you what topics need coverage and which nodes support them. It
does not determine the final paragraph order, subsection titles, or wording.

Do not copy node IDs, local paths, run names, figure-generation labels, export
names, or workflow terminology into the manuscript text. Use those only to verify
coverage and request additional information.

### What to write, and what to leave out

Write about procedures, definitions, assumptions, and analysis choices.

Include details when they define:

* what data or samples were included;
* how measurements were generated;
* how observations were filtered or quality controlled;
* how variables, endpoints, or scores were computed;
* how models, likelihoods, priors, or algorithms were specified;
* how fitting, sampling, optimization, or model selection was done;
* how uncertainty was summarized;
* how simulations were initialized and scored;
* how sensitivity or robustness analyses differed from the main analysis.

Leave out details that are only about:

* file movement;
* table reshaping;
* cache creation;
* figure rendering;
* panel layout;
* plotting scripts;
* local directory structure;
* workflow status;
* package assembly;
* internal run labels;
* provenance bookkeeping.

Mention code or file names only when they are part of a code/data availability
statement or are necessary for reproducibility.

### Handling source files and fact cards

When a file, script, config, or artifact is used, extract the method from it. Do
not describe the file itself as the method.

For example, write about the model equations, fitting rule, threshold,
calibration, or simulation protocol, not about the fact that a particular script
or export existed.

If the spine or fact cards show that a detail is unresolved, do not invent it.
State only what is supported, request more inspection, or leave a clear TODO for
project-owner input.

### Structure and style

Write one coherent Methods section. Do not stitch together independent bucket
summaries.

Use subsection structure only when it helps the reader.

Avoid defensive repetition. If the same analysis population, exclusion, model
set, or sensitivity comparison is relevant in several places, define it clearly
once and refer back to it naturally.

Keep Results and interpretation out of Methods. The Methods should explain what
was done, not argue what the findings mean.

Prefer compact, direct sentences. The first draft should make the structure
clear before trying to be exhaustive.

### Final self-check

Before returning a draft or revision, the writer should check that:

* the reader can tell what was studied;
* the reader can tell how each major reported quantity was generated;
* measurement, data assembly, derived summaries, modeling, fitting, simulation,
  and robustness are not confused with one another;
* important thresholds, filters, assumptions, and endpoint definitions are
  included;
* sensitivity analyses are clearly distinguished from the main analysis;
* unresolved details are marked honestly;
* provenance language has not leaked into manuscript prose.

## Role Turn Contracts

### Methods writer turn contract

This contract applies every time a Methods writer is invoked.

The writer may receive any combination of:

* Methods spine;
* Methods guidelines;
* unresolved-issues list;
* approved fact cards;
* prior writer draft;
* reviewer points;
* review history;
* prior writer notes or project-owner TODOs.

The writer should not freely browse source files unless explicitly authorized.
The writer should use only the supplied packet and should separate manuscript
prose from process notes.

Every writer turn returns:

1. A writer-authored Methods draft.
   * Return the complete Methods draft unless the Overseer explicitly scoped the
     turn to a complete subsection.
   * If a detail is unresolved, write only what is supported and mark the gap as
     a TODO or information request.
2. Information requests, if any.
   * Requests should cover facts the writer needs before the draft can be made
     safer or more complete.
   * If no further source inspection is needed, say so explicitly.
3. Open issues, if any.
   * Identify project-owner TODOs, reviewer points deferred to project-owner
     input, and any reviewer point rejected with brief rationale.

When the packet includes fact cards, reviewer points, or review history, the
writer must handle each material limitation or point by doing one of the
following:

* revising the draft;
* issuing a source-inspection request;
* marking a project-owner TODO;
* explaining briefly why no Methods change is needed.

The writer should not terminate a turn merely because the prose can be improved
later. The writer should terminate only when the current packet has been handled
and any remaining factual gaps have been classified.

Source-inspection request format:

```markdown
## Source-inspection request: `<request_id>`

### Methods area

<Section or topic needing detail.>

### Question

<Specific methodological question.>

### Relevant nodes or files

- `node_or_file`
- `node_or_file`

### Needed output

<Variable definition, model specification, threshold, inclusion rule,
calibration rule, simulation rule, uncertainty rule, etc.>

### Why this matters

<Why the draft cannot safely proceed without this detail.>
```

### File inspector turn contract

Each File inspector receives one request or a coherent request batch.

The inspector returns fact cards only; the inspector does not write manuscript
prose.

Fact card format:

```markdown
## Fact card: `<request_id>`

### Request answered

<One sentence restating the request.>

### Inspected sources

- `node_or_file`
- `node_or_file`

### Methods-relevant facts

- <Terse fact.>
- <Terse fact.>
- <Terse fact.>

### Boundaries / limitations

- <What the inspected sources do not establish.>
- <Any missing provenance or uncertainty.>

### Do not carry into prose

- <Internal labels, file mechanics, plotting details, or workflow terms to avoid.>

### Follow-up needed?

<No, or a specific narrow follow-up request or project-owner question.>
```

### Methods reviewer turn contract

For the first review turn, the Methods reviewer receives only:

* Methods draft;
* Methods guidelines.

For follow-up review turns, the Methods reviewer receives:

* updated Methods draft;
* Methods guidelines;
* full review history.

The review history may include prior reviewer points, prior writer responses or
open issues, and the prior reviewed draft(s). It should contain enough context
for the reviewer to decide whether prior concerns were addressed.

The reviewer does not receive:

* Methods spine;
* provenance table;
* node classification table;
* source-inspection requests;
* fact cards;
* file contents.

On the first review turn, the reviewer returns a concise list of points for the
writer to address. The reviewer does not revise the draft.

On follow-up review turns, the reviewer acts as the gatekeeper for their own
concerns. The reviewer should decide whether prior points have been resolved and
return only remaining points for the writer to address. The reviewer should not
perform a new unrestricted review each round. New points are allowed only when
they are raised directly by modified Methods text or are needed to clarify
whether an original point was resolved.

Reviewer report format:

```markdown
# Methods Reviewer Points

- <Point for the writer to address.>
- <Point for the writer to address.>
```

If the reviewer has no remaining points, return:

```markdown
# Methods Reviewer Points

No remaining reviewer points.
```

## Operating principle

The Methods writer decides what information is needed and owns all manuscript
prose.

The File inspector extracts only requested methodological facts.

The Reviewer judges the draft in a fresh context.

The Overseer controls batching, quality, information hygiene, and role
compliance, but does not write or revise Methods prose.
