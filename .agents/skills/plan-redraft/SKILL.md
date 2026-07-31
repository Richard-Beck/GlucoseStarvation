---
name: plan-redraft
description: "Plan the next manuscript redraft step from a target containing
some or all elements of a valid manuscript, user input in the form of a prompt
or feedback, and optional in-progress redrafting work. Use when Codex needs to
assess the current redraft state, choose the next stage in the redraft sequence,
serve handoff packages to consumers, route terminal manuscript assembly, and
update the active plan."
---

# Plan Redraft

## Purpose

Use this skill to run one planner turn for a manuscript redraft.

The planner receives a **target**, user input in the form of a prompt or
feedback, and optionally in-progress redrafting work. It assesses the current
state, determines the next redraft stage, routes required inputs to the
appropriate stage owner, and updates the active plan.

This skill coordinates redrafting. It does not perform stage-owner work itself:
do not interpret, summarize, classify, or close feedback for a phase owner, run
substantial analyses, regenerate figures, rewrite manuscript sections, edit
claim/evidence files, or replace final manuscript outputs unless the user
explicitly expands scope.

Several planner gates below require fresh-context subagents. If subagents are
unavailable, or a required fresh-context launch fails, hard-escalate to the user
and stop that planner turn rather than performing the gated assessment locally.

## Core Terms

**Target** means the manuscript-facing object being redrafted. It may be a
rendered manuscript, manuscript directory, integration root, figure set, prose
draft, planning root, or other bundle that contains some or all elements of a
valid manuscript.

**Valid manuscript** means a manuscript state with enough internally traceable
pieces to evaluate and improve it as a scientific manuscript. Elements may
include analysis outputs, figures, legends, Results text, Methods text,
Abstract/Introduction/Discussion text, claims, evidence, provenance, references,
supplements, rendered HTML/PDF, and validation or status notes.

**User input** means the prompt or feedback that defines what the redraft should
respond to. Raw user wording is authoritative when available.

**In-progress redrafting work** means an existing redraft root, active plan,
handoff package, stage-owner output, review note, or partial artifact created during
an earlier planning or execution turn.

## Core Contract

Manuscripts are tightly connected documents. Changes made early in the redraft
often alter what later sections should say. The planner should therefore advance
the redraft through a stable order:

```text
analysis -> figures -> claim graph update -> Results -> Methods -> Legend Writing -> Abstract/Introduction/Discussion -> Manuscript Assembly
```

The planner's job on each turn is to determine where the redraft currently sits
in this order, choose the stage owner for the next step, and route the inputs
that owner needs. The stage owner is responsible for interpreting the user input
within its phase scope and determining the executable work. Manuscript Assembly
is the terminal packaging stage: route to it after the content and evidence
stages are complete, skipped by explicit user instruction, or blocked with
accepted exceptions.

The planner must route inputs without prejudice. This means the planner may
identify relevant source artifacts, prior stage outputs, unresolved blockers,
approval constraints, and phase boundaries, but must not convert one consumer's
interpretation, recommendations, or "downstream handling" notes into instructions
for the next consumer. Prior consumer outputs are inputs to inspect, not commands
to obey, unless the user explicitly approved them as binding constraints.

A stage may conclude that no change is needed, but that conclusion should be
recorded in the active plan. Treat completed upstream work as input to all later
stages, without adopting its recommendations as planner instructions. The
handoff should preserve upstream outputs by path, status, and scope, and let the
next consumer decide what those outputs mean within its own phase.

Feedback serving is planner-owned; feedback interpretation is consumer-owned.
When feedback is part of the user input, preserve it as raw source paths, exact
current-turn wording, or opaque feedback IDs. Serve a raw feedback bundle to
every stage owner or consumer that receives a handoff. Do not open or read raw
feedback files as part of planning, summarize them in the plan, assign them to
topics, or decide which items are relevant to a phase. Consumers read and
interpret the feedback within their own scopes.

Manuscript Assembly is protected by a feedback-completeness gate. The planner
must not route a first Manuscript Assembly handoff until every pre-assembly
consumer has a current `signed_off` response on every feedback span and every
active feedback item in every tracking root registered for the redraft. After
Manuscript Assembly has returned, the same gate includes the
`manuscript-assembly` consumer itself. Missing, `needs_followup`, and `reopened`
responses are open feedback and do not pass the gate. A skipped stage or accepted
exception does not waive consumer signoff; that consumer must still inspect the
served feedback and record a signed-off no-change or not-applicable response.

## Consumer-generic Instructions

You are the scope owner for {workflow-stage}. Your responsibilities are to parse
user feedback and instructions, determine all user requested actions both
explicit and implicit that are actionable under your scope, and execute them all
(as far as you are able). Guidelines:
1) Read all feedback and instructions, then enumerate all elements which may require actions within
the scope of {workflow-stage}. Write a short list for yourself containing all
enumerated feedback elements and why you think they might require attention within the scope
of {workflow-stage}. If in doubt about a particular element, err towards inclusion.
Write this list to file, treat it as immutable.
2) Before substantive execution, read `<redraft_root>/run_environment.json` when
present and write a short `parallelization_plan.md` in your return root. Identify
independent work packages, planned subagent waves, CPU slots per active tool
worker, memory-heavy tasks, disjoint output ownership, and unavoidable serial
barriers. Record actual parallel execution in the final return.
3) Delegate independent work within your current stage; do not delegate
responsibility for another workflow stage or advance the planner-owned workflow.
4) You may take multiple turns to complete {workflow-stage}. If your locally scoped protocols require
user inputs, or your tasks are too large to complete in a single turn,
return at an appropriate checkpoint and indicate outstanding work.
5) Upstream outputs can be used to refine and inform response to user feedback, but
should not limit or constrain your work.
6) If your locally scoped protocols include a workflow applicable to an open
piece of user feedback, THAT TASK IS IN SCOPE NOW. It must be attempted.
7) Before conclusion, revisit the immutable list/instructions/feedback.
Does your work address each element to the fullest extent possible? If not,
revise your work or request a new turn.
8) Whenever you directly modify ANY directly manuscript facing components (prose, figures, legends),
deliver complete, canonical, manuscript ready output bundles according to
{workflow-stage} guidance. I.e. deliver a full figure bundle even if only one figure was updated,
a full legend bundle even if only one legend was updated, etc.
9) Before returning, write `compatibility_sidecar.json` according to
`.agents/skills/plan-redraft/references/compatibility_sidecar_contract.md`.
Within your scope, identify and authenticate manuscript-facing facts that
another stage may repeat, define, interpret, or depend upon. The planner owns the
format and later comparison; you own the facts reported from your stage.
10) Write `user_blockers.json` under the planner-owned contract in
`references/user_blocker_policy.md`; you own meaning, not disposition.

## Inputs

Require:

- A target containing some or all valid-manuscript elements.
- User input in the form of a prompt or feedback.

Optional but important when present:

- In-progress redrafting work, especially an existing active plan, prior handoff
  packages, stage-owner outputs, review notes, or partial artifacts.

If the target, user input, or in-progress work is ambiguous, inspect local
context first. Ask the user only when the next step cannot be chosen safely.

## Redraft Root

When no redraft root exists, create one under:

```text
agent-dev/manuscript_redrafts/<run_id>/
```

Infer a short run id from the target and current date. Do not overwrite existing
roots.

When in-progress redrafting work already exists, continue in that root unless
the user explicitly asks for a new one. Keep new planning outputs inside the
redraft root and reference large target artifacts by path instead of copying
them.

When creating a new redraft root, or when continuing an existing redraft root
that does not contain `<redraft_root>/global_contract.md`, launch a fresh-context
subagent before choosing the next stage. Provide this prompt verbatim, replacing
`<redraft_root>` with the actual redraft root path:

```text
Read the plan redraft skill and determine all the subskills contributing (as named consumers/stage owners) to the planning workflow. Match the intake requirements for the manuscript assembly skill to the output contracts for each subskill stage. Thus determine the expected shape of the output bundle from each substage needed to satisfy the assembler. Route output to <redraft_root>/global_contract.md
```

The generated contract must also require the planner-owned return overlays
`parallelization_plan.md`, `compatibility_sidecar.json`, and
`user_blockers.json`, even when the consumer skill does not name them.

The planner may inspect the returned `global_contract.md`, but must not
synthesize or repair it in the planner context. If the subagent cannot be
launched with a fresh context, hard-escalate to the user.

## Planner Workflow

Execute this workflow in a single planner turn:

1. **Intake**
   Identify the target, user input, and any in-progress redrafting work. Record
   target and redraft-work paths, versions, and limitations. Record feedback
   sources as raw paths, exact prompt text already visible in the current turn,
   or existing opaque feedback IDs.
2. **Assess valid-manuscript elements**
   Inventory which manuscript elements are present, missing, stale, partial, or
   already under revision.
3. **Assess the just-returned bundle**
   Read the active plan, most recent planner handoff, and the stage-owner bundle
   returned from that handoff. Determine which stages are complete, blocked,
   awaiting review, or ready for a next handoff. Launch a fresh-context subagent
   to compare only that just-returned bundle against the relevant stage section
   of `<redraft_root>/global_contract.md`, its handoff, and its declared inputs.
   Do not freshly audit other accepted stage bundles. Carry their passing
   assessments forward by exact authenticated bundle hash. If a prior bundle
   changed out of band, treat it as newly returned and audit it. Write the
   return-local report under:

   ```text
   <redraft_root>/contract_assessments/<turn_id>_<stage>_return_assessment.md
   ```

   If no new bundle has returned, reuse existing authenticated assessments and
   record that no return-local audit was due. If the just-returned bundle fails,
   route it back to its owner to repair the nonconforming return. Do not advance
   until the repair is complete, explicitly skipped, or accepted as an
   exception. If a required fresh-context subagent cannot be launched,
   hard-escalate to the user and stop the turn.
   After accepting a returned bundle, read
   `references/user_blocker_policy.md` and reconcile its `user_blockers.json`.
   On `USER_ACTION_REQUIRED`, update the plan and stop without another handoff;
   on `CONTINUE`, carry the ledger and disposition forward.
4. **Choose the next stage**
   Advance through the standard order. Do not skip a stage merely because the
   user's prompt mentioned a later section; upstream state may still determine
   what later text should say. A stage may be marked skipped only when the user
   explicitly instructs that it should be skipped or when the active plan records
   an accepted exception. After all prior stages are complete, skipped, or
   accepted-exception blocked, run the pre-Assembly compatibility gate and then
   the feedback-completeness gate before choosing Manuscript Assembly. If the
   compatibility gate fails, route to the earliest mechanically involved
   consumer as defined below. If the feedback gate fails, route to the earliest
   consumer lacking signoff. Mark the routed stage and every later stage through
   Manuscript Assembly pending replay. Continue through every later consumer in
   workflow order when required by the feedback gate, even when a later stage
   previously returned a conforming bundle. Repeat the applicable gates until
   Manuscript Assembly is reached with compatible sidecars and no open feedback.
5. **Choose the stage owner**
   Identify the owner responsible for the chosen phase. Use the phase scopes
   below to determine ownership, not a hard-coded dependency on another skill.
6. **Guard against instruction propagation**
   Before writing the handoff, review all candidate handoff language for
   imperatives such as "must", "should", "replace", "refresh", "rewrite",
   "adopt", "reject", "promote", or "use". Keep imperatives only when they come
   from user instructions, safety/reproducibility requirements, approval gates,
   feedback-manager rules, or the phase boundary itself. Reword prior-consumer
   recommendations as source artifacts to inspect, not as planner commands.
7. **Route required inputs**
   Write one handoff package that gives the stage owner the target context, user
   input source bundle, in-progress work, upstream outputs, and binding
   constraints needed to decide the phase-specific work. Do not tell the owner
   what outcome to reach or which upstream recommendation to implement.
8. **Update the active plan**
   Record the state assessment, chosen stage, stage owner, handoff path,
   blockers, and expected return artifacts.
9. **Report**
   Tell the user the active plan path, next stage, handoff path, and any blocker
   or user decision needed.

If no stage-owner handoff is appropriate because user approval is needed, write
a review/checkpoint package instead and update the active plan with the pending
decision.

## Redraft Phase Scopes

These scopes define the normal phase owners. They are routing boundaries for the
planner, not instructions for the planner to design the detailed work.

- **Analysis**
  Scope: analysis validity, bug checks, sensitivity checks, model reruns,
  statistic definitions, provenance audits, and explicit no-change decisions
  about analytical support.
- **Figures**
  Scope: figure revisions, panel selection, visual encoding, ordering,
  supplements, panel descriptions needed for figure review, and figure-local
  provenance. Manuscript-facing legend writing is a later phase.
- **Claim Graph Update**
  Scope: reconciling current figure evidence and any provisional legend inputs
  with claim/evidence graph state, updating claim/evidence mappings, graph
  validation, and explicit no-change decisions about claim/evidence structure
  after figure changes.
- **Results**
  Scope: Results logic, claim strength, evidence framing, caveats, figure
  callouts, and section-level interpretation of evidence-facing artifacts.
- **Methods**
  Scope: measurement definitions, preprocessing descriptions, modeling details,
  analysis provenance, reproducibility details, and methods-facing implications
  of analysis/figure/Results changes.
- **Legend Writing**
  Scope: manuscript-facing legend text, panel descriptions, caption-level
  provenance, and consistency with finalized figures, Results framing, and
  Methods outputs. Legend writing should consume Methods outputs as upstream
  source material so measurement definitions, preprocessing details, and model
  descriptions are aligned.
- **Abstract, Introduction, and Discussion**
  Owner: `serve-manuscript-abstract-introduction-discussion`.
  Scope: prepare the current context bundle, carry forward valid prior bundles,
  and serve the current externally supplied sections or their defined empty
  state.
- **Manuscript Assembly**
  Owner: `manuscript-assembly`.
  Scope: assemble a coherent manuscript release or review package from finalized,
  intentionally skipped, or accepted-exception components. This includes the
  rendered draft, editable sources, final assets, evidence basis, traceability,
  review state, validation reports, rebuild instructions, and package status.
  Manuscript Assembly should validate completeness, consistency, renderability,
  and traceability without performing new analysis, figure generation,
  substantive prose drafting, legend writing, or evidence remapping.

The planner should advance through these scopes in order and route the current
handoff to the first scope that is ready and not complete.

## Assembly Gates

Immediately before a first Manuscript Assembly handoff, and whenever later
replay reaches that boundary, read and execute
`references/assembly_gates.md`. Run its pre-Assembly compatibility gate before
its feedback-completeness gate, including its blocker-ledger check. After Assembly returns, run its full eight-
consumer feedback gate before terminal completion. Do not route Assembly or
write `redraft_complete.json` while an applicable gate is blocked.

## Handoff Package

Write handoffs under:

```text
<redraft_root>/handoffs/<turn_id>_<stage_slug>_handoff.md
```

Each handoff should include:

- consumer-generic instructions, verbatim.
- stage and stage owner;
- target and relevant valid-manuscript elements;
- raw feedback source bundle served to the consumer;
- latest feedback-completeness audit and replay-pass identifier when the handoff
  belongs to a feedback-gate replay;
- global contract path and, when present, the just-returned bundle's latest
  return-local contract-assessment report;
- `<redraft_root>/run_environment.json` when present;
- `<redraft_root>/run_policy.json`, the user-blocker contract, ledger, and latest
  reconciliation report;
- the compatibility sidecar contract and any current compatibility-gate report;
- upstream outputs and stated conclusions to inspect, clearly attributed to their
  source consumer;

For a Manuscript Assembly handoff, also include:

- the desired or inferred assembly root, if known;
- current status of each content/evidence component as complete, skipped,
  blocked, or accepted exception;
- manuscript sources, any existing draft, and rebuild inputs;
- final asset, evidence, traceability, review-state, and validation inputs by
  path;
- the blocker ledger, requiring deferred decisions to be written as
  `review_state/deferred_user_decisions.md` and package warnings;
- expected return artifacts: assembly root, package status, rendered draft path
  if produced, validation report, blockers, warnings, and accepted exceptions.

The handoff should not define the stage owner's detailed work plan. It should
route the required inputs and state the phase scope. The stage owner interprets
the feedback, applies the relevant workflow, and determines the exact work.
Do not summarize feedback or include planner-written feedback interpretation in
the handoff. If a feedback archive or feedback-manager index exists, pass its
source paths or IDs as handles without signing off or reclassifying items.

Do not launder upstream-consumer instructions through planner voice. If a prior
stage output contains language like "required downstream handling", include the
artifact path and, if needed, a neutral note such as "the analysis output
contains downstream recommendations for this consumer to evaluate." Do not
rewrite those recommendations as "the next owner must..." or "the next owner
should..." unless the user has explicitly made them binding.

Allowed:

```text
Analysis returned `stage_outputs/analysis/routing_note.md`, which contains
downstream recommendations about FG2. The handoff lists it alongside the target
figure package and feedback packet as source material for the figures
owner.
```

Not allowed:

```text
The figures owner must refresh Figure 2, Figure S3, and Figure S4 from
corrected exports.
```

## Active Plan

Maintain:

```text
<redraft_root>/active_plan.md
```

The active plan is the authoritative planner state. Keep it concise and current.
It should include:

- target summary and valid-manuscript element inventory;
- feedback source register and service ledger, with no semantic summary;
- latest feedback-completeness audit, gate status, replay-pass number, earliest
  consumer missing signoff, and the workflow suffix currently pending replay;
- current run-environment/run-policy paths, blocker ledger/dispositions, and
  latest pre-Assembly compatibility-gate status;
- current stage status table;
- global contract path and latest contract-assessment report path;
- completed handoffs and returned artifacts;
- active handoff, if any;
- blockers and user decisions needed;
- next expected planner action;
- notes on whether any redraft state files need maintenance.

Do not let the active plan become a transcript dump. Link detailed stage outputs,
handoffs, and review notes instead of copying them wholesale.
For feedback, record which raw sources or IDs were served to which consumer and
where the consumer should return its interpretation or signoff, if applicable.

## Optional Review Surface

When the next step is user review rather than stage-owner execution, write a concise
HTML or Markdown checkpoint in the redraft root. It should tell the user:

- what target artifact to inspect;
- what decision or approval is needed;
- what will happen next;
- which parts of the redraft sequence are complete, blocked, or pending.

Free-form user feedback is acceptable. Do not require table editing.

## Validation

Before finishing a planner turn, verify:

- target paths and user-input sources resolve or limitations are recorded;
- in-progress work was inspected when present;
- `<redraft_root>/global_contract.md` exists, or a fresh-context subagent was
  launched to create it and any failure was escalated to the user;
- active-state assessment audits only the just-returned bundle, carries prior
  passing assessments by authenticated hash, and records when no new return was
  due for assessment;
- any just-returned or out-of-band changed bundle that fails its stage contract
  is routed back to that stage owner unless explicitly skipped or accepted as an
  exception by the user;
- the active plan reflects the current stage status;
- the chosen next stage follows the standard redraft order or records a clear
  reason for waiting;
- the handoff package exists, or the active plan clearly states why no handoff
  was created;
- every handoff requires the planner-owned return overlays and routes the run
  environment and run policy when present;
- the handoff routes inputs without prescribing the stage owner's interpretation,
  detailed work plan, or expected conclusion;
- the handoff was checked for instruction propagation: prior-consumer
  recommendations are attributed as source material rather than restated as
  planner directives;
- feedback sources are served as raw paths, exact current-turn wording, or opaque
  IDs without planner-written summaries or classifications;
- before any Manuscript Assembly handoff or terminal completion, a current
  feedback-completeness audit exists and every required consumer has signed off
  every span and active item across every registered tracking root;
- before a first Manuscript Assembly handoff, and whenever later replay reaches
  that boundary, a current compatibility-sidecar gate reports `PASS` or `WARN`,
  with every accepted exception routed to Assembly;
- any failed feedback-completeness gate was routed to the earliest consumer with
  missing signoff and every stage from that point through Manuscript Assembly is
  marked pending replay;
- skipped or accepted-exception stages have feedback-manager signoff rather than
  being treated as exempt from the feedback gate;
- every accepted return has a reconciled blocker sidecar; interactive awaiting
  decisions stop without another handoff, while deferred decisions reach
  Assembly review state as warnings;
- `<redraft_root>/redraft_complete.json` is absent unless Manuscript Assembly has
  returned and a full eight-consumer feedback-completeness audit passes;
- substantial analysis or manuscript-validity blockers are visible;
- redraft state maintenance decisions are recorded.

## Completion Standard

Finish each planner turn with:

- an updated `<redraft_root>/active_plan.md`;
- an existing or newly subagent-generated `<redraft_root>/global_contract.md`;
- a current return-local contract-assessment report from step 3, or an explicit
  record that no new return required assessment;
- a handoff package or review/checkpoint package, unless blocked;
- a concise final response reporting the redraft root, target, user input source,
  current stage, stage owner, next handoff or checkpoint path, and blockers.

For the terminal Manuscript Assembly planner turn, the active plan should record
the assembly handoff or returned assembly root as the final redraft package
destination, and should mark any skipped or accepted-exception stages explicitly
rather than leaving them pending. It must also record a passing full-consumer
feedback-completeness audit. A returned assembly package is not terminal when
any span or active item remains unsigned for any required consumer; reroute and
repeat the workflow suffix until the gate passes at Manuscript Assembly. Write
`<redraft_root>/redraft_complete.json` only after that terminal gate passes.
