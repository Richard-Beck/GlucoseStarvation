# Staged Revision Contract

Use this contract when changing an authenticated manuscript assembly without destroying its sealed state.

## Ownership

Plan Redraft implements this contract through its deterministic controller.
Artifact owners write only to assigned batch locations. Manuscript assembly
receives a reconciled candidate; it does not own the revision workspace.

## Command lifecycle

The controller is a multi-command workspace manager:

1. `authenticate` verifies the existing checksum inventory without writing or
   rendering.
2. `init` repeats authentication, records the baseline identity, and creates
   the workspace.
3. `begin-batch`, `stage-path`, `register-generated`, `register-check`, and
   `delete-path` manage batch state without modifying baseline content.
4. `record` inventories actual overlay changes and exact-hash backreferences. It
   blocks unapproved changes to configured contract-sensitive machinery but
   does not adjudicate manuscript consequences.
5. The agent classifies the observed batch artifacts, traverses the frozen
   manuscript dependency graph, and `register-consequence-audit` validates and
   stores both records.
6. `materialize` composes a review candidate.
7. `preview` invokes the configured renderer against that candidate and writes
   a versioned, visibly unsealed HTML.
8. `set-contract-disposition`, `resolve`, `set-batch-state`, and `status`
   maintain and report revision state.

Use the controller copy stored in the active workspace after initialization.
Authentication and initialization never invoke the preview renderer.

## Required properties

- Treat the sealed assembly as an immutable baseline.
- Exclude the configured revision workspace from the baseline seal.
- Snapshot every baseline file, including checksum-excluded files, and reject
  content or metadata drift while the workspace is active.
- Write edits, generated artifacts, previews, and experimental reports only in the revision workspace.
- Permit multiple ordered revision batches and detect overlapping writes.
- Materialize candidates from the baseline plus selected batches. Normally use
  one candidate per coherent request; record a reason for another candidate.
- Render previews at batch completion or an explicit checkpoint, not after every atomic modification.
- Keep baseline status and validation reports unchanged when a candidate is incomplete or fails.
- Seal only a fully reconciled and terminally validated candidate.

## Conceptual layout

Resolve actual paths from controller configuration:

```text
<revision_workspace>/
  baseline_receipt.json
  working_state.json
  batch_index.tsv
  cumulative_consequences.tsv
  controller/
    baseline_snapshot.tsv
    manuscript_dependency_graph.tsv
  renderer/
  batches/<batch_id>/
    batch.json
    prompts/
      0001.md
      0002_clarification.md
    decisions/
      0001_contract_disposition_<state>.md
    checks/
    scope.tsv
    files/<root-relative edited paths>
    generated/<root-relative generated paths>
    artifact_changes.tsv
    artifact_classification.tsv
    contract_sensitive_changes.tsv
    hash_backreferences.tsv
    consequence_audit.tsv
    consequence_audit_validation.json
    consequence_delta.tsv
    notes.md
  candidates/<candidate_id>/
    composition.json
    assembly/
    validation/
  previews/<preview_id>/
    manuscript_preview.html
    render_receipt.json
    active_batches.tsv
    open_consequences.tsv
```

The workspace owns a versioned controller copy and renderer execution contract.
The renderer entrypoint itself must be part of the authenticated baseline or a
subsequent recorded overlay, so an active revision does not depend on later
skill changes or construction scaffolding.

## Revision batches

Treat one coherent user request as one batch even when it contains multiple
modifications or artifact owners. Preserve the prompt verbatim and append
clarifications rather than rewriting it. Split batches only when modifications
need independent acceptance or abandonment, may block independently, or have a
separate user-review checkpoint.

Register a batch before its first write. Record its baseline, predecessor
batches, prompt identities, scope, state, and assembly-contract disposition.
Record individual modified files and generated artifacts in
`artifact_changes.tsv`; do not create per-edit directories.

Use `preserved` unless the user explicitly approves another contract
disposition. A non-preserved starting disposition or later change requires a
verbatim approval file stored under the batch. Agent inference is not approval.
Contract-sensitive entries name exact control-plane entrypoints and records;
they do not include scientific constructors, package-local generators, or
ordinary content manifests solely because of directory placement.

## Consequences

Freeze the bundled manuscript dependency graph in the workspace. The graph
names general artifact classes and their downstream relationships. Do not
encode individual request types, figure identifiers, project directory names,
or prescribed semantic decisions in either the controller or the graph.

Compare effective artifacts with the parent candidate and record every added,
modified, or deleted artifact. Independently scan effective-parent artifacts
for prior SHA-256 identities of changed files and record every match. Surface
the changed artifacts and hash consumers as the complete input set for the
batch-local classification.

The agent, not the controller, classifies the observed paths and assesses
validity. Assign every changed artifact and exact-hash consumer at least one
graph class with a short rationale; multiple classes are allowed. The special
class `terminal` may be used only with a rationale explaining why that observed
item has no downstream manuscript consequence. The classification is local to
the batch and does not create a persistent path-to-class registry.

Start from every nonterminal class assigned to a changed artifact and traverse
each downstream graph edge. Record one audit row per visited edge with affected
paths, owner, a short rationale, and one decision:

- `invalidated`: the consumer cannot presently be relied upon; continue from
  its class;
- `remains_valid`: the consumer remains usable; terminate that path; or
- `unresolved`: assessment is blocked; terminate that path and record why.

Treat exact-hash consumers as factual dependencies that must appear in the
classification and, unless terminal, in the audit traversal. Deduplicate a
repeated class edge, but do not merge different incoming edges or silently
prune a path. The deterministic validator checks classification coverage,
graph closure, exact-hash coverage, decision/state consistency, and termination
rationales. It never supplies a class or decision. Marking a consumer invalid
does not mandate rewriting; its recorded owner may update it, authenticate it
unchanged, remove it, accept an exception, or block. A batch may preview with a
pending audit, but it cannot become complete without a passing audit.

## Candidates and previews

Materialize candidates only for batch completion, explicit checkpoints,
on-request rendering, or final reconciliation. Compose batches in declared
order, refuse unresolved overlapping writes, and use unchanged baseline
artifacts as fallbacks. Candidates are composition and render targets. Run
mechanical generation from staged batch packages or disposable scratch space;
do not rewrite general assembly machinery merely to accommodate candidate
relocation.

Before deleting disposable execution space, register a batch-local receipt for
every replay or validation result used to justify the edit. Preserve enough
command, input, job/exit, and observed-output identity to audit the claim later;
narrative notes alone are not an execution receipt.

Render despite open consequences or a pending consequence audit; fail only for
a mechanical render error. Mark the HTML outside the manuscript body as an
unsealed preview. Record whether the latest preview matches the working state.
Its receipt must identify the baseline, active batches, effective inputs,
renderer and configuration, output hash, timestamp, open consequences, and
pending audits.
After a successful first preview, pause for review unless the user explicitly
requested reconciliation, promotion, or resealing.

Treat changes to configured contract-sensitive control-plane entrypoints or
records as scope expansion. A batch whose contract disposition is `preserved`
must block until the user approves a new disposition. Terminal validation is
candidate-specific: its receipt must name the materialized candidate assembly
as its target. A baseline-targeted report does not validate a candidate.

## Reconciliation

Resolve consequences through their owning workflows by updating an artifact,
authenticating reuse, recording reviewed no-change, accepting an exception, or
blocking. Route a fully reconciled candidate to terminal manuscript assembly
and validation. Preserve the prior seal and every versioned preview.
