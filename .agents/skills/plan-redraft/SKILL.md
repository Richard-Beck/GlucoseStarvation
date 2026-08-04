---
name: plan-redraft
description: "Configure and route deterministic manuscript revision runs. Use when Codex needs to choose between a new manuscript assembly and staged in-place amendments to a sealed assembly, initialize or continue revision work, register a revision batch, render a revision preview, inspect presumed consequences, or route unresolved consequences to their owners."
---

# Plan Redraft

## Purpose

Act as the thin user-facing wrapper around the deterministic in-place revision
workspace manager bundled with this skill. Select a revision mode, configure
the manager when needed, assess consequences against the bundled manuscript
dependency graph, report revision state, and route invalidated artifacts to
their owners. The manager records facts and validates audit completeness; it
does not decide whether manuscript content remains valid.

For in-place work, read `references/staged_revision_contract.md`; do not repeat
its layout or state model in planner prose. Do not perform stage-owner work
merely because it appears in a consequence ledger.

## Required Inputs

- A manuscript target.
- The user's requested edit or revision objective.
- Optional existing revision configuration or controller state.

Inspect the target before choosing a mode. Ask the user only when the available
state cannot distinguish the modes safely.

## Modes

Choose exactly one:

- `new_root`: build toward a new versioned manuscript assembly. Use when no
  authenticated sealed assembly is available, the user requests a separate
  release, or the intended final package cannot retain the target assembly's
  contract.
- `in_place`: retain an authenticated sealed assembly as an immutable baseline
  and stage edits in its configured revision workspace. A possible contract
  change may be explored in place, but it remains an unresolved consequence
  until final disposition is chosen.

A narrow edit is not by itself a reason to choose `new_root`. In-place work
means staging inside the target root, not overwriting its sealed files.

## Configuration Ownership

For `in_place`, write a project-local controller configuration from
`assets/controller_config.template.json`. Use
`references/manuscript_dependency_graph.tsv` as the default graph. Copy and
extend it only when the manuscript contains a genuinely additional artifact
class; never add graph branches for individual request types. Use
`assets/artifact_classification.template.tsv` only after `record`, for the
paths actually observed in that batch. Record at least:

- mode and target identity;
- baseline seal identity for `in_place`;
- a top-level revision-workspace location within the target assembly;
- exact user instruction or managed-feedback handle when explicitly supplied;
- batch identifier and ordered active batches;
- preview policy, defaulting to `batch_complete` and permitting explicit
  checkpoints or on-request rendering;
- assembly-contract disposition: `preserved`, `possibly_changed`, or
  `replacement_required`; and
- exact project-local paths identifying control-plane entrypoints or
  records whose modification changes assembly, rendering, validation, or seal
  semantics. Do not classify scientific constructors, package-local generation
  scripts, or content-specific manifests as contract-sensitive merely because
  of their location.

Use paths returned by the controller or its configuration. Do not hard-code
project directory names in this skill.

## In-place Manager

Use `scripts/revision_workspace.py`. Initialization, editing, and previewing are
separate operations. The manager performs mechanical state changes; Codex and
the relevant artifact owners perform the edits.

For `in_place` work:

1. Run `authenticate --config CONFIG` when inspecting a prospective baseline,
   or `init --config CONFIG` to authenticate it and initialize the workspace.
   Authentication verifies the existing seal identity only. It must not
   rebuild, render, or run terminal validation.
2. Use the controller copy placed in `<revision_workspace>/controller/` for all
   subsequent commands so an active revision is insulated from later skill
   changes.
3. Run `begin-batch` before the first modification. Preserve the
   originating user prompt verbatim and append later clarifications.
4. Run `stage-path` before editing an existing assembly artifact. Give artifact
   owners the returned batch-overlay path as their only manuscript-facing write
   target. Use `register-generated` for a generated file and `delete-path` for
   an intended deletion. Do not create a directory for every atomic edit.
5. Run `record` after the requested edits. It compares the batch overlay with
   the effective parent and records actual additions, modifications, deletions,
   and exact-hash backreferences. It does not infer semantic consequences.
6. Read the workspace's frozen `controller/manuscript_dependency_graph.tsv`
   completely. Classify every changed artifact and every exact-hash consumer
   reported by `record` in a batch-local `artifact_classification.tsv`. Give
   each observed path/source at least one graph class and a short rationale;
   multiple classes are allowed. Use the special class `terminal` only when an
   observed item intentionally has no downstream manuscript consequence, and
   explain why. Do not build or maintain a repository-wide path registry.
   Starting from every nonterminal class assigned to a changed artifact,
   traverse each downstream edge.
   Write one row in `assets/consequence_audit.template.tsv` for every visited
   edge:
   - `invalidated`: give a short reason and continue from the downstream class;
   - `remains_valid`: give a short reason and stop that path; or
   - `unresolved`: give the blocker and stop that path.
   Include every exact-hash consumer reported by `record`. Deduplicate the same
   class edge while retaining affected artifact paths. Multiple independent
   incoming edges remain separate decisions. Do not silently terminate a path.
   Run `register-consequence-audit` with both the classification and audit
   files; repair them until the closure validator passes. This assessment has
   no intrinsic subagent requirement.
7. Run `materialize` to compose an immutable-baseline-plus-overlay candidate.
   Treat it as a composition and rendering target. Perform mechanical
   generation in the batch's staged package or disposable scratch space, then
   register generated assembly artifacts and record once. Do not modify general
   rebuild wrappers or validators merely to make them run from a nested
   candidate.
   Before deleting scratch space, use `register-check` to preserve any replay,
   validation, or execution receipt cited as evidence. Include the invoked
   entrypoint and input identity, exit status or scheduler job identity, and
   observed output hashes; retain stdout/stderr when they carry material facts.
8. Run `preview` to invoke the configured renderer on a materialized candidate
   at batch completion, an explicit user-review checkpoint, on request, or
   before final reconciliation. Rendering must not wait for consequence
   resolution. Normally create one candidate and one preview for a coherent
   request; intermediate work belongs in scratch space. Create another
   candidate only after a user-requested change, a genuine render failure, or
   an explicitly recorded reason.
9. Run `status` to report whether the latest preview is current with the working
   state. Report pending audits, open `invalidated` or `unresolved` decisions,
   and stale recorded batches.
10. Route only the selected next consequence to its recorded owner. An owner may resolve
   a consequence by changing an artifact, authenticating reuse, recording no
   change after review, or identifying a blocker.
11. Route a reconciled candidate to manuscript assembly for terminal packaging
    and validation. Plan Redraft preserves the original seal and does not perform
    terminal validation itself.

Unless the user explicitly requested reconciliation, promotion, or resealing,
stop after the first current preview and report open consequences for review.
Do not resolve them, update status/checksum/validation artifacts, or continue
into terminal packaging on the user's behalf.

If a batch declared `preserved` changes configured contract-sensitive
machinery, `record` must block. Report the paths and obtain an explicit
disposition before continuing. Supply the user's verbatim approval file to the
controller; never synthesize approval from agent judgment. Do not normalize a
narrow artifact edit into a general renderer, rebuild, validator, status,
checksum, or seal rewrite.

If a controller command blocks, report the exact failed invariant. Do not
emulate it with direct writes into the sealed assembly.

For `new_root`, do not invoke the in-place manager. Configure the new-root
destination and route construction to the first required artifact owner. Do
not recreate a universal stage sequence; route from declared inputs,
consequences, and user choices.

## Preview Semantics

A preview proves only that the composed candidate rendered. It may deliberately
contain baseline consumers whose consequences remain open. Report its exact
path, active batches, whether it reflects the current working state, renderer
receipt, and unresolved consequence count. Never describe an unsealed preview
as valid, passing, or release-ready.

A terminal validation receipt applies to a candidate only when the recorded
validation target is that candidate's materialized assembly root. A report
targeting the immutable baseline cannot validate a candidate, even when it was
produced while working inside the candidate directory.

## Delegation

This routing workflow has no intrinsic subagent requirement. Artifact owners
retain any independent-inspection or fresh-context requirements imposed by
their own contracts.

## Boundaries

- Never modify any sealed-baseline file, including paths intentionally excluded
  from its canonical checksum inventory. The controller authenticates their
  content and metadata for the lifetime of the workspace.
- Never turn a candidate failure into a change to baseline status.
- Never infer that invalidation requires rewriting. The owner may authenticate
  unchanged content after review.
- Never let the controller or an artifact classification make a
  semantic-validity decision.
- Never terminate a graph path without a recorded `remains_valid` or
  `unresolved` rationale.
- Never use ordinary current-turn instructions as managed-feedback state unless
  managed-feedback mode was explicitly enabled.
- Never run terminal validation merely to update the working-state ledger.
- Never run the preview renderer during baseline authentication or workspace
  initialization.

## Completion

Finish a planner turn by reporting:

- selected mode and configuration path;
- baseline or new-root target;
- active batch and candidate identifiers;
- latest preview HTML and whether it is current, when a preview exists;
- consequence-audit state and open decisions grouped by decision;
- selected next owner or user decision; and
- any missing controller capability or blocker.
