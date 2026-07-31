# Redraft Assembly Gates

Run these gates only at a Manuscript Assembly boundary, not on every planner
turn.

## User-Decision Ledger Check

Read the current blocker ledger before the compatibility and feedback gates.
Do not route Assembly while any record is `awaiting_user`. Route every
`deferred_to_assembly` record and its authenticated artifacts to Assembly and
require package `WARN` plus `review_state/deferred_user_decisions.md`. Resolved
records remain audit history and do not warn. A missing ledger is valid only
when every accepted stage return had an empty reconciled blocker sidecar.

## Pre-Assembly Compatibility Gate

Run after all seven pre-Assembly stages have current accepted return bundles and
whenever later replay reaches the Assembly boundary again.

1. Require each return to contain `compatibility_sidecar.json` conforming to
   `references/compatibility_sidecar_contract.md`.
2. Write
   `<redraft_root>/compatibility_gate/<turn_id>_sidecar_registry.json` with the
   seven current consumer/sidecar paths. Include accepted fact exceptions only
   when explicitly authorized by the user or active run policy.
3. Run `scripts/check_compatibility_sidecars.py` and write its JSON and Markdown
   reports beside the registry.
4. Treat `PASS` and `WARN` as passing. Route every `WARN` exception to Assembly
   review state. On `BLOCKED`, route the report and relevant sidecars to the
   earliest workflow consumer named by the mechanical report. The planner does
   not decide which scientific assertion is correct.
5. Re-run the ordinary workflow suffix and this gate after repair.

Do not route a first Manuscript Assembly handoff while this gate is blocked.

## Feedback-Completeness Gate

The ordered feedback-consumer sequence is:

```text
analysis
manuscript-figure-workflow
claim-graph-integration
results-text
method-table-provenance
manuscript-legend-writing
serve-manuscript-abstract-introduction-discussion
manuscript-assembly
```

For a first Assembly handoff, audit the first seven consumers. Once Assembly has
received a handoff or returned, audit all eight. Audit every tracking root in the
active plan using only mechanical feedback-manager state.

For each required consumer, verify every span and every active item has a current
`signed_off` response. Missing, `needs_followup`, and `reopened` responses are
open. Inactive items are not open. Skips, accepted exceptions, no-change, and
not-applicable dispositions pass only with a current feedback-manager signoff.

Write:

```text
<redraft_root>/feedback_gate/<turn_id>_feedback_completeness.md
```

List each tracking root and consumer, unsigned-span and active-unsigned-item
counts, and opaque missing IDs. Do not interpret feedback content.

If the gate fails:

1. Route to the earliest consumer missing any signoff with the raw feedback
   handles and gate audit.
2. Mark that stage and every later stage through Assembly pending replay.
3. Advance through every later consumer in order, requiring a current
   compatibility/no-change or changed-output return even if it previously
   returned a conforming bundle.
4. Re-run both applicable Assembly-boundary gates after the workflow suffix.

Do not route Assembly, accept its package as terminal, or declare completion
while feedback is open.

After Assembly returns and the full eight-consumer gate passes, write
`<redraft_root>/redraft_complete.json` with status `complete`, the assembly root,
package status path, passing feedback-audit path, and completion timestamp. Do
not write it for the first Assembly handoff or a provisional package.
