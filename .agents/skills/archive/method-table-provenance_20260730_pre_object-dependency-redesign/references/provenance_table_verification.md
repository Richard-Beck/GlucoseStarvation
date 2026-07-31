# Provenance Table Verification Workflow

Use this optional workflow to check whether an existing locked
method-table provenance table remains valid against current files or against a
newly supplied manuscript figure/panel artifact package.

This workflow checks for two different invalidation modes:

1. Locked artifacts silently changed, moved, or disappeared on disk.
2. The supplied manuscript figure/panel scope no longer matches the extant
   provenance table, even if the table's previously locked artifacts did not
   change.

Check both modes deterministically at the durable-file boundary. The lock-target
verifier checks the existing graph artifacts, while the endpoint verifier checks
that the declared manuscript-visible panels are exactly the structural leaves
of the provenance graph.

## Inputs

Required:

- An existing provenance table with columns:

```text
id | parent | what | why | comment | lock_target | lock_selector | lock_kind | sha256 | hash_status
```

Optional, but required for scope checks:

- A current manuscript figure/panel artifact package with enough metadata to
  identify manuscript-visible endpoints. This may be a figure-set manifest,
  final image directory, subpanel inventory, package provenance files, semantic
  notes, or equivalent current artifact bundle.

## Required Outputs

Produce:

1. `provenance_lock_verification.md` from the lock-target verifier.
2. `manuscript_endpoint_verification.md` from the endpoint verifier for a
   Methods revision.
3. A targeted upstream escalation report for any flagged endpoints.
4. Unresolved or ambiguous items.

Do not rewrite the provenance table unless the user explicitly asks for a
revision. Do not draft manuscript prose.

## Always Check Locked Artifacts

Use the bundled verifier script for the first invalidation mode:

```bash
python3 .agents/skills/method-table-provenance/scripts/verify_provenance_locks.py \
  path/to/locked_method_table.md \
  --output path/to/provenance_lock_verification.md \
  --details-tsv path/to/provenance_lock_verification.tsv
```

Run this script for every verification task. It parses the locked provenance
table, recomputes SHA256 for rows with stored lock hashes, and reports changed,
missing, ambiguous, unresolved, and unchecked rows.

For rows with `lock_kind = proxy_file`, `representative_file`, or `code`, report
the status as a lock-target status, not as proof that the represented in-memory
object or many-file run set is fully unchanged.

## Endpoint Scope Check

For every Methods revision, first materialize the manuscript-visible endpoints
from the supplied figure/panel package as the Methods-local
`target_figure_set.tsv`, with exactly these tab-separated columns:

```text
endpoint_id | artifact_path | lock_selector | sha256
```

Use one row per panel, with an endpoint id such as `F2#panel_a`, selector
`#panel_a`, and the current final whole-figure path and digest. Then run:

```bash
python3 .agents/skills/method-table-provenance/scripts/verify_manuscript_endpoints.py \
  path/to/methods_root/target_figure_set.tsv \
  path/to/locked_provenance_table.md \
  --root . \
  --output path/to/manuscript_endpoint_verification.md
```

Run this before reusing existing classifications, Methods spines, or Methods
text. The verifier requires exactly one canonical 10-column provenance table,
derives its leaves from `parent` links, and enforces exact leaf/manifest ids,
panel selectors, whole-figure paths and hashes, parent resolution, acyclicity,
unique parent references within each row, current-leaf reachability, and
reconciliation of duplicate direct-file object locks. Shared panel-figure and
proxy hashes remain valid. Do not add or infer a separate endpoint mapping
table.

A nonzero result is a hard provenance verification failure and cannot be
accepted as an exception. Repair or reconstruct the affected provenance table,
then rerun the verifier before Methods validation.

## Targeted Upstream Escalation

Always escalate upstream for every nonmatching endpoint reported by the
verifier.

For each flagged endpoint:

1. Traverse the existing provenance graph upstream through `parent` links.
2. Reuse the mandatory lock-target hash-check results for upstream rows.
3. Identify the nearest upstream rows with `changed`, `missing`, `ambiguous`, or
   `unresolved` lock-target status.
4. If upstream locked artifacts are unchanged but endpoint scope changed, report
   the issue as a scope mismatch rather than a silent artifact-drift problem.

Do not fully retrace new upstream metadata by default. Reconstruct provenance
for each endpoint that fails deterministic validation.

## Subagents

Do not use subagents for the deterministic hash pass.

Do not delegate or manually interpret endpoint mapping. The deterministic
verifier is authoritative. Subagents may be used only to reconstruct provenance
after the verifier identifies affected endpoints.

## Report Format

Use the reports emitted by the bundled verifiers without manually rewriting
their status or endpoint mappings. Add targeted upstream escalation and
unresolved-item sections separately when a verifier reports a failure.
