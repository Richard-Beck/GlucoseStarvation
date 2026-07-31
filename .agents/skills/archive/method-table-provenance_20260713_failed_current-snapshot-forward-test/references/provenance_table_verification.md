# Provenance Table Verification

Use this workflow before reusing, classifying, outlining, or drafting from an
existing provenance table. Verification has two mandatory passes:

1. structural panel-leaf and current-scope validation;
2. SHA256 lock-target validation.

Do not rewrite the table during a verification-only request.

## Inputs

Require one canonical Markdown table with exactly these columns:

```text
id | parent | what | why | comment | lock_target | lock_selector | lock_kind | sha256 | hash_status
```

Materialize a panel-level `target_figure_set.tsv` from the current figure
manifest:

```text
endpoint_id	artifact_path	lock_selector	sha256
F2#panel_a	path/to/F2.png	#panel_a	<whole-figure-sha256>
F2#panel_b	path/to/F2.png	#panel_b	<same-whole-figure-sha256>
```

Use canonical repo-relative final-figure paths. Require one row per current
panel, and require each endpoint id to end with its selector.

## Structural Panel-Leaf Pass

Run:

```bash
python3 .agents/skills/method-table-provenance/scripts/verify_manuscript_endpoints.py \
  path/to/target_figure_set.tsv \
  path/to/locked_provenance_table.md \
  --root . \
  --output path/to/manuscript_endpoint_verification.md
```

The validator derives endpoints from structural graph leaves. It fails unless:

- the document contains exactly one canonical 10-column provenance table and
  no legacy three-column endpoint table;
- ids are unique, all parents resolve, and the graph is acyclic;
- graph leaves and target manifest endpoint ids match exactly;
- each leaf has a Methods-relevant parent, mandatory selector, current final
  whole-figure path, and matching whole-figure hash;
- the number of unique leaf hashes equals the number of current final figures;
- every row is upstream-reachable from a current panel leaf;
- lock kinds, statuses, hashes, and targets are coherent;
- duplicate `computed_self` locks for the same target and selector have been
  reconciled.

Repeated hashes remain valid for panels in one figure and for distinct nodes
using shared proxy, code, or representative locks.

Any nonzero result is a hard structural failure. Do not reuse downstream
classifications, spines, or prose until the current-snapshot graph passes.

## Lock-Target Pass

After structural validation, run:

```bash
python3 .agents/skills/method-table-provenance/scripts/verify_provenance_locks.py \
  path/to/locked_provenance_table.md \
  --root . \
  --output path/to/provenance_lock_verification.md \
  --details-tsv path/to/provenance_lock_verification.tsv \
  --fail-on-problems
```

This recomputes every stored hash. With `--fail-on-problems`, parse errors and
changed, missing, ambiguous, unresolved, or invalid hash checks fail. Honest
terminal, external, missing, ambiguous, or unresolved gap rows with an
un-hashed status remain explicit but do not fail solely for lacking a digest.

Proxy and representative matches prove only that their named lock targets are
unchanged; they do not prove that an unsaved in-memory object or an entire
many-file collection is unchanged.

## Failure Handling

For each failed panel endpoint, traverse its existing parents and use the lock
report to find the nearest changed, missing, ambiguous, or unresolved row. If
locks are unchanged but leaf scope differs, report a current-scope mismatch.

Repair or reconstruct only affected current provenance, but return one complete
replacement snapshot. Do not append a repair overlay or restore disconnected
historical components.

Use the two emitted reports without manually changing their determinations.
Record unresolved semantic trace gaps separately.
