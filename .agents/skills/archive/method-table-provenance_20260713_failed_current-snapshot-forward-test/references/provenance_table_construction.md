# Provenance Table Construction

Use this workflow to construct, revise, or merge Methods-facing provenance for
current manuscript panels. Return one canonical current-snapshot table, not
Methods prose or a history of prior table states.

## Contents

- [Output contract](#output-contract)
- [Current-snapshot boundary](#current-snapshot-boundary)
- [Row semantics](#row-semantics)
- [Lock contract](#lock-contract)
- [Tracing workflow](#tracing-workflow)
- [Multi-panel reconciliation](#multi-panel-reconciliation)

## Output Contract

Use exactly these columns:

```text
id | parent | what | why | comment | lock_target | lock_selector | lock_kind | sha256 | hash_status
```

Use one row per current manifest panel as a structural graph leaf. Give it:

- an id of `<figure>#panel_<panel>`, matching `target_figure_set.tsv`;
- the current final whole-figure path in `lock_target`;
- the mandatory `#panel_<panel>` selector;
- `lock_kind = file`;
- the whole-figure SHA256;
- `hash_status = computed_self`;
- its last Methods-relevant object or transformation as `parent`.

Panels in the same final figure intentionally share `lock_target` and `sha256`.
The selector distinguishes their logical endpoints.

## Current-Snapshot Boundary

Include only the current panel-rooted provenance graph. Every row must be
upstream-reachable from a current panel leaf.

Do not include nodes whose only purpose is:

- plotting, rendering, visual styling, or rasterization;
- panel extraction, layout, composition, or figure integration;
- file copying, packaging, version overlays, or workflow history.

Collapse those display operations into the panel leaf. Preserve superseded
tables or process history outside the canonical Methods bundle when needed.

## Row Semantics

Use `id` for one durable object, named in-memory object, or Methods-relevant
transformation. Prefer stable object names and canonical repo-relative paths.

Use `parent` for immediate upstream inputs. Separate multiple parent ids with
semicolons and backtick ids or paths; do not put prose in this cell. Multiple
parents are valid only when the row directly combines those inputs.

Use `what` for the verified operation and `why` for its analytical role. Use
`comment` for source paths, functions, filters, scopes, thresholds, accepted
exceptions, or the reason a trace terminates.

Prefer boundaries at saved files, named data/model objects, consequential
filters or summaries, joins across analytical streams, model fits, and terminal
raw or external inputs. Do not split one Methods-relevant transformation into
rows that only expose implementation plumbing.

## Lock Contract

Use these `lock_kind` values:

```text
file | code | proxy_file | representative_file | terminal | external | unresolved
```

Use these hashed `hash_status` values:

```text
computed_self
computed_downstream_proxy
computed_upstream_proxy
computed_code_proxy
computed_representative
metadata_checksum
```

Use these un-hashed values where appropriate:

```text
not_applicable | external | missing | ambiguous | unresolved
```

For a durable file or code object, hash the object bytes and use
`computed_self`. For a selected fragment, hash the whole file and record the
selector; do not invent selector-specific hashes.

For an in-memory object or transformation, use the most informative adjacent
durable proxy, in this order:

1. a directly generated downstream artifact;
2. a specific durable upstream input;
3. the generating code or config.

Use `proxy_file` with the corresponding upstream/downstream status, or `code`
with `computed_code_proxy`. Explain non-obvious proxy choices in `comment`.

For a many-file scope, use an honest representative content file and explain
the represented scope, or split the row. Do not substitute a manifest for file
contents unless the row represents that manifest. Do not use recursive
directory hashes by default.

For raw, external, missing, or unresolved boundaries without an honest target,
use `sha256 = NA` and the matching terminal/external/unresolved lock metadata.
Explain why tracing stops.

Repeated proxy, code, representative, or panel-leaf hashes are valid. By
contrast, two `computed_self` rows with the same `lock_target` and selector
claim to be the same durable object and must be reconciled.

## Tracing Workflow

For each current panel:

1. Confirm its exact figure id, panel id, final figure path, and whole-figure
   hash from the current manifest.
2. Create the panel leaf and identify its last Methods-relevant parent from
   panel metadata, semantic notes, and generating code.
3. Trace upstream one verified edge at a time through named analytical objects
   and transformations.
4. Narrow figure-wide source declarations to the inputs actually used by the
   panel.
5. Continue to raw, external, or explicitly unresolved boundaries.
6. Compute locks and return panel-scoped rows for reconciliation.

Read structured metadata and scripts before prose when resolving an edge.
Represent missing support as a terminal gap rather than inventing a parent.

## Multi-Panel Reconciliation

Panel-scoped traces may be temporary working artifacts. Reconcile them into one
replacement current-snapshot table:

1. Normalize ids and current paths.
2. Merge rows only when they represent the same durable object or action and
   their immediate parents are compatible.
3. Preserve distinct logical transformations even when they share a proxy hash.
4. Reconcile aliases in `comment`; do not keep duplicate alias rows.
5. Treat conflicting hashes for one lock target as an error.
6. Remove any row not upstream-reachable from a current panel leaf.
7. Validate the integrated table before treating it as canonical.

Do not append a new overlay to an older canonical table. Existing tables are
evidence to inspect, not historical components that must survive in the new
snapshot.
