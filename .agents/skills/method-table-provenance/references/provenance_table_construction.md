# Recursive Object-Provenance Construction

Use this workflow to build or revise the canonical locked provenance table from
current figure-generating code and analytical objects.

Read `object_resolution.md` before launching Object Resolvers.

## Construction mode

Apply one mode selected under the main skill contract:

- `targeted_update`: initialize the registry with authenticated unaffected rows
  plus only affected current roots and invalid/new objects queued for
  resolution. Preserve retained rows only when current locks, identities, and
  branch connections are authenticated.
- `clean_reconstruction`: initialize from all deduplicated current
  figure-generating scripts and known generating objects. Do not import prior
  graph rows.

`verified_reuse` does not enter this construction workflow.

## Outputs

The canonical output remains one table:

```text
id | parent | what | why | comment | lock_target | lock_selector | lock_kind | sha256 | hash_status
```

Use these construction artifacts while building it:

```text
figure_object_roots.tsv
object_registry.tsv
dependency_edges.tsv
object_resolution_index.tsv
object_fact_cards/
unresolved_methods_issues.md
```

### Figure/object roots

Write `figure_object_roots.tsv` with:

```text
endpoint_id | figure_id | generator_path | generator_selector | generated_object | evidence
```

Allow multiple rows for an endpoint when multiple scripts or objects directly
participate. Deduplicate generator paths before worker assignment. Preserve the
endpoint mapping even when one script generates several figures.

### Object registry

Write `object_registry.tsv` with:

```text
object_id | object_type | canonical_path | selector | sha256 | resolution_status | transformation | methods_role | evidence
```

Use these object types:

```text
panel_endpoint
script
function
action
named_object
data_file
config
raw_source
external
unresolved
```

Use these resolution states:

```text
queued
resolved
terminal
unresolved
```

### Dependency edges

Write `dependency_edges.tsv` with:

```text
child_id | parent_id | relation | dependency_status | evidence | comment
```

Use these dependency states:

```text
candidate
confirmed
unresolved
```

Relation names may be concise and object-appropriate, such as `generated_by`,
`consumes`, `sources`, `configured_by`, `calls`, or `selects_from`.

## Identity and deduplication

Deduplicate before assigning work and again when each fact card returns.

Use:

- canonical repo-relative path plus selector for durable or selected objects;
- canonical script path plus function/object name for in-memory actions and
  named objects; and
- a stable explicit identifier for external or terminal objects.

Do not merge objects by SHA256 alone. Identical bytes at different paths can
have distinct provenance roles. Treat copied artifacts as aliases only when
structured evidence identifies one canonical source.

Do not collapse distinct selected objects inside one file. Hash the whole file
and retain selectors.

## Closure algorithm

1. For clean reconstruction, register all deduplicated figure-generating
   scripts and known generating objects as `queued`. For targeted update,
   register only affected current roots and invalid/new objects as `queued` and
   retain authenticated unaffected objects as resolved.
2. Assign each canonical queued object to one Object Resolver.
3. Validate the returned fact card against the assigned identity and controlled
   object/status vocabularies. Reject contract-invalid cards rather than
   silently normalizing invented types.
4. Register its confirmed local transformation and direct dependencies.
5. Omit loaded or sourced but unused, display-only, administrative, runtime, or
   implementation-only candidates.
6. Add new unresolved direct parents to the queue after deduplication.
7. Repeat until no queued object remains.
8. Preserve explicit `unresolved` and legitimate `terminal` boundaries.

Do not let a worker recurse. The router must see and deduplicate every newly
discovered object before another assignment.

## Script objects

For scripts or functions, enumerate direct scientific:

- sourced local scripts;
- data reads and object loads;
- configs, manifests, command-line inputs, and dynamically resolved paths;
- called project functions that perform Methods-relevant work; and
- saved analytical outputs.

Distinguish:

- declared or loaded dependencies;
- dependencies actually consumed by the assigned object;
- display-only dependencies; and
- unused dependencies.

Return and register only consumed, candidate, or unresolved scientific
dependencies. Do not inventory rejected dependencies in Object Resolver fact
cards.

Do not register or return languages, runtimes, environment managers, installed
packages, libraries, namespaces, shell commands, system executables, or generic
I/O, serialization, table-manipulation, plotting, and layout tools. Preserve
their scientifically meaningful effect as an action or transformation when
needed. Reserve `external` for genuine external scientific sources or opaque
methodological boundaries, not software.

Apply the same filter to project-local code. Omit path resolvers, config readers,
argument handlers, orchestration wrappers, logging, output writers, and
execution-only validation/status gates unless they themselves define scientific
values, populations, transformations, or choices. Retain the scientific config,
resolved model code, or generated object rather than the implementation helper
that locates or loads it.

Software identities and versions belong to a separately requested environment
inventory or later Methods reporting, never the recursive object graph by
default.

## Data objects

For a durable or named data object:

1. Establish its exact path, selector, and contents relevant to the consumer.
2. Prefer structured provenance, execution records, and actual write calls when
   identifying its generating code.
3. Identify and enqueue the immediate generating code or action.
4. Let the later resolver for that code or action identify its own direct
   inputs; do not flatten those inputs into direct parents of the data object.
5. Preserve any embedded upstream paths as candidate hints until that generator
   is resolved.
6. Stop explicitly at raw, external, manually curated, or unresolved
   boundaries when no defensible generator exists.

A filename mention is not evidence that a script generated the object.

## Semantic resolution

File closure is construction evidence, not the canonical Methods graph.

Represent a transformation when it changes or defines scientific values,
populations, comparisons, or uncertainty, including:

- filtering or inclusion/exclusion;
- calibration or normalization;
- joins that define analytical populations;
- summaries and endpoint construction;
- model fitting, selection, inference, or prediction;
- simulation protocols;
- thresholds and classification;
- variable-to-axis, group-to-facet, or other scientifically meaningful display
  mappings.

Omit pure layout, typography, palette, copying, integration, serialization,
and raster rendering.

Use one node per durable object, named in-memory object, or Methods-relevant
action boundary. Do not combine distinct transformations merely because they
occur in one function.

## Canonical reconciliation

The router converts resolved registry content into one current-snapshot graph.

- Include confirmed Methods-relevant dependencies.
- Retain unresolved boundaries honestly.
- Exclude runtime, package, library, shell-tool, and generic implementation
  dependencies even if a worker returned them.
- Exclude implementation-only project helpers and execution-only validation
  records even if they are directly called.
- Reconcile aliases by canonical identity and evidence.
- Keep conflicts visible until resolved.
- Attach each current panel leaf to its last Methods-relevant object or action,
  not automatically to a whole figure script.
- Prune every row not upstream-reachable from a current panel leaf.

In targeted update, the prior locked table contributes only rows independently
authenticated as unaffected. In clean reconstruction it is not a merger input;
compare it with the new table only after current reconciliation when a
regression or change report is useful.

## Lock fields

For durable files:

- use a canonical repo-relative `lock_target`;
- hash the complete file;
- use a selector only for a named fragment or selected object;
- set `lock_kind` to `file` or `code`; and
- set `hash_status` to `computed_self`.

For in-memory objects or actions, prefer a direct downstream artifact, then a
specific upstream input, then generating code as an honest proxy. Record the
proxy direction in `hash_status` and explain it in `comment`.

For many-file objects, use an honest representative content file only when it
can reveal relevant drift. Otherwise split the node or mark the lock
unresolved. Do not substitute a manifest hash for content it does not contain.

Use `NA` only for legitimate external, terminal, missing, ambiguous, or
unresolved boundaries, with an explicit status and explanation.

## Audit and validation

Partition the reconciled table into non-overlapping batches of at most 75 rows.
Auditors receive assigned rows plus direct-parent definitions and report
defects only:

- unsupported generator or immediate-parent claims;
- missing Methods-relevant transformations;
- loaded dependencies incorrectly treated as consumed;
- duplicate identities or aliases;
- display-only nodes in the canonical graph;
- premature stopping at derived reports or manifests; and
- unresolved conflicts hidden as certainty.

Resolve reported defects centrally. Then:

1. Validate `object_registry.tsv` and `dependency_edges.tsv` with
   `scripts/validate_object_registry.py --require-closed`.
2. Prune unreachable canonical rows.
3. Run `scripts/verify_provenance_locks.py`.
4. Run `scripts/verify_manuscript_endpoints.py`.

Do not draft Methods prose during construction.
