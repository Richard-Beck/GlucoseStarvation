# Object Resolver Contract

Use this contract for one canonical object at a time.

## Worker boundaries

Give the worker:

- one assigned object id and type;
- its canonical path and selector when known;
- the current figure/package metadata that names it;
- candidate direct dependencies when mechanically discoverable; and
- any bounded structured provenance needed to interpret the object.

Do not give the worker the prior canonical provenance table. Do not ask it to
reconcile global duplicates.

The worker may inspect the assigned object, its defining code, and only the
additional files needed to resolve its direct dependencies. It must not recurse
into newly discovered parent objects.

## Required reasoning

For the assigned object:

1. Verify identity, path, selector, and SHA256 when directly hashable.
2. Enumerate direct scientific data, project-code, config, action, and source
   candidates.
3. Determine which candidates actually contribute to the object.
4. Identify the immediate generating code or action for data objects.
5. Describe the local scientific transformation and Methods role.
6. Omit unused, display-only, administrative, runtime, and implementation-only
   candidates from the return.
7. Return uncertain relationships as unresolved; do not guess.
8. List newly discovered direct parent objects for router-controlled queuing.

Use only these object types everywhere a fact card requests `object_type` or
`dependency_type`:

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

Do not invent narrower synonyms such as `code_file`, `data_table`,
`validation_record`, `R_package`, or `file`. Put the narrower description in
`relation` or `comment`.

## Scientific provenance boundary

Return scientific lineage, not the execution environment.

Do not return runtime or implementation machinery anywhere in the fact card:

- languages, interpreters, or environment managers;
- installed packages, libraries, or namespaces;
- shell commands or system executables;
- generic I/O, serialization, table-manipulation, plotting, or layout tools;
- project-local helpers used only for path resolution, config parsing,
  argument handling, orchestration, logging, or output writing;
- validation, audit, or status records that only permit or stop execution
  without defining scientific values or populations;
- session metadata or software versions.

Represent the scientifically meaningful operation instead. For example, a
shell command that prefilters a CSV contributes the input data and row-selection
rule, not the shell executable. A generic path helper contributes the scientific
config and resolved model code, not a helper node. A plotting package
contributes no parent object.

Project-local code is a parent only when it generates an object or defines a
scientifically meaningful transformation, model, measurement, calibration,
selection rule, or summary. A config is a parent only when its values control
such a scientific choice, not merely a path, output location, or runtime option.

Mention the scientific method in `transformation` when needed. Software
identity and version are outside the fact-card contract and can be recovered
later from code or session records if the user separately requests a
software/environment inventory.

Reserve `external` for a genuine external scientific data source or opaque
methodological boundary without a project-local generator. It never means an
installed package or command-line tool.

Directness depends on the assigned object:

- for a script or function, direct parents are the files, configs, project
  functions, or objects it directly sources, reads, calls, or consumes;
- for a durable data object, its direct parent is its immediate generating
  code or action, not the inputs of that generator;
- for a named in-memory object, its direct parent is its immediate defining
  function or action; and
- the resolver for that generator, function, or action discovers its own
  direct parents in a later router-controlled wave.

Structured provenance embedded in a data object may be used to prove its
generator. Other source paths embedded there are useful candidate hints, but
are not direct parents of the data object when an intervening generator is
represented.

## Fact-card output

Return exactly:

```markdown
# Object fact card

- object_id: `<assigned id>`
- object_type: `<allowed type>`
- canonical_path: `<repo-relative path or NA>`
- selector: `<selector or NA>`
- sha256: `<64 lowercase hex characters or NA>`
- resolution_status: `resolved|terminal|unresolved`
- generating_code_or_action: `<direct generator id/path or NA>`
- transformation: <one operational statement>
- methods_role: <one local Methods-facing statement or "none">
- known_consumers: `<ids or NA>`

## Direct dependency relationships

| dependency_id | dependency_type | relation | dependency_status | evidence | comment |
|---|---|---|---|---|---|

## Objects to enqueue

| object_id | object_type | canonical_path | selector | reason |
|---|---|---|---|---|

## Unresolved issues

- None.
```

Use only these dependency states:

```text
confirmed
candidate
unresolved
```

Every confirmed or unresolved direct parent must appear under `Objects to
enqueue` unless it is already identified as a genuine terminal external
scientific boundary. Use only the allowed object types in both dependency
tables. A contract-invalid type or dependency state, or any returned
runtime/software/tooling dependency, makes the fact card invalid; do not expect
the router to normalize it silently.

Evidence must identify a structured provenance record or a precise code
location. A file search hit alone is insufficient evidence of generation.

## Prompt template

```text
You are an Object Resolver performing a bounded Methods-provenance task.

Do not use skills. Do not read any prior canonical Methods/provenance table.
Do not edit files. Do not recurse into parent objects.

Assigned object
- object_id: {OBJECT_ID}
- object_type: {OBJECT_TYPE}
- canonical_path: {CANONICAL_PATH}
- selector: {SELECTOR}
- known consumers: {CONSUMERS}

Allowed object and dependency types
panel_endpoint, script, function, action, named_object, data_file, config,
raw_source, external, unresolved

Use those exact values. Do not emit synonyms such as code_file, data_table,
validation_record, R_package, or file.

Return scientific lineage only. Do not return languages, interpreters,
environment managers, installed packages/libraries/namespaces, shell commands,
system executables, generic I/O/serialization/table-manipulation/plotting/layout
tools, implementation-only project helpers (path/config readers, argument
handling, orchestration, logging, output writing), execution-only validation or
status records, session metadata, or software versions anywhere in the fact
card. Represent their scientifically meaningful effect as a transformation
instead. Return project code or config only when it generates an object or
defines a scientific transformation or choice. Reserve external for genuine
external scientific data or opaque methodological boundaries, never software
or command-line tools. Omit unused, display-only, administrative, runtime, and
implementation-only candidates entirely; do not inventory rejected candidates.

Permitted starting evidence
{EVIDENCE_PATHS}

First mechanically enumerate the assigned object's direct scientific data,
project-code, config, action, and source candidates. Then resolve which
candidates are actually consumed, the immediate generating code/action, the
local scientific transformation, and the local Methods role. Omit rejected
candidates. When generator identity or consumption cannot be demonstrated,
mark it unresolved.

Return exactly the Object fact-card format supplied below. Do not trace newly
discovered parents; list them for the router to enqueue. For a durable data
object, resolve its immediate generator only; leave the generator's inputs for
the generator's later resolver. For a named in-memory object, resolve its
immediate defining function or action only.

{FACT_CARD_FORMAT}
```

Use an efficient low-reasoning model for well-bounded extraction when
available. Escalate dynamic dispatch, indirect generators, or ambiguous
scientific transformations rather than lowering the evidence standard.
