---
name: serve-manuscript-abstract-introduction-discussion
description: "Create and carry forward ordered manuscript A/I/D bundles. Each bundle contains a required context document and may contain a literature map and externally authored Abstract, Introduction, and Discussion source. Use when manuscript context must be packaged for external authorship or the newest available external A/I/D text must be served verbatim for manuscript injection."
---

# Serve Manuscript Abstract, Introduction, and Discussion

## Bundle contract

A bundle contains:

```text
<bundle>/
  context.txt                 # required
  manifest.json               # required
  literature_map.json         # optional
  aid_source.md                # optional
```

`manifest.json` records a stable bundle ID and the path, role, byte size, and
SHA-256 of every bundle payload file. A literature map and A/I/D source belong
only to the context in the same bundle.

## Input contract

Require:

- `output_root`.
- `context_sources`: a nonempty ordered list of `{path, role}` records from
  which to build the new current bundle.

Optionally accept:

- `previous_output_root`: a prior output satisfying the output contract.
- `context_separator`, defaulting to one blank line.

Do not accept loose literature maps or A/I/D sources outside a bundle.

## Output contract

Return:

```text
<output_root>/
  current_bundle/
  old_bundle_1/
  old_bundle_2/
  ...
  served/
    abstract.md
    introduction.md
    discussion.md
    references.md              # when supplied
    manifest.json
  status.json
```

Concatenate `context_sources` in declared order using the recorded separator and
create `current_bundle/context.txt` and its manifest. If a prior output is
supplied, validate its complete output contract. When valid, copy its
`current_bundle` to `old_bundle_1` and each prior `old_bundle_N` to
`old_bundle_(N+1)`, preserving bundle contents and IDs. When invalid, ignore the
entire prior root, copy none of its bundles, and continue. Record its path,
`ignored_invalid` disposition, and exact validation failures in `status.json`.

To serve text, inspect bundles in order: `current_bundle`, `old_bundle_1`,
`old_bundle_2`, and so on. Use the first `aid_source.md` found. Mechanically
extract its Abstract, Introduction, Discussion, and References blocks without
changing body text. If no bundle contains `aid_source.md`, create zero-byte
Abstract, Introduction, and Discussion files and omit References. Record the
source bundle ID and served-file hashes in `served/manifest.json`.

## Boundaries

- Do not search literature or create, merge, repair, or assess literature maps.
- Do not interpret, review, or revise A/I/D text or alter citations or formatting.
- Do not choose bundles by timestamps or filenames outside the declared order.
- Do not render or assemble the manuscript.

## Validation

Require every bundle to contain a valid context and manifest, verify all hashes,
and reject orphan maps or A/I/D sources. Verify that bundle ordering is complete
and unique and that served files are either byte-preserving extracts from the
recorded bundle or the explicit zero-byte empty state.
