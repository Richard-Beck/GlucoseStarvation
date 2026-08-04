# Target workspace policy (DRAFT)

Status: draft for review. This document describes a target state; it does not
authorize path migrations or deletion of current working material.

## Purpose

The repository should distinguish maintained project source, scientific
dependencies, and disposable workflow scaffolding. A current manuscript must
be rebuildable from maintained source and declared scientific dependencies
without relying on temporary agent or figure-development workspaces.

## Workspace classes

### Maintained project source

This class belongs in Git and includes maintained analysis code, configuration,
workflow entrypoints, documentation, tests, and project-local Codex skills and
references under `.agents/`.

Maintained source should use stable project paths or documented path resolvers.
It should not depend on `junk/` or on timestamped temporary workspaces.

### Scientific dependencies

This class contains raw observations, curated inputs, derived analytical
objects, and validated releases needed to reproduce the current scientific
outputs. Large scientific dependencies may remain outside Git, but their
location, identity, role, and availability must be documented.

The fact that an object is not consumed directly by a plotted panel is not
sufficient reason to discard it. Raw-source boundaries must be reviewed
separately from Methods-facing terminal objects. In particular, received raw
assay workbooks and their documented transformations belong to the scientific
source record even when figures consume processed tables.

### Disposable scaffolding

Agent workspaces, figure-drafting packages, temporary renders, exploratory
packages, caches, and task-specific audit machinery are disposable scaffolding.
Examples include `agent-dev/` and `manuscript_figures/`.

Scaffolding may record useful lineage, but it must not become a live dependency
of an accepted manuscript assembly. Approved outputs and required generation
code must be promoted out of scaffolding before the scaffolding is removed.

## Manuscript independence gate

A manuscript assembly is independent of scaffolding only when all of the
following are true:

1. It declares every live scientific and runtime dependency in package-local
   manifests, including identity checks where practical.
2. Its preparation, figure rebuild, manuscript render, and validation
   entrypoints all succeed in an isolated environment containing the assembly
   and only those declared dependencies.
3. Executable assembly code has no live references to `agent-dev/`,
   `manuscript_figures/`, `tmp/`, or `junk/`.
4. Historical paths retained for provenance are explicitly classified as
   origin or lineage records rather than live dependencies.
5. Dependency checks resolve symlinks and include the external targets required
   at runtime.

An isolation test that omits an assembly-preparation entrypoint does not prove
ground-up assembly independence.

## Manuscript figure packages

`manuscript_figures/` is a figure-development staging area, not a durable
release location. The current manuscript already contains assembly-local
figure-generation packages. Historical `manuscript_figures/` paths may remain
in traceability records, but the assembly's executable rebuild must use its
local generation code and declared scientific inputs.

The staging packages should be eligible for removal after the complete
manuscript independence gate passes. Until then they remain recoverable
scaffolding rather than authoritative sources.

## Agent workspaces

`agent-dev/` is temporary scaffolding and has no special archival or release
role. Accepted prose, figures, provenance, and other manuscript components must
be copied or promoted into the manuscript assembly or another agreed maintained
location.

The current manuscript still contains hard-coded `agent-dev/` references in
assembly preparation, validation, and contextual bundle code. Ongoing
scaffolding-independence work is intended to remove and test those dependencies.
`agent-dev/` must remain available until that work passes the complete gate, but
it should not be promoted into a durable project layer.

## Target data taxonomy

The target organization under `data/` is:

```text
data/
  raw/       # received or acquired source observations
  curated/   # reviewed annotations, mappings, and cleaned inputs
  derived/   # reproducible intermediate analytical objects
  releases/  # validated, immutable datasets and model outputs
```

Report exports should either be reproducibly generated outputs or part of a
named release. External microscopy and other sources that should not be copied
into this repository require an inventory recording their canonical path,
identity, scientific role, and availability expectations.

No data migration should occur until current consumers have been inventoried
and a path-transition plan has been reviewed.

## Git, local artifacts, and junk

Git should contain maintained source and compact documentation. Large local
scientific dependencies may be Git-ignored, but ignored status must not be
confused with disposability.

`junk/` is a Git-ignored, repository-local quarantine that preserves original
repository-relative paths. Material placed there is outside the active project
and must not be referenced by active code. Cleanup-specific allowlists, audits,
and reports may be retained there or recovered from Git history; they are not
ongoing project policy.

## Current transition state

- The existing cleanup keep-root allowlist and its audit machinery were
  one-time cleanup scaffolding, not a durable project convention.
- `manuscript_figures/` appears to be replaceable by the assembly-local figure
  packages, but removal should wait for a complete isolated rebuild.
- The current isolated rebuild is not yet passing. Its first failure reflects
  an external target of an `all_raw` symlink that was not mounted in the
  isolated environment, rather than evidence that the underlying TIFF was
  deleted.
- Direct `agent-dev/` references remain and assembly preparation is not yet
  covered by the isolation entrypoint set.

## Decisions required before adoption

- Confirm the final data category names and migration strategy.
- Choose the maintained home for raw-source boundary inventories.
- Decide whether manuscript releases must be fully standalone or may declare
  stable project-local scientific dependencies.
- Define the retention period for disposable scaffolding after an independence
  gate passes.
