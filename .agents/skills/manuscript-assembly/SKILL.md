---
name: manuscript-assembly
description: "Assemble and validate a coherent manuscript release or review package from finalized components. 
Use when Codex needs to create or update the final folder for a manuscript, assemble and validate sources."
---

# Manuscript Assembly

## Purpose

Use this skill to produce a durable manuscript package according to defined guidelines,
or clearly block manuscript assembly if the required inputs are incomplete.

This skill owns package assembly, consistency checks, final-folder structure,
and assembly status. If assembly reveals missing, stale, or misorganized content, stop and clearly report the issue.

## Assembly Root

When the user supplies an assembly root, use it. Otherwise create a new
versioned `<assembly_root>` following the local project convention. Infer its
identifier from the source manuscript state and current date. Do not overwrite
an existing assembly root.

## Terminal Dependency Closure

A completed assembly must rebuild and validate using only assembly-local
contents, explicitly declared scientific dependencies, and explicitly declared
runtime dependencies. Inputs used only to construct the assembly may be
retained as lineage records, but terminal entrypoints must not require them.

Keep these roles distinct:

- **assembly-local content**: copied or promoted manuscript sources, assets,
  rebuild code, receipts, and validation contracts;
- **live dependencies**: external scientific data or maintained code required
  by a terminal entrypoint;
- **runtime dependencies**: interpreters, containers, libraries, and system
  facilities required only to execute terminal entrypoints;
- **lineage sources**: origins used for navigation or one-time import and never
  exposed during terminal validation.

Do not infer roles from path names. Classify each dependency from its use. A
recorded source or lineage path is not implicitly a live dependency.

## Package Contract

Organize the final folder around durable manuscript concepts, not workflow
history. Use this structure unless the local project already has an equivalent
layout:

```text
<assembly_root>/
  README.md
  status.json
  draft/
  source/
  assets/
  evidence/
  traceability/
  review_state/
  validation/
  rebuild/
```

Required contents:

- `README.md`: human-facing index, package status, source state, and first files
  to inspect.
- `status.json`: machine-readable status with `pass`, `warn`, or `block`,
  assembly timestamp, lineage sources, and validation report paths.
- `draft/`: rendered manuscript draft and any rendered supplement or review
  companion.
- `source/`: editable manuscript text, captions, references, and sidecar source
  files used to render the draft.
- `assets/`: final figures, tables, supplemental assets, and an inventory of
  where each appears in the manuscript.
- `evidence/`: concise claim/conclusion support record and any current
  evidence-state files needed to audit claim strength.
- `traceability/`: separate lineage, live-dependency, runtime-dependency, and
  copied-input identity records with paths, versions, checksums, and roles.
- `review_state/`: what feedback or requested changes the package responds to,
  what remains open, and what was skipped, deferred, or accepted as an exception.
- `validation/`: human-readable and machine-readable final assembly validation.
- `rebuild/`: renderer scripts, figure-asset rebuild scripts, configs, and
  command notes sufficient to regenerate the rendered draft and assembled final
  manuscript assets from package inputs.

One-time import or preparation code is lineage material, not a terminal rebuild
entrypoint. Keep it outside the terminal rebuild surface when it is retained.

### Figure Asset Rebuild Contract

The assembly package must include an assembly-local figure rebuild package:

```text
<assembly_root>/
  assets/
    figures/
      F1.png
      F2.png
      S1.png
  rebuild/
    figures/
      README.md
      run_all_figures.sh
      figure_rebuild_manifest.tsv
      package_scripts/
	F1.R
	F2_S1.R
  validation/
    figure_rebuild_validation.tsv
```

`rebuild/figures/run_all_figures.sh` must regenerate every PNG in `assets/figures/` from scripts present in `assets/figures/package_scripts/`. These scripts must generate and 
assemble figures from source data (e.g. .csv/.tsv/.Rds) except in narrow cases where use of raster images is unavoidable (e.g. when image components are schematics, cartoons, or biological images).
Generally all plots (barplots, scatterplots, pie charts, violin plots, etc.) that can be generated from analysis-derived datasets must be. Aim for one generating script per figure or closely connected 
group of figures. Assembly scripts should avoid performing significant reanalysis, permitted operations usually include sub-setting, merging & joining, computation of simple summary statistics etc.


`figure_rebuild_manifest.tsv` must have one row per assembled final figure with:
`figure_id`, `asset_path`, `rebuild_output_path`, `source_package`,
`polish_root`, `polishing_script`, `rebuild_command`, `direct_inputs`,
`dependency_paths`, `approved_raster`, `expected_sha256`, and
`accepted_exception`.

Treat source-package and polishing-root fields as lineage. The promoted script,
assembly-local assets, direct inputs, and dependency paths define the rebuild
closure. Resolve symlinked inputs far enough to declare the paths that must
actually be exposed during replay.

## Assembly Tools

Use the bundled project-independent scripts with an assembly-local JSON config:

- `scripts/assemble.py` materializes a reviewable candidate live-dependency
  manifest from configured manifests, code roots, and explicit additions.
- `scripts/validate_scaffolding_independence.py` independently replays terminal
  entrypoints with only the locked live and runtime manifests exposed.
- `assets/assembly_config.template.json` is the project-independent config
  template; copy and specialize it inside the assembly.

The candidate manifest is discovery output, not an approved allowlist. Review,
deduplicate, classify, and lock it before validation. The validator must never
add dependencies or revise the locked manifests in response to a failure.
Project-specific paths, runtime details, entrypoints, and lineage manifests
belong in the local config, not in this skill or the bundled scripts.

## Workflow

1. **Identify sources**
   Determine the manuscript state to assemble and upstream content roots.

2. **Classify components**
   Mark each required component as finalized, stale, missing, out of scope. BLOCK
   and close out immediately if all required components are not present and finalized.

3. **Create or update the assembly root**
   Build the package layout. Copy small, durable manuscript-facing artifacts.
   Declare required large scientific inputs as live dependencies. Record
   construction-only origins and workflow logs as lineage without making them
   terminal dependencies.

4. **Materialize and lock dependencies**
   Configure the assembly tool, materialize candidate live dependencies,
   reconcile them against the actual terminal entrypoints, and write separate
   locked live and runtime manifests. Record construction origins separately.

5. **Render or collect the draft**
   Use the project renderer when available. If rendering requires substantive
   text, figure, or legend changes, stop and route to the owner. If rendering is
   mechanical, run it and place the result under `draft/`.

6. **Build package-level indexes**
   Write the asset inventory, evidence support summary, upstream input register,
   review-state summary, rebuild notes, figure rebuild manifest, `README.md`,
   and `status.json`.

   The asset inventory must link each figure asset to its corresponding
   `figure_rebuild_manifest.tsv` row and record whether the figure is rebuilt
   from polishing scripts or carried as an approved immutable raster exception.

7. **Validate**

   Validate figure rebuildability by running
   `rebuild/figures/run_all_figures.sh` into a temporary validation output
   directory, comparing every regenerated PNG against `assets/figures/`, and
   writing `validation/figure_rebuild_validation.tsv`. Treat checksum mismatch,
   missing rebuild commands, missing scripts, undeclared inputs as `BLOCK` unless the user explicitly accepts an
   exception.

   Then run the standalone dependency-closure validator. It must copy the
   assembly into an isolated project view, expose only locked live and runtime
   dependencies, leave lineage sources unavailable, and run every terminal
   rebuild and validation entrypoint. Any required undeclared input is `BLOCK`.

8. **Close out**
   Report the assembly root, rendered draft path, status, blockers, warnings,
   accepted exceptions, and whether a project navigation document should be
   updated.

## Validation Checklist

At minimum, validate:

- Figure rebuildability.
- Internal consistency - are figures named consistently across all inputs?
- Captions or legends exist for every rendered figure/table asset that needs
  one.
- Manuscript cross-references, anchors, figure/table callouts, and bibliography
  links resolve when tooling permits checking them.
- The evidence support record identifies all figures produced by package_scripts/ and all datasets consumed by package_scripts/, and fully traces the computational steps required to generate these datasets.
- Rendered draft files are regenerated.
- Rendered HTML, if produced, embeds or links assets according to the requested
  delivery format and has unique anchors.
- Rebuild instructions are sufficient for another agent to regenerate the draft.
- Terminal rebuild and validation pass with lineage sources unavailable.
- Terminal validators authenticate copied content from assembly-local receipts
  rather than reopening source copies.
- Journal-facing draft and captions do not expose internal paths, commands, or
  workflow labels except in an audit/provenance section.

Use `PASS`, `WARN`, and `BLOCK` consistently:

- `PASS`: assembly is coherent and ready for the next review/submission step.
- `WARN`: assembly is usable, but carries explicit exceptions or
  nonblocking risks. Contract conformance for every consumed required input must
  still pass.
- `BLOCK`: assembly should not be treated as final because required source,
  asset, evidence, review, render, or validation state is missing or
  contradictory.

## Boundaries

Do not:

- rerun scientific analyses;
- perform new figure design or sfigure revision;
- use assembly to repair missing provenance or invent rebuild commands;
- overwrite upstream outputs while validating rebuild commands;
- treat a lineage record as permission to read that source during terminal
  rebuild or validation;
- allow the isolation validator to discover, add, or repair dependencies;
- revise scientific prose beyond mechanical packaging fixes;
- write or edit figure legends or captions;
- edit claim/evidence mappings;
- reinterpret raw feedback as a package owner;
- bury missing required content as a warning without user approval;
- assemble from mismatched, ambiguous, or failed inputs
  by applying local fixes, inferred mappings, or manual substitutions.

Do:

- preserve finalized inputs;
- run mechanical figure rebuild validation in an assembly-local or temporary
  output directory and compare checksums;
- require clean upstream outputs before treating a required component as
  assembly-ready;
- keep workflow logs out of the main manuscript package unless they are needed
  as linked audit material;
- make the package inspectable and regenerable;
- preserve a project-local assembly config and locked dependency manifests;
- keep assembly, reconciliation, and validation under the main agent unless the
  user explicitly requests bounded delegation.

## Completion Standard

Finish with:

A populated assembly root following the package contract OR a BLOCK report detailing reasons you cannot proceed.
