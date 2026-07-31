---
name: manuscript-assembly
description: "Assemble and validate a coherent manuscript release or review package from finalized components. 
Use when Codex needs to create or update the final folder for a manuscript, assemble and validate sources."
---

# Manuscript Assembly

## Purpose

Use this skill to produce a durable manuscript package according to define guidelines,
or clearly block manuscript assembly if the required inputs are incomplete.

This skill owns package assembly, consistency checks, final-folder structure,
and assembly status. If assembly reveals missing, stale, or misorganized content, stop and clearly report the issue.

## Assembly Root

When the user supplies a target folder, use it. Otherwise create a new root
under:

```text
manuscripts/<run_id>/
```

Infer `run_id` from the source manuscript state and current date. Do not
overwrite an existing assembly root. Reference large upstream artifacts by path
and checksum unless the final package needs a local copy for portability.

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
  assembly timestamp, source roots, and validation report paths.
- `draft/`: rendered manuscript draft and any rendered supplement or review
  companion.
- `source/`: editable manuscript text, captions, references, and sidecar source
  files used to render the draft.
- `assets/`: final figures, tables, supplemental assets, and an inventory of
  where each appears in the manuscript.
- `evidence/`: concise claim/conclusion support record and any current
  evidence-state files needed to audit claim strength.
- `traceability/`: upstream input register with paths, versions, checksums, and
  rationale for copied versus referenced artifacts.
- `review_state/`: what feedback or requested changes the package responds to,
  what remains open, and what was skipped, deferred, or accepted as an exception.
- `validation/`: human-readable and machine-readable final assembly validation.
- `rebuild/`: renderer scripts, figure-asset rebuild scripts, configs, and
  command notes sufficient to regenerate the rendered draft and assembled final
  manuscript assets from package inputs.

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

## Workflow

1. **Identify sources**
   Determine the manuscript state to assemble & upstream content roots.

2. **Classify components**
   Mark each required component as finalized, stale, missing, out of scope. BLOCK
   and close out immediately if all required components are not present and finalized.

3. **Create or update the assembly root**
   Build the package layout. Copy small, durable manuscript-facing artifacts.
   Reference large analysis artifacts, raw data, and workflow logs by path and
   checksum unless portability requires copies.

4. **Render or collect the draft**
   Use the project renderer when available. If rendering requires substantive
   text, figure, or legend changes, stop and route to the owner. If rendering is
   mechanical, run it and place the result under `draft/`.

5. **Build package-level indexes**
   Write the asset inventory, evidence support summary, upstream input register,
   review-state summary, rebuild notes, figure rebuild manifest, `README.md`,
   and `status.json`.

   The asset inventory must link each figure asset to its corresponding
   `figure_rebuild_manifest.tsv` row and record whether the figure is rebuilt
   from polishing scripts or carried as an approved immutable raster exception.

6. **Validate**

   Validate figure rebuildability by running
   `rebuild/figures/run_all_figures.sh` into a temporary validation output
   directory, comparing every regenerated PNG against `assets/figures/`, and
   writing `validation/figure_rebuild_validation.tsv`. Treat checksum mismatch,
   missing rebuild commands, missing scripts, undeclared inputs as `BLOCK` unless the user explicitly accepts an
   exception.

7. **Close out**
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
- keep assembly, reconciliation, and validation under the main agent unless the
  user explicitly requests bounded delegation.

## Completion Standard

Finish with:

A populated assembly root following the package contract OR a BLOCK report detailing reasons you cannot proceed.
