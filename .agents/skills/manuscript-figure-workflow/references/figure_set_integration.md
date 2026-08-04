# Figure Set Integration Workflow

## Revision Policy

When revising figures from an existing integration root:
1) Create a new root
2) Copy forward all figure generation/assembly scripts (e.g polish_figure.R or equivalent) from old root. Do this for all figures whether they are modified or not.
3) Make any required modifications as necessitated by the revision context.
4) Run all figure generation assembly scripts in the new root, irrespective of whether they are  modified or not. Do not copy final figures from old polishing roots or integration roots! Use polishing scripts to rebuild these figures from source data and inject these rebuilt figures into the new integration root.
5) Compare hashes of newly regenerated figures to prior figures. Hash matching figures do not require updated semantic interpretation or repeated visual QC. Hash fails do, unless the reason is trivial (e.g. random jittered points).

## Purpose

Use this module after manuscript figure packages have passed drafting/polishing.
The module promotes package-local figure-generation code into the integration
target, reruns that code there, verifies regenerated PNG identity against the
polishing-root final images, and then normalizes the verified outputs into a
single manuscript figure set with a canonical figure index, traceable rebuild
records, semantic-interpretation links, and package-local or temporary aliases.

This module is a packaging and rebuild-verification gate. It may also produce a
derived prior-subpanel lineage audit when a downstream task explicitly needs
prior-draft correspondence. It should not regenerate analyses, redesign figures, or
write/modify manuscript prose.

## Required Inputs

Establish these inputs before writing outputs:

- **Scope/version:** requested manuscript or redraft version, run id, and output
  root.
- **Input packages:** active polished package roots.
- **Package artifacts:** package-local contract, manifest/index, final PNGs,
  provenance, validation report, generation code/command, direct inputs,
  dependencies, notes, accepted warning rationale, and legend paths when already
  available. The package-local figure-generation script, helper files, and fixed
  configs needed to regenerate the package final PNGs must be identifiable before
  integration can complete.
- **Figure placement rules:** current package-to-manuscript mapping,
  main/supplement placement, and normalized figure ids.
- **Semantic interpretations:** panel-level interpretation files and any index
  mapping them to source package, source output, figure id, and panel id.
- **Prior manuscript draft or prior integration, when prior lineage is requested
  or downstream-required:** previous figure index/manifest, final images,
  package aliases, legend paths, package manifests, manuscript draft figure
  callouts, or other structured source that can identify earlier subpanel
  placement.

If prior lineage is requested and a prior manuscript draft/integration is not
available, proceed and record `not_available` for prior-linkage fields rather
than inventing lineage. If prior lineage is not requested, do not create a
placeholder lineage file solely to satisfy the packaging contract.

## Main-agent responsibility

The main agent owns package inventory, code promotion, package replay,
cross-package normalization, and final validation. Subagents are not required,
and their absence must not pause or block integration.

Optional subagents may assist with bounded read-only package inventory or
independent QC only when permitted by active instructions and useful for the
specific task. They must not own shared integration outputs or replace the main
agent's direct inspection and reconciliation.

For each active package, the main agent must:

- inventory final outputs, panels, legend paths when present, provenance,
  validation status, and warnings;
- identify generation scripts, commands, direct inputs, dependencies, and
  checksums;
- copy/promote the package-local figure-generation script(s), helpers, fixed
  configs, and rebuild notes into the integration target without modifying the
  polishing root;
- rerun the promoted script(s) from the integration target using the
  project-approved runner;
- write regenerated package PNGs into an integration-local package rebuild
  folder, hash them, and compare those hashes against the polishing-root
  `final_images/` PNG hashes;
- report promoted script paths, integration-local rerun commands, source PNG
  hashes, regenerated PNG hashes, hash-match status, warnings, approval
  references, and blockers;
- link package panels to semantic-interpretation files;
- compare package panels to prior draft/integration material where recoverable
  when prior lineage is requested;
- report missing files, ambiguous lineage, and blockers.

When optional inventory help is used, the main agent validates the return and
still performs cross-package normalization and final lineage decisions.

## Workflow

1. **Establish scope**
   - Identify requested version, active packages, prior
     draft/integration source, and output root.
   - Record any inference, especially when choosing between multiple manuscript
     generations.

2. **Build package inventories**
   - Inspect every input package and write a compact inventory with final
     outputs, panel rows, legend paths when present, provenance, validation
     status, semantic-interpretation links, prior-linkage hints, blockers,
     approved warnings, and approved raster-input references.
   - Preserve package inventories under the output root when practical.

3. **Promote and verify package rebuilds**
   - For each active polished package, copy package figure-generation
     entrypoints, required helpers, fixed configs, and package-local rebuild
     notes under `final_figure_scripts/package_scripts/<source_package>/`.
   - Rerun each promoted package script from the integration root or an
     integration-local package working directory.
   - Write regenerated package outputs under `package_rebuilds/<source_package>/`.
   - Hash every regenerated whole-figure PNG.
   - Compare each regenerated PNG hash to the corresponding polishing-root
     `final_images/` PNG hash.
   - Treat missing scripts, missing declared inputs, failed reruns, missing
     regenerated PNGs, and hash mismatches as blockers.
   - Record package-level rebuild results in the package inventory,
     `figure_rebuild_manifest.tsv`, and `figure_byte_identity_report.tsv`.

4. **Normalize figure identities**
   - Assign normalized `F<N>` ids for main figures and `S<N>` ids for
     supplemental figures from the requested placement rules.
   - Ensure each manuscript figure has exactly one integrated whole-figure PNG.
   - Do not create panel-specific integrated PNG filenames.
   - Do not create integrated PNG filenames with suffixes, panel letters, words,
     page markers, or other adornments.
   - If package material cannot be represented by one whole-figure PNG, assign
     the excess material a separate normalized `F<N>` or `S<N>` id.

5. **Promote verified outputs into final images**
   - Copy only verified regenerated package PNGs into `final_images/` using
     normalized manuscript figure names.
   - Do not copy unverified polishing-root PNGs directly into `final_images/`.
   - User-approved raster inputs must already have been incorporated through the
     package polishing workflow under the global raster policy in `SKILL.md`;
     integration does not add raster inputs or bypass package rebuild
     verification.
   - Preserve the source polishing-root path, integration-local regenerated path,
     normalized final path, hashes, and commands.

6. **Build the canonical figure index**
   - Follow the shared figure-object vocabulary.
   - Emit one row per manuscript-visible subpanel, or `panel_id=whole_figure`
     only for single-panel figures.
   - Include naming aliases, source roots, image paths/hashes, panel-source
     fields, semantic interpretation paths, rebuild hash status, and minimal
     prior correspondence.
   - Do not encode row roles, no-output dispositions, generic status fields, or
     approval-bypass columns in the figure index.

7. **Write the manuscript-level rebuild wrapper**
   - Write `final_figure_scripts/run_all_figures.sh`.
   - The wrapper must rerun every promoted package-generation command needed to
     regenerate the integration-local package PNGs and then copy/normalize the
     verified outputs into `final_images/`.
   - The wrapper must not perform expensive analysis reruns, model fitting,
     segmentation, classifier training, or data rebuilding unless explicitly
     requested.
   - The wrapper is the single transferable command downstream packages should
     promote or rebase to rebuild all integrated final figures from promoted
     scripts and fixed inputs. The integration root itself is a staging handoff,
     not an implied terminal dependency of a downstream manuscript package.

8. **Link semantic interpretations**
   - Require one semantic interpretation per manuscript evidence panel.
   - Check that every linked interpretation has the expected target,
     fresh-context visual description, visual claims, source clarification,
     provenance chain, caveats, and not-assessed note.
   - Do not evaluate whether the interpretation is correct or whether it
     supports a manuscript claim.
   - Mark missing or incomplete interpretations as blockers.

9. **Optionally compare to prior manuscript draft or integration**
   - Populate prior correspondence fields in the figure index when known.
   - Create a derived subpanel lineage audit only when prior-lineage output is requested,
     prior material exists, or a downstream plan explicitly requires the audit.
   - Document panels that are:
     - `unchanged`: same apparent subpanel identity and no material
       file/provenance change detected;
     - `moved_within_package`: same apparent subpanel identity, same package,
       different figure/panel placement;
     - `moved_across_packages`: same apparent subpanel identity, different
       package or figure-generation family;
     - `removed`: present in prior material but absent from current integrated
       figure set;
     - `new`: present in current integrated figure set but not found in prior
       material;
     - `changed`: same apparent role but materially different image, source,
       transform, panel composition, or provenance;
     - `uncertain`: insufficient structured evidence to choose a stronger label.
   - Use structured indexes/manifests, panel labels, short content, provenance
     paths, source package names, checksums, legend paths, and
     semantic-interpretation summaries as linkage evidence.
   - Do not perform detailed scientific interpretation in this step. The output
     is a lineage/navigation aid for downstream claim and manuscript review.

10. **Document omitted packages and validation warnings**
   - Record intentionally omitted package material in notes or validation output,
     not as figure-index rows.
   - Record warnings, blockers, and approval references in validation output.

11. **Validate packaging and rebuild identity**
   - Run `final_figure_scripts/run_all_figures.sh` from the integration root or
     into a temporary validation output directory.
   - Compare regenerated hashes to the current `final_images/` hashes.
   - Check required files, normalized names, figure-index columns, path
     existence, unique keys, semantic interpretation paths, package rebuild
     status, wrapper rebuild status, repeated figure-level field consistency,
     and prior-linkage
     completeness only when lineage was requested.
   - Emit deterministic validation results and unresolved blockers.

## Output Contract

Write outputs under a timestamped or versioned root, usually:

```text
<integration_root>/
  final_images/
    F1.png
    F2.png
    S1.png
  final_figure_scripts/
    run_all_figures.sh
    package_scripts/
      FG1_measurement_foundation/
      FG2_direct_feature_rebuild/
  package_rebuilds/
    FG1_measurement_foundation/
    FG2_direct_feature_rebuild/
  package_inventories/
  figure_set_manifest.csv
  figure_rebuild_manifest.tsv
  figure_byte_identity_report.tsv   # optional derived validation detail
  subpanel_lineage.csv              # optional derived lineage audit
  omitted_packages.md               # optional notes
  figure_set_validation_report.json
  figure_set_integration_report.md
```

`final_images/` must contain integrated whole-figure PNGs only. Name each PNG
exactly from its normalized manifest `current_figure_name`: `F<N>.png` for main figures
and `S<N>.png` for supplemental figures, where `<N>` is a positive integer. Do
not write panel-specific or suffixed integrated images such as `F2a.png`,
`S10_extra.png`, or `Figure_2.png`.

### figure_set_manifest.csv

This is the canonical figure index. Follow the shared figure-object vocabulary.
It replaces separate numbering-crosswalk, semantic-index, image-artifact, and
panel-source contract files.

### figure_rebuild_manifest.tsv

Follow the shared rebuild-record vocabulary. Join to the figure index by
natural figure/stage keys, not opaque ids.

### figure_byte_identity_report.tsv

Optional derived validation detail. Do not make it a primary contract object.

### semantic_interpretation_index.csv

Do not require this as a handoff object. Semantic coverage is canonical in
`figure_set_manifest.csv` via `semantic_interpretation_path`. A compact semantic
coverage table may be generated as a derived review view.

### subpanel_lineage.csv

Optional prior-lineage audit. Produce this only when a prompt, plan, prior-draft
comparison, claim-graph update, Results/legend revision, or manuscript review
explicitly needs prior subpanel correspondence. Within
`manuscript-figure-workflow`, this file is produced by the figure-set integration
module and is not consumed by polishing, semantic interpretation, drafting, or
ideation.

When produced, derive it from manifest prior-correspondence fields and available
prior-integration evidence. It is not consumed by polishing, semantic
interpretation, drafting, or ideation.

## Deterministic Validation Targets

A validation pass should be able to check:

- Required output files exist.
- Every active package has promoted generation scripts under
  `final_figure_scripts/package_scripts/`.
- Every active package script reruns successfully from the integration target.
- Every regenerated package PNG has a SHA256 recorded in the figure index or
  rebuild manifest.
- Every regenerated package PNG hash matches its polishing-root final PNG hash.
- `final_figure_scripts/run_all_figures.sh` regenerates all normalized
  `final_images/` PNGs byte-identically.
- `final_images/` contains exactly one PNG for each manifest `current_figure_name` that
  has an image-bearing row.
- Every image-bearing manifest row points inside `final_images/`.
- Integrated PNG basenames match normalized `current_figure_name` values.
- Every image-bearing `current_figure_name` and integrated PNG basename matches exactly
  `F<N>` or `S<N>`, where `<N>` is a positive integer.
- No integrated PNG is panel-specific or suffixed.
- Manifest required columns exist, keys are unique, repeated figure-level fields
  are consistent, and referenced paths exist.
- Active package outputs do not point to rejected drafts, superseded draft roots,
  or obsolete final-PNG folders.
- Every active input package has a package inventory.
- Omitted package material is documented outside the figure index when relevant.
- Every manuscript-visible panel has a semantic interpretation path in the
  figure index.
- Every required semantic-interpretation section is present.
- If prior lineage was requested, `subpanel_lineage.csv` exists and records
  either prior-linkage status for current/prior panels or a clear
  `not_available` status.
- `figure_set_integration_report.md` records selected inputs, optional
  delegation used or `none`,
  package inventory paths, promoted script paths, package rebuild summary,
  wrapper rebuild summary, validation summary, blockers, warnings,
  and unresolved decisions.

## Completion Standard

Finish only after active polished figure packages have been rebuilt from
integration-local promoted scripts, regenerated package PNG hashes match
polishing-root final PNG hashes, normalized final images are generated by
`final_figure_scripts/run_all_figures.sh`,
semantic interpretation paths are recorded in the figure index, requested prior
lineage is recorded or explicitly unavailable, and blockers are resolved or
clearly reported.

The final response should list the figure-set integration root, active packages,
semantic interpretation coverage status, prior-linkage status when requested,
validation status, and blockers or approved warnings.
