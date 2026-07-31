---
name: claim-graph-integration
description: "Build or refresh a manuscript claim graph from authoritative claims and figure-level semantic descriptions. Use when Codex must run independent fresh-context audits per new or modified figure, reconcile support and undermining observations across figures, construct a clean graph without migrating inherited non-authoritative structure, and validate or render the result."
---

# Claim Graph Integration

## Purpose

Construct a clean claim graph from:

- an authoritative claim set;
- current figure bundles containing panel-level semantic descriptions;
- optional authenticated audits from a prior successful return.

Treat each figure as the semantic audit unit. Reconstruct the canonical graph
from the complete current audit set; do not migrate claims, evidence nodes, or
relationships from a prior graph.

## Core Rules

- Preserve every authoritative claim exactly as supplied.
- Use only `support` and `undermine` relationships. Do not create qualification
  relationships or qualification aliases.
- Keep fresh-context figure auditors independent of prior graphs, lineage,
  manuscript prose, skills, and repository context.
- Treat prior returns only as an audit cache for unchanged figures.
- Do not create a claim node for every panel observation.
- Leave downstream consumer updates outside this skill.

## Required Inputs

- **Authoritative claims:** stable IDs and exact claim text, plus any additional
  fields the user declares locked.
- **Current figure bundles:** one stable figure ID per bundle, the ordered panel
  IDs, and all semantic descriptions the figure auditor may read.
- **Figure index:** a canonical JSON index used to authenticate each figure
  bundle.
- **Output root.**
- **Optional prior return:** prior authentication, figure audits, and result
  graph for exact reuse checks.

The figure index and authoritative-claims index are internal workflow artifacts.
Build them only from the inputs supplied to this skill.

## Fresh-Context Figure Audits

Read `references/figure_audit_prompt.md` before launching auditors.

Launch one fresh-context subagent for every current figure lacking an
authenticated reusable audit. A new or modified figure, changed authoritative
claim set, changed audit prompt contract, missing audit, or checksum mismatch is
a cache miss.

Each auditor receives only:

- one complete figure bundle;
- the complete authoritative claim set;
- the prompt contract.

The auditor must not read skills, a prior graph, manuscript prose, code,
feedback, lineage, project documentation, or other files. It must return only a
`### Observation–claim relationships` table using:

- relation: `support` or `undermine`;
- strength: `strong`, `moderate`, or `weak`.

If the figure has no material relationship to an authoritative claim, require
an explicit no-relationship result.

Save audits under `figure_audits/` and validate them with
`scripts/validate_figure_audits.py`.

## Updating an Existing Graph

Do not patch or migrate the prior graph.

1. Authenticate the current figures, authoritative claims, prompt contract, and
   any prior audit receipts using `references/input_authentication.md`.
2. Reuse audits only for authenticated unchanged figures.
3. Audit every cache-miss figure in fresh context.
4. Reconstruct the graph from the authoritative claims and the complete current
   audit set.

A changed prior graph alone is not a reason to rerun figure auditors because
they never receive it.

## Main-Agent Reconciliation

After all current figure audits pass validation, reconcile them into the
canonical graph. The main agent owns this step and may:

- consolidate overlapping observations;
- attach evidence directly to authoritative claims;
- group related observations under useful umbrella claims;
- adjudicate relationship strength across figures;
- preserve conflicting support and undermining observations.

Do not use prior non-authoritative claims as defaults. Keep reconciliation
decisions explicit in `claim_reconciliation.md`.

Use `references/claim_graph_spec.md` when constructing or editing graph JSON.

## Authentication and Reuse

Follow `references/input_authentication.md`.

Audit reuse depends on:

- the canonical figure-bundle hash;
- the authoritative-claims hash;
- the figure-audit prompt-contract hash;
- an authenticated audit receipt.

Complete-return reuse additionally requires exact current inputs and an
authenticated prior result graph. Deterministic validation remains mandatory.

## Required Outputs

Write:

```text
claim_graph_integrated.json
claim_graph_integrated.dot
claim_graph_integrated.png
claim_graph_validation_report.txt
claim_graph_integration_report.md
claim_reconciliation.md
observation_claim_relationships.csv
figure_audits/
figure_audit_index.json
current_figure_index.json
authoritative_claims.json
input_authentication.json
input_reuse_report.json
```

`observation_claim_relationships.csv` must contain:

- `figure_id`
- `panels`
- `evidence_id`
- `observation`
- `claim_id`
- `relation`
- `worker_strength`
- `final_strength`
- `reason`
- `reconciliation_note`

## Validation

Validate figure audits:

```bash
python3 scripts/validate_figure_audits.py \
  --figure-index <output_root>/current_figure_index.json \
  --claims-index <output_root>/authoritative_claims.json \
  --audit-index <output_root>/figure_audit_index.json
```

Validate and render the graph:

```bash
python3 scripts/validate_claim_graph.py \
  <output_root>/claim_graph_integrated.json

python3 scripts/plot_claim_graph.py \
  --graph <output_root>/claim_graph_integrated.json \
  --output <output_root>/claim_graph_integrated.png \
  --dot-output <output_root>/claim_graph_integrated.dot
```

Validation must confirm:

- exact authoritative-claim preservation;
- no qualification fields or relationships;
- valid evidence and claim relationship endpoints;
- allowed relation and strength values;
- every evidence node participates in at least one relationship;
- one valid fresh or authenticated audit per current figure;
- exact input and audit authentication.

## Completion

Finish only after all current figures have valid audits, reconciliation is
complete, the clean graph and relationship table agree, authentication is
finalized, and deterministic graph and package validation pass.
