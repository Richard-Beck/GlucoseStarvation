# Claim Graph Prototype Specification

This folder holds a lightweight manuscript claim/evidence/dependency graph. The canonical object is
`manuscript_claim_evidence_graph_v3_prototype.json`.

## Top-Level Object

Required top-level keys:

- `metadata`: graph-level identifiers, allowed values, conventions, and user-fixed claim contract.
- `claims`: claim nodes.
- `evidence`: evidence nodes keyed to manuscript figure panels or supplemental panels.
- `methods`: method/procedure nodes.
- `assumptions`: assumption/risk nodes.

All node IDs must be unique within their section. Claim references must point to existing claim IDs, and evidence,
method, and assumption references must point to existing nodes in their own sections.

## Claim Nodes

Required claim fields:

- `id`: stable claim identifier.
- `text`: human-readable claim text.
- `tier`: one of `metadata.claim_tiers`.
- `type`: free-text class such as `biological_claim`, `prediction_claim`, or `model_result`.
- `status`: one of `metadata.status_values`.
- `supported_by_claims`: claim IDs with support edges into this claim.
- `qualified_by_claims`: claim IDs with qualification/caveat edges into this claim.
- `supports`: claim IDs this claim supports.
- `qualifies`: claim IDs this claim qualifies.
- `undermines`: claim IDs this claim undermines.
- `evidence`: evidence IDs.
- `depends_on_methods`: method IDs.
- `depends_on_assumptions`: assumption IDs.
- `notes`: short implementation or interpretation note.

Optional claim fields:

- `user_fixed`: boolean. When `true`, the claim is locked by explicit user instruction.

Current tiers:

- `tier_0`: top-level manuscript claims. Multiple tier-0 claims may coexist; no single root claim is required.
- `tier_1`: primary supporting claims supplied by the user.
- `supporting_result`: data/model/simulation results that support or interpret tier-0/tier-1 claims.
- `boundary_or_undermining`: scope limits, caveats, exceptions, or constraints.

## User-Fixed Claims

Claims with `user_fixed: true` are governed by `metadata.user_fixed_claim_contract.claims`.

Rules:

- A user-fixed claim must not be removed from `claims`.
- A user-fixed claim's `text`, `tier`, and `status` must not be changed unless the user explicitly requests that change.
- Every claim listed in the contract must have `user_fixed: true`.
- Every claim with `user_fixed: true` must be listed in the contract.
- Edges, evidence links, method links, assumption links, type, and notes may be edited if they improve graph accuracy and
  do not change the locked claim text/tier/status.

The validator enforces presence of user-fixed claims and equality of the locked fields against the metadata contract.

## Edge Semantics

- `supports` and `supported_by_claims` encode positive support.
- `qualifies` and `qualified_by_claims` encode caveats, scope limits, or tension that modifies interpretation without
  directly refuting the claim.
- `undermines` encodes direct conflict or refutation.

For consistency, claim-to-claim relationships should normally be represented on both sides when practical, but the
plotter and validator tolerate either directional source.

## Evidence, Methods, and Assumptions

Evidence nodes should include:

- `id`
- `kind`
- `source`
- `notes`

Method nodes should include:

- `id`
- `kind`
- `source`
- `summary`

Assumption nodes should include:

- `id`
- `text`
- `risk`

## Maintenance

Run validation after edits:

```bash
python3 .agents/skills/claim-graph-integration/scripts/validate_claim_graph.py \
  <output_root>/claim_graph_integrated.json
```

Regenerate the visual graph after claim or edge edits:

```bash
python3 .agents/skills/claim-graph-integration/scripts/plot_claim_graph.py \
  --graph <output_root>/claim_graph_integrated.json \
  --output <output_root>/claim_graph_integrated.png \
  --dot-output <output_root>/claim_graph_integrated.dot
```
