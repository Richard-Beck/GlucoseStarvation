# Clean Claim Graph Specification

## Top-level object

```json
{
  "metadata": {},
  "claims": [],
  "evidence": [],
  "relationships": []
}
```

Required metadata:

```json
{
  "schema_version": "claim-graph/v4",
  "relation_values": ["support", "undermine"],
  "strength_values": ["strong", "moderate", "weak"],
  "authoritative_claim_contract": {
    "claims": []
  }
}
```

The graph is reconstructed from current figure audits. It has no migration
contract for claims, evidence, or edges from an older graph.

## Claims

Required fields:

- `id`: unique stable claim ID.
- `text`: claim text.
- `user_fixed`: boolean.

Optional fields may include `kind`, `tier`, `status`, and `notes`.

Every claim in `metadata.authoritative_claim_contract.claims` must appear in
`claims` with `user_fixed: true`. Every field supplied in the contract entry
must match the claim exactly. Every claim with `user_fixed: true` must appear in
the contract.

Non-authoritative claims are reconciliation products, normally umbrella claims.
They must use `user_fixed: false`.

Do not put relationships inside claim nodes.

## Evidence

Each evidence node represents one reconciled observation or a coherent group of
observations from a current figure audit.

Required fields:

- `id`: unique evidence ID.
- `figure_id`: current figure ID.
- `panels`: nonempty list of contributing panel IDs.
- `observation`: evidence-grounded observation.
- `source`: source pointer or compact source description.

Optional fields may include `kind`, `provenance`, and `notes`.

Every evidence node must be the source of at least one relationship.

## Relationships

Relationships are typed, directed records:

```json
{
  "id": "R1",
  "source_type": "evidence",
  "source_id": "E_F5_e_1",
  "target_claim_id": "C11",
  "relation": "support",
  "strength": "moderate",
  "reason": "Higher-resource regions show a low-ploidy growth advantage."
}
```

Required fields:

- `id`: unique relationship ID.
- `source_type`: `evidence` or `claim`.
- `source_id`: an evidence ID or claim ID matching `source_type`.
- `target_claim_id`: claim ID.
- `relation`: `support` or `undermine`.
- `strength`: `strong`, `moderate`, or `weak`.
- `reason`: concise justification.

Claim-to-claim self-relationships are invalid.

There are no qualification relationships. Do not use fields or values such as
`qualifies`, `qualified_by_claims`, or a `qualify` relationship.
Caveats belong in evidence, reasons, strength adjudication, or reconciliation
notes.

## Authoritative claims

The authoritative contract is the only inherited semantic contract. A prior
graph is not an input for claim migration.

If evidence undermines an authoritative claim, retain the exact claim and record
the undermining relationship. Do not silently rewrite or remove the claim.

## Reconciliation

Figure auditors map observations only to authoritative claims. The main agent
may attach evidence directly to those claims or create non-authoritative
umbrella claims and connect those claims onward.

The canonical relationship table and graph must agree on evidence-to-claim
relations and final strengths.
