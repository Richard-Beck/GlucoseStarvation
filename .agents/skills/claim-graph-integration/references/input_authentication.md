# Figure Audit Input Authentication

Authenticate figure audits independently of any prior claim graph.

## Figure index

```json
{
  "schema_version": 1,
  "figures": [
    {
      "figure_id": "F5",
      "canonical_input": {
        "panels": [
          {"panel_id": "a", "semantic_sha256": "..."}
        ]
      }
    }
  ]
}
```

`canonical_input` may contain any JSON value needed to represent the complete
supplied figure bundle. Derive it only from declared inputs.

## Authoritative claims index

```json
{
  "schema_version": 1,
  "claims": [
    {"id": "C0", "text": "Exact authoritative claim text."}
  ]
}
```

Include every field that must be locked in the final graph.

## Audit index

```json
{
  "schema_version": 1,
  "audits": [
    {"figure_id": "F5", "path": "figure_audits/F5.md"}
  ]
}
```

## Snapshot and compare

```bash
python3 scripts/authenticate_claim_graph_inputs.py snapshot \
  --figure-index <current_figure_index.json> \
  --claims-index <authoritative_claims.json> \
  --audit-contract <figure_audit_prompt.md> \
  --evidence-input semantic_package=<path> \
  --output <current_input_authentication.json>

python3 scripts/authenticate_claim_graph_inputs.py compare \
  --current <current_input_authentication.json> \
  --prior <prior_input_authentication.json> \
  --output <input_reuse_report.json>
```

After reconciliation, rerun `snapshot` with `--audit-index` and
`--result-graph` to write final `input_authentication.json`.

## Reuse rules

An audit is reusable only when:

- its figure hash is unchanged;
- the authoritative-claims hash is unchanged;
- the audit-contract hash is unchanged;
- its recorded file exists and matches its recorded SHA-256.

Modes:

- `reuse_complete`: every current figure audit is authenticated, all declared
  inputs are unchanged, and the prior result graph is authenticated.
- `reuse_partial`: at least one but not all current figure audits are reusable.
- `fresh_full`: no current audit is reusable, including whenever authoritative
  claims or the audit contract changed.

If any figure changes, rerun only that figure audit, then reconcile a new graph
from the complete current audit set. Do not patch the prior graph.
