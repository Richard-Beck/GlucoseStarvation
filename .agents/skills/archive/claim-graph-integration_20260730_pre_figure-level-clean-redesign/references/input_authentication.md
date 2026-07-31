# Claim Graph Input Authentication

Authenticate only the inputs declared to Claim Graph Integration. Do not require
or discover a particular upstream artifact type.

The stage supplies:

- its starting claim graph;
- every path constituting the canonical current evidence package;
- an internal current-evidence-item index when the package exposes stable item
  records;
- an optional prior successful `input_authentication.json` and output bundle.

The item index format is:

```json
{
  "schema_version": 1,
  "items": [
    {
      "item_id": "F1_a",
      "canonical_input": {
        "exact": "consumer-normalized payload from the supplied package",
        "declared_source_sha256": "..."
      }
    }
  ]
}
```

`canonical_input` may contain any JSON value needed to represent that supplied
evidence item. It must be derived only from the declared Claim Graph inputs. The
skill must not seek semantic interpretations, legends, source tables, or other
producer artifacts unless they are explicitly part of the supplied evidence
package.

Create and compare snapshots with:

```bash
python3 scripts/authenticate_claim_graph_inputs.py snapshot \
  --starting-graph <graph.json> \
  --evidence-input manifest=<path> \
  --item-index <current_evidence_item_index.json> \
  --output <current_input_authentication.json>

python3 scripts/authenticate_claim_graph_inputs.py compare \
  --current <current_input_authentication.json> \
  --prior <prior_input_authentication.json> \
  --output <input_reuse_report.json>
```

After writing the integrated graph, rerun `snapshot` with `--result-graph` and
`--audit-index` to write the final `input_authentication.json`. The separate
audit index has this shape:

```json
{
  "schema_version": 1,
  "audits": [
    {"item_id": "F1_a", "path": "evidence_item_audits/F1_a.md"}
  ]
}
```

Audit paths are resolved from the repository root. The snapshot records their
actual SHA-256 values. Comparison authenticates those files again and writes the
verified prior path and checksum under `reused_audits`; an absent or changed
receipt makes that item a cache miss.

Reuse rules:

- `reuse_complete`: materialize the prior complete return, rerun deterministic
  graph/package validation, and launch no evidence-item subagents.
- `reuse_partial`: reuse prior audits only for listed `reusable_item_ids` and
  launch one subagent for every current miss.
- `fresh_full`: use the ordinary one-subagent-per-current-item workflow.

Every reused audit must retain its prior path and SHA-256 in the new return. Any
missing or mismatched prior audit is a cache miss. Authentication never replaces
final graph validation or the user-fixed claim contract.
