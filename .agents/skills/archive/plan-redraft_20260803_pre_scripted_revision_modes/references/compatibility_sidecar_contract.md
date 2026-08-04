# Compatibility Sidecar Contract

Use this contract for the planner-owned pre-Assembly consistency gate. The
planner owns the common format and comparison; each stage consumer owns which
cross-stage manuscript facts are relevant within its scope and how those facts
are authenticated.

Each returned stage bundle must contain `compatibility_sidecar.json`:

```json
{
  "schema_version": 1,
  "consumer": "results-text",
  "bundle_id": "stable return identifier",
  "completion_attestation": true,
  "facts": [
    {
      "fact_id": "S13a.qc_pass_count.MDA-MB-231",
      "value": 5,
      "unit": "selected model fits",
      "scope": "posterior sampling QC; denominator 5",
      "status": "validated",
      "source_artifacts": [
        {"path": "path/to/source.csv", "sha256": "..."}
      ],
      "manuscript_artifacts_using_fact": ["path/to/Figure_S13.md"],
      "notes": ""
    }
  ]
}
```

Rules:

- Set `completion_attestation: true` only after considering all exact counts,
  units, definitions, panel identities, category meanings, terminology, and
  interpretive boundaries in the stage output that another stage may repeat,
  define, interpret, or depend upon.
- Use stable upstream identifiers in `fact_id` when available: figure/panel IDs,
  claim IDs, method IDs, or another canonical subject identifier plus predicate.
- Record only cross-stage manuscript facts, not every statement inspected.
- Use `status: validated` when the consumer authenticated the value from the
  listed source artifacts. Use `status: unresolved` when it could not.
- Report a fact even when another stage may disagree. The consumer does not
  reconcile or suppress another owner's assertion.
- Every source artifact must exist. Supply its SHA-256 whenever it is a regular
  file.
- An empty `facts` list is valid when the consumer has considered its complete
  return and has no cross-stage facts to report.

Before a first Manuscript Assembly handoff, the planner writes a registry:

```json
{
  "schema_version": 1,
  "sidecars": [
    {"consumer": "analysis", "path": "path/to/compatibility_sidecar.json"}
  ],
  "accepted_exceptions": [
    {
      "fact_id": "optional.fact.id",
      "reason": "Exact user-authorized or run-policy reason",
      "authorized_by": "user"
    }
  ]
}
```

The registry must include exactly one current sidecar for each of the seven
pre-Assembly consumers: `analysis`, `manuscript-figure-workflow`,
`claim-graph-integration`, `results-text`, `method-table-provenance`,
`manuscript-legend-writing`, and
`serve-manuscript-abstract-introduction-discussion`.

Run:

```bash
python3 .agents/skills/plan-redraft/scripts/check_compatibility_sidecars.py \
  --registry <registry.json> \
  --report-json <report.json> \
  --report-md <report.md>
```

The checker compares `value`, `unit`, and `scope` for matching `fact_id` values.
Missing or malformed sidecars, unresolved unaccepted facts, and unaccepted
cross-consumer contradictions block Assembly. Accepted exceptions remain
visible warnings.
