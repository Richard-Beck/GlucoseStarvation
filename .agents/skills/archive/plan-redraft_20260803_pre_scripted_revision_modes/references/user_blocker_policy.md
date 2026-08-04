# User-Decision Blocker Policy

Use this planner-owned contract for user decisions that permit provisional
workflow continuation. It does not convert missing mandatory inputs, unsafe or
destructive actions requiring new authority, failed deterministic validation,
or inability to construct a required return bundle into warnings.

## Run Policy

The runner writes `<redraft_root>/run_policy.json`:

```json
{
  "schema_version": 1,
  "generated_at": "ISO-8601 timestamp",
  "user_blocker_mode": "interactive",
  "source": "runner"
}
```

`interactive` is the default. It pauses once for unresolved user decisions.
`defer_to_assembly` authorizes workflow continuation, not a scientific
conclusion; unresolved decisions remain Assembly warnings.

## Consumer Sidecar

Every stage return writes `user_blockers.json`:

```json
{
  "schema_version": 1,
  "consumer": "claim-graph-integration",
  "bundle_id": "stable return identifier",
  "completion_attestation": true,
  "blockers": [
    {
      "blocker_key": "S13a.choose-support-status",
      "question": "Which support status should be assigned?",
      "blocking_scope": "user-decision",
      "relevant_artifacts": [
        {"path": "path/to/evidence.json", "sha256": "..."}
      ],
      "provisional_handling": "Retain the current status pending review",
      "options": []
    }
  ]
}
```

The consumer owns blocker meaning, the stable `blocker_key`, question, and
provisional handling. The planner does not invent or reinterpret them. Use only
`blocking_scope: user-decision`. Report operational or contract failures through
the ordinary stage status instead. An empty blocker list is valid.

Every relevant artifact must be a regular file with its current SHA-256. The
planner utility fingerprints the contract version, consumer, blocker key, and
sorted artifact hashes. Question rewording does not create a new blocker;
changed relevant content does.

After writing the sidecar, obtain its deterministic fingerprints without
changing planner state:

```bash
python3 .agents/skills/plan-redraft/scripts/reconcile_user_blockers.py fingerprint \
  --sidecar <return>/user_blockers.json \
  --output <return>/user_blocker_fingerprints.json
```

`user_blocker_fingerprints.json` is a required return artifact whenever the
sidecar contains a blocker.

In `defer_to_assembly` mode, cite the matching full fingerprint in any
feedback-manager `deferred_by_run_policy` signoff. In interactive mode, return
the blocker and leave affected feedback `needs_followup`; after a user
resolution, the planner routes that consumer back for ordinary signoff. If an
awaiting decision is later deferred, route the consumer back for the bounded
policy signoff before advancing.

## Reconciliation

After accepting a returned bundle, run:

```bash
python3 .agents/skills/plan-redraft/scripts/reconcile_user_blockers.py reconcile \
  --sidecar <return>/user_blockers.json \
  --ledger <redraft_root>/user_blockers/ledger.json \
  --run-policy <redraft_root>/run_policy.json \
  --user-action-file <redraft_root>/user_action_required.json \
  --report-json <redraft_root>/user_blockers/<turn>_report.json \
  --report-md <redraft_root>/user_blockers/<turn>_report.md
```

The ledger dispositions are `awaiting_user`, `deferred_to_assembly`, and
`resolved`. Repeated fingerprints retain their disposition and increment their
occurrence count. Switching an awaiting blocker to a `defer_to_assembly` run
converts it to deferred; prior deferred or resolved dispositions remain durable.

In interactive mode, `user_action_required.json` has status `awaiting_user` and
consolidates every awaiting fingerprint. Write no new stage handoff. The runner
exits cleanly after observing it. In defer mode, continue the workflow.

Record later user dispositions with:

```bash
python3 .agents/skills/plan-redraft/scripts/reconcile_user_blockers.py resolve \
  --ledger <redraft_root>/user_blockers/ledger.json \
  --resolution-file <resolution.json> \
  --user-action-file <redraft_root>/user_action_required.json
```

The resolution file contains `schema_version: 1` and a `resolutions` list. Each
entry requires a fingerprint, `authorized_by: user`, a nonempty resolution, and
disposition `resolved` or `deferred_to_assembly`.

## Assembly

Route the ledger to Manuscript Assembly. Require
`review_state/deferred_user_decisions.md` listing every deferred fingerprint,
consumer, question, provisional handling, and authenticated artifact. Deferred
records produce package `WARN`, not `BLOCKED`. Assembly must not describe them
as scientifically resolved.
