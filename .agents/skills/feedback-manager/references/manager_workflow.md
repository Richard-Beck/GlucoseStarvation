# Feedback Manager Workflow

## Core Rule

The feedback-manager tracks state. It does not interpret feedback meaning.

Do not infer targets, summarize user intent, classify directive strength, create
feedback items, assign semantic owners, or decide what a consumer should do.
Store raw feedback pointers and let consumer skills read the raw feedback and
relevant manuscript context. Only consumers create or update feedback items.

## Project Defaults

- Prefer storing feedback-manager outputs under the root whose manuscript state is
  being reviewed or revised.
- For integrated manuscript review feedback, prefer:

```text
<integration_root>/feedback_archive/<feedback_id>/
<integration_root>/feedback_tracking/<feedback_id>/
```

- For active redraft execution feedback tracking, prefer:

```text
<redraft_root>/feedback_tracking/<feedback_id>/
```

- Use date-based feedback IDs such as `20260701_manuscript_draft_review` unless
  the user provides a specific identifier.
- Treat archives as immutable once checksums are written. Put mutable tracking
  files in `feedback_tracking/`, not inside the archive packet.

## Archive Intake

When asked to archive new feedback, create a thin immutable packet:

```text
feedback_archive/<feedback_id>/
  README.md
  raw_feedback/
  reviewed_artifact/
  metadata/archive_manifest.json
  metadata/git_context.txt
  SHA256SUMS
```

The packet should preserve direct/raw user feedback files without interpretation,
the exact reviewed artifact when available, source paths, mtimes, sizes,
checksums, current git HEAD, dirty worktree status, whether the reviewed artifact
is self-contained, and any accepted reason a reviewed artifact is unavailable.

Do not copy large source trees, separate figure PNGs, legends, or section sidecars
when the reviewed artifact is a self-contained HTML containing the reviewed
content. Record the rationale in the manifest.

## Tracking Outputs

Canonical tracking state always contains spans. Items are optional and appear
only after consumers create trusted work units:

```text
feedback_tracking/<feedback_id>/
  feedback_spans.jsonl
  feedback_items.jsonl  # optional until a consumer creates an item
```

No separate source list is canonical. Raw feedback source paths are already
recorded on each span and can be derived from `feedback_spans.jsonl`.

Generated `served_feedback/` directories may exist under a tracking root when
`serve_feedback` is run without `--output_dir`; they are consumer packet outputs,
not tracking ledgers.

### feedback_spans.jsonl

One JSON object per raw-feedback span. Required fields:

```json
{
  "span_id": "SPAN-20260701-014",
  "raw_source": "agent-dev/.../raw_feedback/user-feedback.txt",
  "start_line": 42,
  "end_line": 68,
  "source_sha256": "optional-file-or-span-checksum",
  "excerpt_hash": "optional-short-hash",
  "created_by": "feedback-manager",
  "created_at": "ISO-8601 timestamp",
  "consumer_responses": []
}
```

Span records are mechanical pointers. Do not add target tags, target confidence,
locator notes, summaries, routing labels, or semantic descriptions.

Consumers sign off spans after reading the raw feedback and deciding what, if
anything, to do within their scope. A span signed off by a consumer is not
reopened; it remains hidden from that consumer on later serve calls. Use
`status: "needs_followup"` only when a consumer has inspected a span but it
should still be served again to that same consumer.

### feedback_items.jsonl

This file is absent or empty until a consumer creates an item. When present, it
contains one JSON object per durable consumer-created feedback work unit. A
feedback item groups one or more spans and stores consumer responses for that
item.

The feedback-manager does not create items. A consumer creates an item when it
determines that one or more spans require action, signoff, or a documented
no-change/not-applicable response within that consumer's scope. A consumer may
also update an existing item when it determines the item already represents the
relevant work unit.

Required fields:

```json
{
  "item_id": "ITEM-20260706T123456-001-results-text",
  "source_spans": ["SPAN-20260701-014"],
  "created_by_consumer": "results-text",
  "created_at": "ISO-8601 timestamp",
  "active": true,
  "consumer_responses": []
}
```

Optional fields:

```json
{
  "updated_by_consumer": "manuscript-integration",
  "updated_at": "ISO-8601 timestamp",
  "item_update_note": "Consumer-owned update note.",
  "supersedes": ["ITEM-..."],
  "related_items": ["ITEM-..."]
}
```

Do not store target fields, item summaries, directive classifications, global
pending state, downstream queues, or final closure fields on feedback items.

Consumer responses are appended to `consumer_responses` by the signoff helper:

```json
{
  "response_id": "RESPONSE-20260706T123456-001-results-text",
  "consumer_name": "results-text",
  "responded_at": "ISO-8601 timestamp",
  "status": "signed_off",
  "response_type": "scope_satisfied|not_applicable|no_change_decision|validation_only|artifact_changed|partial_change|deferred_by_run_policy",
  "response": "Consumer-owned response text.",
  "checked_artifacts": [
    {
      "path": "agent-dev/.../manuscript_draft.html",
      "sha256": "artifact checksum when practical",
      "role": "checked artifact"
    }
  ],
  "evidence_artifacts": [],
  "changed_artifacts": [],
  "residual_risk": "",
  "response_file": "path/to/feedback_response.json"
}
```

The manager does not infer whether another consumer should respond because of a
consumer response. Existing active items are served to a consumer until that
consumer has a `status: "signed_off"` response on the item. Raw spans are served
to a consumer until that consumer has a `status: "signed_off"` response on the
span.

Unlike spans, signed-off items can be reopened. When a consumer modifies an item
record itself, such as changing its span membership, active state, or item update
note, the response-intake helper marks prior `signed_off` item responses as
`reopened` before appending the modifying consumer's new response.

Use `status: "needs_followup"` only when migrating or preserving an older
consumer response that was not a current signoff. Such responses remain visible
to that consumer when feedback is served.

Use `status: "reopened"` only on item responses that were previously signed off
but were invalidated because the item itself was modified. Reopened item
responses remain in history, but they no longer suppress serving that item to
the consumer.

## Serving Feedback

Consumers receive all spans for which they do not have a signed-off response and
all active existing items for which they do not have a signed-off response. Do
not match by inferred target. Do not filter by tags. Do not decide that a
consumer is irrelevant.

Run:

```bash
.agents/skills/feedback-manager/scripts/serve_feedback \
  --consumer_name <consumer> \
  --tracking_root path/to/feedback_tracking/<feedback_id>
```

The served packet includes:

- a short README;
- machine-readable served metadata;
- every unsigned raw feedback span and raw excerpt when the raw source is
  available;
- every active existing item not already signed off by that consumer;
- a response template the consumer can fill in.

Do not show a span or existing item to the requesting consumer when that consumer
already has a signed-off response on it. Use `--include_signed` only for explicit
review or repair.

## Consumer Response Intake

Consumers submit response files with:

```bash
.agents/skills/feedback-manager/scripts/feedback_signoff \
  path/to/served_feedback_dir \
  path/to/feedback_response.json
```

The signoff script updates only `feedback_spans.jsonl` and `feedback_items.jsonl`.
It appends span signoffs, creates new consumer-owned items, updates existing
consumer-owned items, reopens item signoffs when an item is modified, and appends
`consumer_responses` entries. It does not write separate response, signoff,
routing, pending, or status-summary ledgers.

Consumers should not need to know the JSONL ledger schema.

The response file has four consumer-facing sections:

- `span_signoffs`: close out spans the consumer has inspected. Signed-off spans
  are not reopened.
- `new_items`: create a consumer-owned work unit from one or more spans.
- `item_signoffs`: append a completion, no-change, not-applicable, or follow-up
  response to an existing item without changing the durable item record.
- `item_modifications`: change an existing item's span membership, active state,
  or durable work-unit note. A nonempty reason and actual before/after change are
  required, and every modification reopens prior item signoffs.

An existing item may appear in only one of `item_signoffs` or
`item_modifications` in a single response.

Use `deferred_by_run_policy` with `status: "signed_off"` only when the routed
run policy authorizes continuation past a user decision. Cite the corresponding
planner blocker fingerprint in `residual_risk`. This records completion under
the current policy without changing the item or claiming that the underlying
decision was scientifically resolved.

## Validation Checks

Before finishing feedback-manager work, verify:

- `feedback_spans.jsonl` exists and is valid JSONL;
- if `feedback_items.jsonl` exists, it is valid JSONL;
- all `source_spans` in `feedback_items.jsonl`, when present, resolve to span IDs;
- each item, when present, has `created_by_consumer`, not
  `created_by_skill: "feedback-manager"`;
- all span `raw_source` paths exist unless an accepted exception is recorded;
- no span or item contains target-inference fields such as `weak_targets`,
  `target_confidence`, `locator_note`, or `item_scope_note`;
- each span `consumer_responses` value is a list;
- each `consumer_responses` value is a list;
- each consumer response has `consumer_name`, `responded_at`, `status`,
  `response_type`, and `response`;
- each response artifact path exists unless an accepted exception is recorded;
- checksums exist for immutable archive files when an archive packet was created.
