---
name: feedback-manager
description: "Archive, track, serve, and record durable manuscript-feedback state across workflow stages. Use when the user explicitly requests feedback archival or tracking, an existing feedback archive, tracking root, or served packet is supplied, or a workflow explicitly enables managed-feedback mode. Do not use for ordinary current-turn instructions, inline comments, approvals, corrections, or immediate artifact revisions."
---

# Feedback Manager

## Purpose

Use this skill to keep manuscript feedback auditable across redrafting workflows.

The feedback-manager owns feedback state, not feedback meaning. It preserves raw
feedback, stores mechanical span records, serves unsigned spans and existing
consumer-created items to consumer skills, and records consumer-created or
consumer-updated items.

The feedback-manager must not infer feedback targets, summarize user intent,
classify directive strength, create feedback items, assign semantic ownership, or
decide what a consuming skill should do. Consumers inspect the raw feedback and
relevant manuscript context within their own scope. Only consumers create or
update feedback items.

## Feedback Consumer Orientation

To get managed feedback, pass its tracking root explicitly:

```bash
.agents/skills/feedback-manager/scripts/serve_feedback \
  --consumer_name <active-skill-name> \
  --tracking_root path/to/feedback_tracking/<feedback_id>
```

Do not infer or auto-select a tracking root.

If `--consumer_name` is omitted, the script reads `CODEX_SKILL_NAME` or
`CONSUMER_NAME`. When neither environment variable is available, pass the active
skill name explicitly.

The script writes a served-feedback directory containing:

```text
README.md
served_feedback.json
feedback_response_template.json
```

The consumer should read the raw feedback spans, existing items, and relevant
manuscript artifacts. If the consumer determines feedback requires work,
signoff, or a no-change/not-applicable record within its scope, it fills out the
response template to sign off spans, create a new item from one or more spans,
sign off an existing item, or explicitly modify an existing item, then submits
it with:

```bash
.agents/skills/feedback-manager/scripts/feedback_signoff \
  path/to/served_feedback_dir \
  path/to/feedback_response.json
```

Consumers should not edit `feedback_tracking/*.jsonl` directly. They only need the
served feedback packet, the response template, and the signoff script.

Consumer rules:

- A span is raw feedback. Sign off a span after deciding what, if anything, it
  requires in your scope. Signed-off spans are not reopened.
- An item is a consumer-created work unit linking one or more spans. Create an
  item only when the span requires action, signoff, or a documented
  no-change/not-applicable decision in your scope.
- Put stage completion, no-change, or not-applicable text in an item signoff's
  `response`; an item signoff cannot change the durable item record.
- Modify an existing item only when its span membership, active state, or durable
  work-unit note genuinely changes. Every actual item modification reopens prior
  item signoffs.
- Put each existing item in either `item_signoffs` or `item_modifications`, never
  both in one response.
- When a routed run policy explicitly defers a user decision, use a signed-off
  `deferred_by_run_policy` response and cite its blocker fingerprint in
  `residual_risk`; do not modify the item merely to record the deferral.
- Sign off an item when your scope is complete for the current item state.

## Feedback Manager Orientation

For archive intake, tracking schema, script behavior, and validation checks, read:

```text
.agents/skills/feedback-manager/references/manager_workflow.md
```

## Completion Standard

Finish feedback-manager work only after raw feedback is archived or explicitly
recorded, spans are mechanical records without target inference, and any
consumer-created or consumer-updated items have been recorded without manager
interpretation.
