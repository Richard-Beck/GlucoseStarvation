# Fresh-context figure audit prompt

Use this template for one figure at a time. Replace bracketed placeholders with
the supplied inputs.

---

Perform a strict clean-room semantic assessment.

Do not use or read any skill, prior claim graph, manuscript, project
documentation, feedback, figure code, source data, lineage record, provenance
outside the supplied semantic descriptions, or any other file. Do not launch
subagents. Do not write or edit files.

You may read only these semantic descriptions for figure `[FIGURE_ID]`:

`[SEMANTIC_INPUTS_OR_EXCERPTS]`

Assess only these authoritative claims:

`[AUTHORITATIVE_CLAIMS]`

Identify every distinct observation in the supplied figure that supports or
undermines one or more authoritative claims.

There is no qualification category.

- If an observation genuinely supports a claim, use `support`.
- If it genuinely conflicts with a claim, use `undermine`.
- If it merely adds context, scope, uncertainty, heterogeneity, or a caveat
  without supporting or undermining the claim, omit it.
- Do not invent claim wording or infer an unstated biological identity.
- Do not treat model output as empirical observation. Model-derived evidence may
  support or undermine a model or prediction claim when the supplied semantics
  justify the relationship.
- Cite the exact contributing panel or panel range.

Use these strength labels:

- `strong`: direct and clear evidence for the stated claim.
- `moderate`: meaningful evidence with a material limitation in directness or
  scope.
- `weak`: a genuine but limited relationship.

Return only:

### Observation–claim relationships

| Panel(s) | Exact observation | Claim | Relation | Strength | Reason |
|---|---|---|---|---|---|

If no material relationships are present, return the heading and table header,
then write:

`No material relationships identified.`

Do not add an empty placeholder row. Do not use vertical-bar characters inside
table cells.

---
