# Subpanel Semantic Interpretation Workflow

## Purpose

Interpret manuscript subpanels from independent visual reads while keeping
source clarification separate. Every subpanel requiring a new interpretation
must receive its own fresh-context visual inspector. The main agent coordinates
inspectors, performs source clarification, and assembles the final artifacts.

## Mandatory visual-inspector requirement

Launch exactly one fresh-context visual-inspector subagent per subpanel that
lacks an authenticated reusable visual interpretation. Never assign multiple
subpanels to one inspector and never replace the inspector with main-context
image interpretation.

An unchanged subpanel may reuse a prior visual interpretation only when the
displayed pixels, visual-inspector prompt contract, panel identity, and retained
output are authenticated as unchanged.

If callable subagents are unavailable or prohibited, stop and tell the user that
fresh subpanel interpretation cannot proceed.

## Main-agent workflow

1. Enumerate subpanels from user instructions, the canonical figure index, or
   other non-semantic routing metadata.
2. Freeze each target's stable figure/panel ID and exact rendered image.
3. Build one minimal visual-inspection packet per subpanel.
4. Launch one fresh-context visual inspector per packet. The main agent launches
   inspectors directly; do not create an intermediate worker-agent layer.
5. Preserve every inspector return verbatim and check only its identity and
   output shape before semantic use.
6. After all visual reads are frozen, perform source clarification in the main
   context. Group clarification across panels that share a figure, generator,
   dataset, or encoding rather than repeating the same trace per subpanel.
7. Launch an optional clarification subagent only when a bounded ambiguity
   genuinely benefits from separate inspection. Clarification delegation is
   never required per subpanel and must not change the frozen visual read.
8. Assemble one interpretation artifact per subpanel, keeping visual material
   and source clarification in separate sections.

## Visual-inspection packet

Give the inspector only:

- the stable figure/panel ID;
- the rendered subpanel as an image item; and
- when necessary, a supplied visual-only shared legend or key image needed to
  decode marks visible in the target.

Do not give it repository paths, filenames, skills, prior conversation, figure
semantics, source code, provenance, legends as prose, manuscript claims,
feedback, or drafting history. Tell it not to use skills, browse files, or launch
subagents.

Use this prompt contract:

```text
Inspect only the supplied subpanel pixels. Do not infer scientific identity,
correctness, provenance, workflow intent, or manuscript-claim support. Return
only the requested Markdown structure. Keep observations concrete and concise;
report unreadable or ambiguous elements instead of guessing.

### Target
- panel_id: <supplied id>

### Visible structure
- display type:
- axes and units:
- groups or facets:
- color, shape, line, or image encodings:
- visible annotations:

### Visible observations
1. <visible pattern, ordering, difference, trend, or spatial feature>
2. <visible observation, if present>
3. <visible observation, if present>

### Visual limitations
- unreadable or ambiguous elements:
- overplotting, clipping, cropping, or resolution limits:
- comparisons not resolvable from the pixels:
```

Inspectors may return fewer than three observations. They must not pad the
output or convert visible patterns into mechanistic or manuscript-level claims.

## Main-agent source clarification

Start only after the relevant visual-inspector returns are frozen. Inspect only
the bounded existing artifacts needed to address:

- **Visual-observation checks:** quantify, qualify, or falsify observations made
  by the inspector using the plotted or source data where recoverable.
- **Visible-element resolution:** resolve colors, shapes, labels, units,
  denominators, transforms, facets, abbreviations, or shared keys that were
  visually ambiguous.
- **Provenance chain:** identify plotted values, plotting code, generated tables,
  transformations, filters, and upstream sources as far as the available
  package supports.

Record `not determined` rather than guessing. Do not rerun analyses, fit models,
judge raw-data correctness, inspect raw feedback, or perform manuscript-level
evidence review under this module.

Source clarification may test an inspector observation, but it does not become
visual evidence. Keep corrections and quantitative results in the clarification
section.

## Per-subpanel output contract

Write one compact artifact containing:

- **Target:** stable figure/panel ID and exact image identity.
- **Fresh-Context Visual Description:** the inspector return verbatim.
- **Visual Claims:** concise statements supported only by that return, with a
  short cited excerpt or tight paraphrase.
- **Visual-Only Interpretation:** the apparent message supported by visible
  axes, labels, marks, annotations, and patterns.
- **Source Clarification:** visual-observation checks, visible-element
  resolution, and provenance chain, or `not requested`/`not determined`.
- **Visual Caveats:** visible ambiguities and limitations.
- **Not Assessed:** raw feedback, prior interpretations, manuscript prose,
  source-data correctness, rerun analysis, and manuscript-claim support.

Use this form for visual claims:

```markdown
**V1. [Visual claim stated plainly.]**
Support: fresh visual inspector observed "[short quote or tight paraphrase]."
```

Semantic coverage remains canonical in `figure_set_manifest.csv`. A separate
interpretation index may be generated as a derived review view.

## Hard rules

- Require one independent fresh visual inspector for each newly interpreted
  subpanel.
- Do not insert a worker agent between the main coordinator and the inspector.
- Do not combine multiple subpanels in one inspector assignment.
- Do not ask the inspector to read skills or repository material.
- Preserve inspector output verbatim.
- Do not create new visual observations from main-context image inspection.
- Do not present source clarification as visual evidence.
- Keep source clarification main-agent-owned by default and group shared traces.
- Use clarification subagents only for bounded ambiguities, never by default per
  subpanel.
- Do not claim visible labels, legends, or annotations are correct; treat them
  as visible text until clarified separately.
