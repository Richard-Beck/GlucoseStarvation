---
name: archived-subpanel-semantic-interpretation
description: "Archived historical copy of subpanel-semantic-interpretation, superseded by manuscript-figure-workflow on 2026-07-06. Use only when explicitly auditing or comparing the pre-promotion subpanel interpretation workflow."
---

# Subpanel Semantic Interpretation

## Purpose

Use this skill to interpret manuscript subpanels while keeping visual interpretation separate from source clarification.

The discipline is controlled delegation. Visual claims must come only from a fresh-context visual inspector. Panel-scoped source clarification may be added only after the visual read, by a separate clarification subagent, and must stay separate from the visual interpretation.

## Roles

Choose the role before doing any figure work:

- **Overseer:** default when asked to interpret multiple subpanels. Locate subpanel targets, choose a sensible output location, launch one worker subagent per subpanel, and collect their output paths. Do not inspect images, read rendered figures, or interpret panels.
- **Worker:** default for a single subpanel, or when an overseer tells you that you are a worker. Run the worker workflow below for exactly the assigned subpanel.

## Mandatory Subagent Gate

Before applying either role, verify that a callable subagent / multi-agent tool is available.

If no subagent facility is available, stop immediately and tell the user that this skill cannot be applied because it requires delegated fresh-context work. Do not substitute your own image inspection, main-context visual read, or self-contained guess.

## Overseer Workflow

1. **Enumerate targets**
   - Identify the requested subpanels and labels from user instructions, file listings, manifests, or other non-visual routing information.
   - Do not open image viewers, inspect rendered pixels, infer scientific meaning from filenames, or read figure content.

2. **Choose output locations**
   - Use a user-specified output path if supplied.
   - Otherwise create or select a concise folder near the requested figure package, such as `subpanel_interpretations/`, with one worker output file per subpanel plus an optional index.
   - Ensure output names identify figure/panel labels without encoding interpretive conclusions.

3. **Launch workers**
   - For each subpanel, launch a subagent with no forked conversation context.
   - Instruct the subagent to use this `subpanel-semantic-interpretation` skill and tell it: "You are a worker."
   - Pass only the target image path, user-supplied figure/panel label, requested output path, and any non-interpretive routing constraints.
   - Tell the worker that any subagents it launches must not use skills.

4. **Collect results**
   - Wait for each worker.
   - Verify that each expected output file exists and appears to follow the worker output contract.
   - Do not assess whether the visual interpretation is right. The overseer may check completeness, target labels, paths, and obvious missing sections only.

## Worker Workflow

1. **Freeze the target**
   - Record the exact image path and any figure/panel label supplied by the user.
   - Use the path only as an identifier. Do not infer meaning from directory names, filenames, dates, or package names.

2. **Fresh-context visual read**
   - Launch an image-inspection subagent with no forked conversation context.
   - Tell the image inspector not to use any skills and not to launch further subagents.
   - Give it only the target rendered image, preferably as a local image item. Do not give surrounding paths except what is technically required to load the image.
   - Tell it not to read or infer from surrounding project files, filenames, legends, code, metadata, raw feedback, prior conversation, or manuscript context.
   - Ask for a detailed description of visible pixels only: axes, axis labels, tick labels, legends visible inside the image, colors, marks, shapes, approximate positions, trends, annotations, panel letters, and visual caveats.
   - Ask it not to make scientific, workflow, data-lineage, or manuscript-level interpretations beyond what can be inferred from the image alone.
   - Wait for the image inspector to return before launching any clarification subagent.
   - Preserve the image inspector output verbatim in the interpretation record.

3. **Panel-scoped source clarification**
   - Launch a second subagent only after the visual read is complete.
   - Tell the clarification subagent not to use any skills and not to launch further subagents.
   - Give it the target path, the user-supplied panel label, and the fresh-context visual description returned by the image inspector.
   - Ask it to produce exactly three clarification sections:
     - **Visual-observation checks:** for trends, correlations, group differences, orderings, or other observations made by the image inspector, inspect existing plotted/source data where recoverable and quantify, qualify, or falsify the observation. Include the source artifact and method used for each check.
     - **Visible-element resolution:** clarify missing, ambiguous, or visually unresolved elements directly observed by the inspector, such as color/shape/line mappings, labels, units, denominators, transforms, facets, or abbreviations.
     - **Provenance chain:** trace the plotted values as far as possible through existing local plotting code, manifests, provenance files, generated tables, intermediate files, algorithms, transforms, filters, and upstream source datasets. Enumerate each recoverable derivation step, not only the final plotted table.
   - Prohibit raw user feedback, prior interpretations, manuscript drafts, broad manuscript context, rerunning analyses, new model fitting, or judging whether the underlying data are correct.
   - Require explicit "not determined" entries for any observation check, visible element, or provenance step it cannot identify from existing local artifacts.
   - Preserve the clarification output separately from the visual read.

4. **Extract visual claims**
   - State only claims supported by the image inspector's visual description.
   - Distinguish visible text/encoding from inferred data identity. For example, an axis label is visible; whether the plotted values correctly implement that label is not assessed.
   - Use cautious language for patterns that are approximate from the rendered image, such as "appears", "visually", "roughly", or "the visible points suggest".
   - Do not use the clarification subagent's output as support for a visual claim.

5. **Interpret visually, not evidentially**
   - Interpret the apparent visual message of the panel only from visible axes, labels, marks, annotations, and patterns.
   - Keep quantitative checks, visible-element resolutions, provenance, source identity, transformations, denominators, sample scope, and statistical summaries in the separate source-clarification section.
   - Do not make claims about workflow intent, drafting history, panel swaps, correctness of labels, or support for manuscript claims.

6. **Name limits explicitly**
   - List visual ambiguities, unreadable labels, overplotting, missing definitions, and assumptions that cannot be checked visually.
   - State that source clarification was panel-scoped and that no raw feedback review, manuscript-level review, rerun analysis, or source-data correctness audit was performed.

## Output Contract

### Overseer Output

Return a compact coordination report with these sections:

- **Role:** Overseer.
- **Targets:** subpanel labels and image paths assigned to workers.
- **Output Location:** folder or files where worker outputs were written.
- **Worker Results:** one line per subpanel with worker status and output path.
- **Not Assessed:** state that the overseer did not inspect figures or interpret panels.

### Worker Output

Return or write a compact interpretation with these sections:

- **Target:** exact image path and user-supplied figure/panel label, if any.
- **Role:** Worker.
- **Fresh-Context Visual Description:** verbatim image-inspector output or a clearly marked excerpt plus path to full retained output if saved.
- **Panel-Scoped Source Clarification:** verbatim clarification-subagent output or a concise extract with the required subsections: Visual-observation checks, Visible-element resolution, and Provenance chain. Clearly mark unresolved items as "not determined".
- **Visual Claims:** concise claim-led bullets about visible axes, labels, encodings, marks, annotations, and apparent trends. Each claim must cite a short excerpt from the fresh visual description.
- **Visual-Only Interpretation:** cautious paragraph or bullets describing the apparent message of the panel from pixels alone. Use phrases such as "visually appears to show" or "the visible pattern suggests".
- **Clarification Notes:** terse bullets derived only from the clarification subagent, clearly separated from visual claims. Include any quantified visual-observation checks, resolved visible elements, and provenance-chain summary.
- **Visual Caveats:** visible limitations and ambiguities.
- **Not Assessed:** a standard note that no raw user feedback, prior interpretations, manuscript drafts, source-data correctness validation, rerun analysis, or manuscript-level evidence review was performed.

For visual claims, use this shape:

```markdown
**V1. [Visual claim stated plainly.]**
Support: fresh visual inspector observed "[short quote or tight paraphrase]."
```

## Hard Rules

- If asked to interpret multiple subpanels, act as an overseer by default.
- Overseers must not inspect figures, describe visual content, or interpret panels.
- Workers must not skip the fresh-context image-inspection subagent.
- Workers must wait for the image inspector before launching the source-clarification subagent.
- Workers must tell all subordinate subagents not to use skills.
- Workers must not directly inspect the image in the main context.
- Do not infer visual meaning from path names, filenames, package names, or prior conversation.
- Do not present provenance, code, legends, comments, filenames, or drafter summaries as evidence.
- Do not mix the source-clarification addendum into visual claims or visual interpretation.
- Quantitative checks are allowed only for observations made by the fresh-context visual inspector and only from existing plotted/source artifacts.
- Provenance tracing must stay panel-scoped and go as far as possible through existing local artifacts; explicitly mark unresolved gaps.
- Do not rerun analyses, fit new models, validate raw source-data correctness, or perform manuscript-level evidence review under this skill.
- Do not claim that visible labels, legends, or annotations are correct; treat them as visible text only.
- Do not make manuscript-level claims from this skill's output. Use the output as a visual and source-clarification input to a separate evidence audit if needed.
