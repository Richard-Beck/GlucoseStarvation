# Manuscript Figure Style

Use this shared reference when drafting, polishing, or visually reviewing glucose-starvation manuscript figures. Workflow-specific skills may add stricter requirements, but should not contradict this file without documenting a user-approved exception.

## Figure Images

- Keep review and final figure deliverables PNG-first unless a workflow explicitly requires another format.
- Do not put figure-level titles, subtitles, captions, result sentences, or `Figure N` / `Figure SN` headers on figure images. Put manuscript prose in separate legends or notes.
- Axis labels, facet/strip labels, scale bars, colorbars, legends, and direct annotations that are part of a data panel are allowed.
- Size journal-style final figures at about 3.25 to 3.5 inches wide for single-column figures and about 7 inches wide for double-column or full-page figures. Maximum final height should generally be 9 to 9.5 inches.
- Export final color figures at 300 dpi or higher, mixed raster/line figures at 600 dpi or higher, and black-and-white line-art figures at 1000 dpi or higher when relevant.
- Use 8 to 12 point Arial, Helvetica, or Times New Roman text in final figure images unless a project-specific upstream plotting script makes a documented exception necessary.

## Panels And Layout

- Label every displayed subpanel with lowercase labels: `a`, `b`, `c`, and so on.
- Keep panel labels outside plotted data regions where possible, and use consistent size, weight, placement, and ordering across a figure.
- Arrange panels in the intended reading order. Avoid layouts that make panel order ambiguous.
- Keep spacing, gutters, margins, and alignment deliberate and compact enough for manuscript use without crowding labels or plotted content.
- Avoid excessive in-image prose and redundant evidence. If interpretation needs a sentence, prefer the legend or manuscript text.

## Color And Annotation

- Prefer palettes that separate scientific groups clearly and remain readable in grayscale or for common color-vision deficiencies.
- Define all colors, symbols, line styles, abbreviations, error bars, and statistical notation in the legend.
- Keep direct annotations close to the data feature they explain and short enough that they do not dominate the panel.
- Use consistent color, line, symbol, and annotation meanings across related figures unless a change is explicitly documented.

## Figure Legends

- Store legends separately from figure images.
- Begin each legend with `Figure N,` or `Figure SN,` followed by a succinct descriptive title.
- Include a brief overview of what is shown, including method, sample size, subject/context, and conditions where applicable.
- Describe each panel label and define all labels, symbols, colors, line styles, and abbreviations used in the figure.
- Keep legends concise, usually three or fewer sentences unless the figure is complex.
- Use sentence case and avoid moving full Results or Discussion interpretation into the legend.

## Visual QC

Before marking figure outputs ready for the next workflow stage, inspect the rendered PNGs directly. Check for:

- visible figure-level titles, subtitles, captions, or headers
- clipped panel labels, axis text, legends, colorbars, scale bars, or plotted data
- unreadable text or symbols at intended print size
- overlapping plot elements, labels, or annotations
- awkward or excessive spacing
- incorrect panel order or inconsistent panel labels
- excessive density, excessive in-image prose, or redundant panels
- reused raster material that preserves known defects

If a serious review-visible defect remains, regenerate an alternate, revise the plotting/layout code, or document the defect as a blocker or user-approved exception.
