# Independent visual QC

Status: **PASS** after source-level correction and complete rerender on
2026-07-29.

- `figure_1_all_timecourses.png`: PASS. Panel order a→b is clear; no clipping
  or overlap; dense facets remain readable at the 7-inch target width.
- `figure_2_confluence_robustness.png`: PASS. Panel order a,b / c,d is clear;
  axes and legends are readable; vertical facet strips are tight but legible.
- `figure_3_focused_distributions_and_same_2n.png`: PASS. Panel order a,b / c,d
  is clear; legends, axes, and ratio annotations are readable and unclipped.
- `figure_4_cytoplasmic_signal_and_multimodal_fields.png`: PASS. Panel order
  a,b,c over d is clear. The first independent inspection found clipping in the
  rightmost microscopy-column header. The header was wrapped in the R
  constructor and the complete subpanel/layout/final chain was rerun. The
  refreshed header and all row labels are readable and unclipped.

Raster-policy audit: PASS. The final constructor does not read or copy an
approved-reference, audit-subpanel, or final PNG. Figure 4d reads raw channel
and mask TIFFs, then computes crops, normalization, outlines, and composites in
memory. All four finals are assembled from live panel objects.

Non-blocking inherited style notes: a small amount of facet/annotation text is
6.5–7.5 pt, and the microscopy cyan/amber class outlines differ from the
blue/red quantitative encoding. Both remain legible and are explicitly
defined.
