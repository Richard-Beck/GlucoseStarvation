# Implemented user-requested F5 revisions

Implemented in the seven-panel polished Figure 5 revision on 2026-07-29.

1. Apply the accepted signed pseudo-log glucose display to the strategy-path
   panel:
   `asinh(glucose / 1e-4)`, with an exactly regenerated zero-inclusive surface
   grid and display-only negative padding of the zero-boundary tiles.
2. Replace panel A's support-ranked strategy axis with a fixed scientific order:
   ascending initial glucose concentration, then strategy family in the order
   `XXX`, `XAA`, `XCC`. Continue to omit `0AA`.
3. Promote the three reviewed panel-A alternatives under
   `tmp/f5_panel_a_draft_options_20260729/` to panels b-d in the same left-hand
   column as panel a.
