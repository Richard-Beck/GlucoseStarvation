---
section_id: results_selection_simulations
title: "Fitted simulations predict that starvation selects high ploidy and sustained abundance selects low ploidy"
primary_claims:
  - "C0: Ploidy modulates response to glucose starvation."
  - "C3: High ploidy cells are more robust to starvation."
  - "C4: Starvation selects for high ploidy."
  - "C11: Abundance selects for low ploidy."
main_figures: ["Figure 5"]
supplemental_figures: ["Figure S11", "Figure S12"]
input_versions:
  figure_description: "20260729 v03 figure-set manifest and semantic interpretations for F5a-g, S11a, and S12a-f"
  claim_graph: "20260730 v03 clean accepted claim graph; user-fixed C0, C3, C4, and C11"
---

## Results Text

We next asked whether ploidy-specific growth responses would translate into selection when low- and high-ploidy populations shared the same modeled resource. We initialized the simulated populations at equal fractions and compared 20 displayed six-day glucose schedules spanning seven starting concentrations and refresh, add-glucose, or carry actions at days 2 and 4. Selection was summarized by posterior support and the posterior median log ratio for high-ploidy abundance exceeding low-ploidy abundance. These are predictions from fitted mixed-population simulations, not measurements from empirical competition assays.

At day 6, posterior support for high-ploidy enrichment ranged from approximately 25% to 94% across glucose schedules and seeding levels. The zero-glucose carry schedule consistently favored high ploidy, with approximately 82% support at all three seeding levels, whereas many other schedules moved from below 50% support at 0.5× seeding to above 50% at 1× and 2× seeding (Figure 5A). The direction also changed with time and resource handling. Support under zero-glucose carry remained approximately 82–84% from days 2 to 6, while support under the 25-mM refresh schedule remained below 50% at day 6. In contrast, the 25-mM add-glucose and carry schedules rose from less than 20% support at day 2 to approximately 60% at day 6 (Figure 5B). Thus, sustained starvation favored high ploidy and maintained glucose abundance favored low ploidy, but high starting glucose alone did not determine the outcome.

The predicted response also depended on model structure and cell line. At 1× seeding on day 6, model-specific support across the schedule set ranged from roughly 90–130 to nearly 400 of 400 draws, although several models approached maximal support under the zero-glucose strategies (Figure 5C). MCF10A and MDA-MB-231 frequently approached 500 of 500 supporting draws across models, whereas SNU668 was generally shifted toward low-ploidy enrichment and fell below the 250-draw midpoint for multiple schedules (Figure 5D). Independent single-line fits likewise produced line-specific high-to-low posterior contrasts in yield, uptake, resource response, and maintenance or death, with large opposing parameter shifts in MCF10A and smaller or reversed shifts in several other lines (Figure S11A). These differences identify cell line and model specification as material determinants of the simulated selection response.

The fitted growth landscapes provided a mechanistic basis for this resource dependence. Posterior-median high-minus-low instantaneous growth was typically positive in low-glucose, lower-latent-resource regions and negative in higher-glucose, higher-resource regions (Figure 5E). Three representative schedules traversed these landscapes differently across cell lines and models, generating ploidy-specific separations in both live-cell abundance and instantaneous net growth over time (Figure 5F,G). The predicted enrichment therefore arose from the interaction between ploidy-specific growth and the resource trajectory imposed by each schedule, rather than from starting glucose concentration alone.

Comparing fit contexts preserved the broad resource-dependent reversal while exposing its limits. Across four-line, all-line, and single-line fits, low-glucose schedules often enriched high ploidy, whereas high-glucose schedules often enriched low ploidy at days 4 and 6 (Figure S12A–F). However, these directions were not universal: SNU668 frequently favored low ploidy even under low glucose, while SUM-159-fuse and many MDA-MB-231 configurations retained high-ploidy enrichment under high-glucose schedules. Taken together, the fitted simulations predict that starvation selects for high ploidy and sustained resource abundance selects for low ploidy, while showing that both predictions depend on schedule, time, seeding density, cell line, model structure, and fit context.

## Revision Notes

1. C4 and C11 are stated directly as model-supported predictions and are explicitly distinguished from empirical competition results. Current S12 counterexamples are retained locally rather than averaged into a single cross-context direction.
2. Figure S11A does not measure selection and is used only to establish line-specific heterogeneity in the fitted high-to-low parameter contrasts.
3. The prior 448-schedule objective ranking, exact selected-schedule endpoints, arbitrary retained-model filter, raw-score fallback, leave-one-line objective analysis, and starting-composition claims are absent from the current F5/S11/S12 semantic package and were therefore removed.
4. Approximate numerical support values are included only where the current semantic interpretations and accepted figure audit independently resolve them.
