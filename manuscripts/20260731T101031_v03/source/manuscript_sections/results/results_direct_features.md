---
section_id: results_direct_features
title: "High ploidy preserves live-cell time while increasing glucose drawdown per live-cell hour"
primary_claims:
  - "C0: Ploidy modulates response to glucose starvation."
  - "C1: High ploidy increases glucose consumption rate."
  - "C2: High ploidy decreases yield (per cell)."
  - "C3: High ploidy cells are more robust to starvation."
main_figures: ["Figure 2"]
supplemental_figures: ["Figure S3", "Figure S4"]
input_versions:
  figure_description: "20260729 v03 figure-set semantic interpretations and final images for F2_a-c, S3_a, and S4_a"
  claim_graph: "20260730 v03 clean figure-level claim graph; user-fixed C0, C1, C2, and C3 reviewed"
---

## Results Text

The trajectories established that ploidy altered glucose-starvation responses, but they did not distinguish whether this reflected growth, persistence, glucose use per unit live-cell time, or cell production. We therefore extracted complementary features from the observed time courses. A representative trajectory illustrates the derivative tangent, live-cell AUC, and total-cell change over the glucose-measurement window (Figure 2B). The corresponding MDA-MB-231 profiles show how growth, glucose drawdown, and cell yield varied across starting glucose concentrations in low- and high-ploidy cultures (Figure 2C).

We first asked whether high ploidy preserved viable-cell time during the most severe glucose restriction. At 0 mM glucose, live-cell AUC was higher in high-ploidy cultures in all four cell lines in the standardized summary (Figure 2A). The raw all-line comparison showed the same ordering in all five cell lines, including SUM-159-fuse (Figure S3). Thus, despite heterogeneity in other trajectory features, these data show that high-ploidy cultures were more robust to starvation in the direct sense that they maintained more live-cell hours in the absence of supplied glucose.

We next quantified glucose consumption relative to the live-cell time available to consume it. For each cell line and ploidy state, we fitted an intercept-bearing regression of glucose drawdown against live-cell AUC at starting glucose concentrations of 1 mM or less; the slope measures glucose drawdown per live-cell hour. The fitted coefficient was higher in high ploidy in each of the four cell lines summarized in Figure 2A. These data show that high ploidy increased the fitted glucose-consumption rate per live-cell hour across the four-line primary comparison. This coefficient is a culture-level relationship between glucose drawdown and live-cell AUC, rather than a direct measurement of single-cell uptake or metabolic flux. Consistent with broader context dependence, SUM-159-fuse showed the opposite ordering (Figures S3 and S4).

The remaining features defined the limits of these shared effects. At 1 mM glucose, net total-cell yield was lower in high ploidy in all four Figure 2A comparisons, although the MDA-MB-231 values were nearly identical; SUM-159-fuse again reversed this ordering (Figure S3). Because this feature is an unnormalized change in total-cell number, it supports reduced cell output at fixed low glucose but does not by itself establish a per-cell or volumetric yield difference. Growth features were also heterogeneous: low-glucose and high-glucose growth effects crossed zero among cell lines, and the live-AUC intercept did not shift uniformly (Figures 2A and S3). Across the full glucose series, growth and yield profiles likewise varied by cell line and ploidy, and several yield trajectories peaked before the highest glucose condition (Figure S4). Taken together, the direct features identify a shared high-ploidy phenotype of greater live-cell persistence without glucose and, in four cell lines, greater culture-level glucose drawdown per unit live-cell AUC, while growth and total-cell output remain strongly context dependent.

## Revision Notes

1. The current figure descriptions support the direction of the fitted glucose-drawdown effects but do not establish the prior section's exact raw slope coefficients; those coefficients were therefore omitted.
2. C0, C1, and C3 are stated directly within the supported scope. The SUM-159-fuse reversal locally limits the generality of C1, and the prose distinguishes the culture-level glucose-drawdown coefficient from direct single-cell uptake or metabolic flux.
3. The 1 mM feature is unnormalized net total-cell gain. It supports lower cell-number output in the four-line primary comparison but does not directly establish the per-cell yield claim in C2 or a volumetric-yield effect.
4. Figure 2B is described as illustrating total-cell change rather than net total-cell yield because the current semantic interpretation identifies its max-delta annotation as the observed total-cell range within the glucose-measurement window.
