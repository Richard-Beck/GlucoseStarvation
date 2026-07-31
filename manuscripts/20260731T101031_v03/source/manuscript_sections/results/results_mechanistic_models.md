---
section_id: results_mechanistic_models
title: "Richer resource models improve fits, while ploidy effects transfer unevenly"
primary_claims: []
main_figures: ["Figure 3"]
supplemental_figures: ["Figure S5", "Figure S6", "Figure S7", "Figure S8", "Figure S9"]
input_versions:
  figure_description: "20260729 v03 figure-set semantic interpretations for Figure 3A-C and Figures S5A-C, S6A-E, S7A-D, S8A-B, and S9A"
  claim_graph: "20260730 v03 clean figure-level integrated claim graph"
---

## Results Text

The direct trajectory features established that ploidy alters glucose-starvation responses, but they did not identify the dynamic processes that generate those differences. We therefore compared a family of ordinary differential equation models that linked measured glucose to growth and allowed additional structure through a hidden second resource, allocation between energy and biomass production, maintenance, and optional waste effects. Crossing four resource configurations with three waste modes produced 12 candidate structures spanning measured-glucose-only, two-resource, and waste-enabled models (Figure 3A).

Comparing deviance against effective parameter count supported a set of alternatives rather than a single resolved biochemical mechanism. Across the five Pareto-front candidates, deviance decreased as the effective parameter count k increased, from the 1R model at k = 24 to the flexible two-resource model with multiplicative waste at k = 66 (Figure 3B). The flexible two-resource waste model also had the smallest displayed delta AIC among the 12 assessed structures, whereas the 1R model had the largest (Figure S5B). This advantage was not explained simply by optimization stability: the displayed fraction of valid starts within 1 varied substantially among models and did not track Pareto membership (Figure S5A). The added complexity was nevertheless explicit and structured, arising from additional resource-kinetic, allocation, and waste-kinetic effects at the line and ploidy levels (Figure S5C). We therefore retained a mechanistically diverse five-model set for subsequent comparisons rather than selecting one latent mechanism as uniquely established.

We next asked whether the advantage of richer structures was confined to a small subset of conditions. Figure 3C reports condition-specific deviance loss per well relative to the best of the five retained models. This loss varied jointly with model, cell line, ploidy state, and starting glucose in the analysis excluding SUM-159-fuse (Figure 3C). The 1R model incurred substantial losses in several conditions, whereas the flexible two-resource candidates, with and without multiplicative waste, more often approached the best local fit, although each retained isolated regions of higher loss. Thus, hidden-resource and waste structure improved the simultaneous description of heterogeneous trajectories, but no candidate was uniformly best across every condition.

The fitted trajectory atlases and latent-state summaries show where these models differ while also defining the boundary of mechanistic interpretation. Across all five lines, the models produced distinct time courses in the blocks labeled N, R1, W1, and R2, with separations depending on ploidy and starting glucose (Figure S6A-E). Within the models, R1 denotes measured glucose, whereas R2 and W1 are model-implied states rather than directly measured environmental variables. Accordingly, the high-minus-low growth surfaces and phase-space displays visualize line-, model-, and glucose-dependent latent-state relationships, but they do not identify R2 as a particular nutrient or W1 as a particular waste product (Figure S8A; Figure S9A). The fraction-positive summaries were likewise heterogeneous across lines, with SUM-159-fuse displaying a visibly distinct, nearly uniformly positive surface (Figure S8B). These latent-state patterns provide hypotheses for how richer models encode trajectory heterogeneity, not direct measurements of the proposed hidden processes.

We next tested whether fitted ploidy effects generalized to an opposite-ploidy state withheld from a target line in arbitrary transfer contexts. The analysis compared held-out log likelihood under a fitted ploidy effect with a null model lacking that effect (Figure S7A). Transfer was not uniformly successful: most line-model-direction scores were negative, although outcomes varied strongly among lines and models (Figure S7B). The highlighted cases illustrate this heterogeneity. MCF10A low-to-high transfer with the 2R(f),W(m) model improved on the null by approximately 494 log-likelihood units, whereas SNU668 high-to-low transfer with the 2R(f) model was worse by approximately 802 units (Figure S7B-D). Thus, a fitted ploidy effect can improve prediction in particular contexts but is not a reliable universal transfer rule.

Taken together, the model comparison supports hidden-resource and waste structure as a better description of glucose-starvation dynamics than measured glucose alone. The latent variables remain mechanistic hypotheses, and held-out transfer shows that fitted ploidy effects can generalize in selected contexts but are not universally transportable.

## Revision Notes

1. No `user_fixed: true` claim is primarily owned by Figure 3 or Figures S5-S9, so `primary_claims` remains empty. Figures S6 and S8 provide secondary model-based support for the locked claim that ploidy modulates glucose-starvation responses.
2. Figure S7 scores are aggregate held-out log-likelihood differences without per-well or per-observation normalization. Its highlighted curves are deterministic trajectories from selected optimization draws, not posterior uncertainty summaries.
3. R2 and W1 are treated as model-implied latent states. Figures S6, S8, and S9 do not establish their biological identities, so these panels support bounded, hypothesis-generating interpretation rather than direct identification of physical resources or waste products.
4. Figure S14 is omitted because its current content is outside this section's assignment and no longer supports the v02 prediction-transfer discussion.
