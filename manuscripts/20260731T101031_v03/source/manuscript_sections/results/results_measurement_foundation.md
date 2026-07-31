---
section_id: results_measurement_foundation
title: "Ploidy generally increases cell size and alters starvation trajectories in a context-dependent manner"
primary_claims:
  - "C0: Ploidy modulates response to glucose starvation."
  - "C3: High ploidy cells are more robust to starvation."
  - "C7: Cell size is proportional to ploidy."
main_figures: ["Figure 1"]
supplemental_figures: ["Figure S1", "Figure S2", "Figure S15", "Figure S16", "Figure S17", "Figure S18"]
input_versions:
  figure_description: "20260729 v03 integrated semantic interpretations for F1, S1, S2, and S15-S18"
  claim_graph: "20260730 v03 clean figure-level integrated claim graph"
---

## Results Text

To determine how ploidy changes the response to glucose limitation, we profiled matched lower- and higher-ploidy states from five cell-line backgrounds across a range of starting-glucose concentrations. Time-resolved microscopy and extracellular-glucose measurements provided paired live-cell, dead-cell, and glucose trajectories, while object segmentation and classification converted the images into cell counts (Figure 1a; Figure S1a). Manual alive-point counts closely tracked manual alive-object labels: 64 of 90 validation frames matched exactly, and 82 differed by no more than two objects (Figure S2a). Target-scoped classifier performance supported live-cell calling (ROC AUC = 0.86; balanced accuracy = 0.78; precision = 0.86; recall = 0.91; Figure S2b). Measurement coverage was nearly complete across cell line, ploidy, and glucose conditions; censoring was confined to six SUM-159-fuse condition-by-glucose combinations (Figure S2c,d). Glucose calibration slopes were positive in every calibration set, although residual error varied and was greatest for SUM-159-fuse (Figure S2e-g). Together, these checks established an adequate measurement framework for comparing starvation responses while identifying SUM-159-fuse glucose measurements as the principal calibration-sensitive subset.

We next asked whether the matched ploidy states exhibited the expected relationship between ploidy and cell size. Figure 1b showed clearly larger median cell area at higher ploidy in three of five backgrounds: MCF10A, SNU668, and SUM-159-chem. Nuclear area generally also increased with ploidy, and condition-level changes in cell and nuclear area were positively associated (Figure 1c,d). MDA-MB-231 had only a minor ploidy difference and correspondingly little cell-area separation, whereas SUM-159-fuse was the exception of interest: despite an approximately twofold ploidy contrast, its labeled high/4N cells were smaller than its labeled low/2N cells (Figure 1b-d). Thus, cell size generally increased with ploidy, while SUM-159-fuse uncoupled higher ploidy from a larger-cell phenotype.

The SUM-159-fuse exception complicated a straightforward biological interpretation and raised a sample-identity concern. Across the shared 24-48 h and 0.1-1 mM glucose comparisons, fuse C00 and I00 low/2N cells were larger than their high/4N counterparts, whereas SUM-159-chem showed the opposite ordering; this reversal persisted across overlapping field cell counts and total segmented areas (Figures S15a, S16a,c, and S17a). Moreover, the nominally identical fuse and chem low/2N populations were phenotypically inconsistent: fuse low/2N median cell area was 1.50-1.83-fold that of chem low/2N, although nuclear-area differences were weaker and nonuniform and therefore did not support a clean label swap (Figures S15b, S16b,d, and S17b-d). Independently, a cytoplasmic-red measurement of unresolved biological meaning had opposite high-versus-low ordering in fuse and chem, and this reversal remained after cell-size stratification (Figure S18a-c). These anomalies do not establish a swap, but they make an unresolved label or sample-identity problem plausible. SUM-159-fuse is therefore generally withheld from main-text analyses, with supplementary analyses testing whether its exclusion changes the conclusions.

Across the four retained backgrounds, live-cell counts generally accumulated when starting glucose was abundant, peaked and declined at intermediate glucose, and declined without added glucose, while extracellular glucose declined over time (Figure 1e,f). Higher ploidy improved persistence in some starvation trajectories, but the direction and magnitude of ploidy separation nevertheless varied among cell lines and glucose conditions (Figure 1e,f). Taken together, the measurements show that ploidy modulates growth and glucose-starvation responses, while the direction and magnitude of the raw phenotype depend on cellular background and nutrient context.

## Revision Notes

1. Locked claim C3 is stated with a context qualifier because Figure 1e supports improved higher-ploidy persistence in some starvation trajectories but also shows that the direction and magnitude of separation vary by cell line and glucose condition.
2. Locked claim C7 is stated as strongly as the assigned evidence permits but is locally qualified because the figures establish a general positive association rather than exact proportional scaling; MDA-MB-231 has little cell-area separation and SUM-159-fuse reverses the expected ordering.
3. Figures S15-S18 support a plausible label or sample-identity problem through opposing fuse-versus-chem cell-area and cytoplasmic-red orderings, but do not prove a swap or identify which lineage is affected. Nuclear-area evidence is weaker and nonuniform, and the mixture classifier is not an error-free reference.
4. The manuscript policy now states explicitly that SUM-159-fuse is generally withheld from main-text analyses and that supplementary analyses test whether this exclusion changes conclusions.
5. Figure 1f establishes culture-level glucose depletion but does not estimate or test a uniform per-cell glucose-consumption rate; no such claim is made here.
