# Results Style Guide

Use this guide before drafting Results prose, and read it again before final revision. The goal is to reproduce the argumentative style: conclusion-first headings, direct framing of the biological or modeling question, compact methodological context where needed, figure-linked evidence, and synthesis sentences that state what the evidence implies.

## Write conclusion-first headings

Results headings should state the main finding, not merely name the topic or assay. A heading such as:

> Multiple-hitting CTLs exhibit heterogeneous late onset killing

is stronger than a topic label such as “Simulation of CTL killing,” because it tells the reader what the section establishes. Likewise:

> IFNγ-mediated cell-cycle arrest is sufficient for tumor control

works because it makes a specific, testable claim. The heading should be the shortest accurate version of the result’s central conclusion.

A good heading should answer: *What should the reader know after reading this section?*

## Open with the question the figure resolves

Do not begin by walking through panels. Begin by explaining why the experiment, model, or analysis was needed. The opening sentence should identify the uncertainty, hypothesis, expectation, or biological problem that motivates the result.

For example:

> We first sought to establish whether the multiple-hitting hypothesis was a feasible explanation for the heterogeneous, delayed onset, “burst” killing kinetics observed...

This sentence frames the result as an answer to a specific mechanistic question. Similarly:

> In most circumstances, the fittest karyotype is expected to eventually dominate the population. However, we asked whether this outcome could be altered simply by changing the missegregation rate.

This opening establishes the default expectation, then identifies the reason to test it. The figure becomes necessary because it resolves a tension.

A strong Results section usually begins by saying, in effect: *Here is the expectation or uncertainty; here is what we tested.*

## Include methods only when they sharpen interpretation

Results prose should include enough methodological context for the reader to understand the logic of the evidence, but it should not become a Methods section. Methodological details belong in the Results only when they affect interpretation.

For example:

> In these Monte Carlo simulations, CTLs hit targets at a constant rate λ, then targets died after receiving η hits.

This is not procedural clutter. It defines the model in the minimum terms needed to interpret the result. Likewise:

> We built approximate transition matrices over the charted region of each ALFA-K–fitted landscape and determined the steady state distribution of karyotypes across missegregation rates.

This sentence gives the reader the analytical premise without describing every implementation detail.

A useful test is: *Would omitting this method detail make the result harder to understand or easier to misinterpret?* If yes, include it. If not, leave it for Methods.

## Make each panel serve the argument

Do not describe figures in an A-to-Z sequence unless the argument genuinely requires it. Each panel should be used for what it establishes: setup, measurement, contrast, mechanism, sensitivity, exception, validation, or synthesis.

For example, a panel citation should support a claim:

> When η = 1, the hazard experienced by contacted targets does not change with time. In contrast, when η>1 the hazard experienced by contacted targets increases over time...

The point is not merely that a panel exists. The point is that the panel establishes a difference in hazard structure between single-hit and multiple-hit killing.

Similarly:

> There was no substantive difference between the single-hit and multihit scenarios in terms of tumor growth (Fig. 4A), or number of CTLs (Fig. 4B).

Here, the cited panels support a specific negative result. The sentence is useful because it tells the reader what comparison was made and what conclusion follows.

Panel citations should function as evidence anchors, not as coverage markers.

## Prefer interpretation over visual description

Results prose should rarely dwell on visual encodings, annotations, colors, boxes, overlays, or layout. Those details matter only when they are needed to understand the result. The prose should say what the figure establishes, not simply what it looks like.

A weak sentence says:

> The lower panel shows a bimodal distribution.

A stronger sentence says:

> In the latter case, a bimodal distribution occurred, which could be interpreted as a subpopulation of high-rate killers, yet importantly such a population did not exist in our simulations.

The stronger sentence explains the interpretive significance of the visual pattern. It also prevents a misleading conclusion. The reader learns not only what appeared in the figure, but what that appearance means.

## Build paragraphs as steps in an argument

Each paragraph should move the reader’s understanding forward. Good Results prose often progresses through a sequence: rationale, setup, comparison, exception, mechanism, synthesis.

Transitions such as “first,” “next,” “therefore,” “however,” “in contrast,” and “taken together” are useful when they mark real argumentative movement. For example:

> We next extended our Monte Carlo simulations to allow CTL:target interactions in a 1:n ratio...

This transition tells the reader that the analysis is not merely adding another panel; it is relaxing a prior assumption and testing whether the result still holds.

Likewise:

> However, for η = 2, the variance approached the mean and for η = 10 far exceeded the mean...

The contrast matters because it identifies the parameter regime in which the phenomenon emerges.

A paragraph should not merely add observations. It should refine, test, extend, qualify, or explain the central claim.

## State the mechanism or implication

Strong Results prose does not stop at “X increased,” “Y decreased,” or “the groups differed.” It explains what the evidence implies biologically, mechanistically, or conceptually.

For example:

> Multihitting CTL populations initially killed at a low rate, because targets had generally not acquired enough damage to die. Subsequently, targets accumulated damage, and the manifested killing rate per conjugated CTL rose above the killing rate for the single-hit scenario.

This passage connects the observed temporal pattern to a mechanistic explanation. The result is not just that killing rates changed over time; the result is that accumulated sublethal damage can create delayed apparent killing.

Similarly:

> These results demonstrate that missegregation rate can interact with both karyotype fitness landscape shape and ploidy to determine clonal dominance in culture.

This synthesis sentence states the conceptual result of the section. It tells the reader what the evidence changes about the interpretation of karyotype dominance.

The final sentence of a Results section should usually state what the evidence demonstrates, explains, rules out, or motivates next.

## Use caveats locally, not defensively

Caveats should be precise and placed where they affect interpretation. Avoid generic hedging that weakens supported conclusions.

For example:

> Although tumors were rapidly controlled in our model with IFNγ, they were not entirely eradicated.

This caveat is useful because it identifies a specific boundary of the result. It does not weaken the main conclusion that IFNγ-mediated cell-cycle arrest can control tumors in the model.

Similarly:

> Although heterogeneous infiltration may lead to strong spatial variability in killing rate, we conclude that temporal variation in killing is likely large, especially when CTLs cooperate.

This sentence acknowledges an alternative source of variability while preserving the supported inference. The caveat is local; the conclusion remains clear.

Use caveats to define the scope of the result, not to apologize for it.

## Group supporting panels when they share a role

Not every subpanel needs its own sentence. Panels can be grouped when they perform the same argumentative function. For example, multiple panels may together establish that a result is robust across metrics, that a contrast holds across examples, or that an exception is local rather than general.

A grouped citation should still support a specific claim. For example:

> UMAP projections in karyotype space show clear separation between the groups, supporting their interpretation as quasispecies.

The grouped evidence matters because it supports the interpretation of the groups, not because the panels need to be mentioned.

Avoid panel-by-panel prose that reads like a legend. Use panels selectively to build the section’s argument.

## End with synthesis

A Results section should end by consolidating the evidence into the strongest accurate conclusion. Good synthesis sentences often begin with “Taken together,” “These results imply,” “These results demonstrate,” or a similarly direct construction.

For example:

> Taken together, these results imply that multiple-hitting CTLs could explain both heterogeneous and delayed onset killing among clonal CTL populations.

This sentence integrates the section’s observations into a single mechanistic interpretation. It is stronger than ending with another panel description.

A good final sentence should answer: *What has the reader learned that they did not know at the start of the section?*

## Practical drafting checklist

Before finalizing Results prose, check that the section:

* has a conclusion-first heading;
* opens with the question, expectation, or uncertainty being resolved;
* includes only the methodological context needed to interpret the result;
* cites panels where they support specific claims;
* groups panels by argumentative function rather than walking through them mechanically;
* states what each major result means, not merely what changed;
* uses caveats only where they clarify scope or interpretation;
* ends with a synthesis sentence that states the implication of the evidence.

Write the Results as a guided scientific argument. The reader should always understand why the next result appears, what each figure establishes, and how the evidence changes the interpretation.
