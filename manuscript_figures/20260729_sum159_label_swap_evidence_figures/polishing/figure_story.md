# SUM-159-fuse label-swap analysis: proposed figure narrative

## The short version

Two observations raise the possibility that the SUM-159-fuse labels are
reversed.

First, the cells labelled fuse 2N are often larger than those labelled fuse 4N,
whereas the chem experiment has the expected 4N-greater-than-2N ordering. This
is especially clear between 24 and 48 hours at 0.1, 0.5, and 1 mM glucose. It is
also difficult to reconcile the absolute size of fuse 2N with chem 2N, because
those are putatively the same starting cell line.

Second, cytoplasmic mCherry-channel signal has opposite high/low ordering in the
fuse and chem experiments. Fuse low cells are brighter than fuse high cells;
chem low cells are dimmer than chem high cells. This reversal remains after
comparing cells within the same size strata, so it is not simply a consequence
of fuse low cells being larger.

Neither result is conclusive by itself. Nuclear area is less consistently
swap-like than cell area, and the mixture classifier is an imperfect reference
rather than a true gold standard. Nevertheless, morphology and cytoplasmic
signal are largely independent measurements, and both identify an unexplained
reversal between fuse and chem. Taken together, they make a label or sample
identity problem plausible enough to justify a direct provenance check.

## Reading order

1. **Figure 1** shows all cell- and nuclear-area time courses before any
   analytical subset is selected.
2. **Figure 2** asks whether the differences in the selected interval can be
   explained by cell density or total segmented area.
3. **Figure 3** examines the selected 24–48 h, shared-glucose subset in detail.
4. **Figure 4** asks whether cytoplasmic red supplies independent evidence for
   a fuse-versus-chem reversal.

## Figure 1 — Start with the complete time courses

Figure 1a shows segmented-cell area and Figure 1b nuclear area for all four
contexts and every G0 concentration available in each experiment. Empty panels
represent condition/context combinations that were not acquired.

The broad view gives an initially swap-like impression. In both C00 and I00,
the blue low/2N cell-area curve is above the red high/4N curve through much of
the experiment, particularly at 0.5 and 1 mM (Figure 1a, first two rows). Chem
shows the opposite and more conventional ordering, with high/4N cells larger
than low/2N cells across almost every condition (Figure 1a, third row). The
mixture assay is less consistent: its two classes are closer, cross with time,
and change ordering with glucose (Figure 1a, fourth row).

Nuclear area provides weaker support for a simple swap. Chem again has a clear
high-greater-than-low separation (Figure 1b, third row). I00 generally retains
a smaller high-greater-than-low nuclear difference, while C00 is close to equal
or changes ordering over time (Figure 1b, first two rows). The fuse phenotype
therefore cannot be summarized as “everything about low looks like 4N.”
Instead, the strongest anomaly is cell area, with a less complete nuclear
counterpart.

The gold band marks 24–48 hours. This is the most useful interval for detailed
comparison because it follows the early attachment and seeding transition but
precedes much of the late confluence, starvation-associated deterioration, and
cell loss visible in several trajectories. It is not claimed to remove all
density effects; Figure 2 addresses those directly.

Stars mark 0.1, 0.5, and 1 mM, the only glucose concentrations represented in
all four contexts. Restricting the main comparison to these columns prevents
the apparent fuse/chem difference from being driven by comparing concentrations
that were only acquired in one experiment. The unselected concentrations remain
visible in Figure 1 so that the reader can see that the chosen subset was not
selected solely because it produced the preferred result.

**Conclusion from Figure 1:** some of the complete time-course data look
consistent with a fuse low/high swap, especially the cell-area trajectories,
but the evidence is not uniform across nuclear area, time, glucose, or the
mixture assay. Figures 2 and 3 therefore focus on a prespecified comparable
window rather than claiming that every panel supports the same conclusion.

## Figure 2 — The selected differences are not simply confluence effects

Figure 2 uses only the interval and glucose conditions highlighted in Figure 1:
24–48 h and 0.1, 0.5, and 1 mM. Each point represents one image×ploidy summary.
Panels a–b use the number of segmented cells in the field as the x-axis;
panels c–d use total segmented area. Cell area is shown in panels a and c, and
nuclear area in panels b and d.

Cell and nuclear size do vary with density within several groups. This is
expected and is particularly visible in chem and in some C00 conditions.
However, the low/high differences generally remain at overlapping field counts
and total segmented areas. For example, the C00 and I00 low cell-area points
remain above their high counterparts at 0.5 and 1 mM over the observed density
range (Figure 2a,c, first two rows). Chem retains the opposite ordering
(Figure 2a,c, third row). The nuclear differences are smaller, but chem again
maintains high-greater-than-low values after accounting for the x-axis
(Figure 2b,d).

These plots do not prove that confluence has no effect. In a few panels the
low/high density ranges only partly overlap, and the fitted lines are
descriptive rather than a causal adjustment. They do show that density alone is
not sufficient to explain the lineage-specific reversal: changing field count
or total segmented area does not collapse fuse and chem onto a common high/low
relationship.

**Conclusion from Figure 2:** the 24–48 h area differences are not an artefact
of choosing fields at one particular confluence. This supports using the
highlighted subset for the distributional comparisons in Figure 3.

## Figure 3 — Detailed interpretation of the comparable subset

### Cell-area distributions (Figure 3a)

Figure 3a contains the clearest morphology-based evidence for a swap. At
0.5 and 1 mM, both fuse batches have low/2N cell-area distributions displaced
well to the right of high/4N. The median differences are large:

- C00 low is 44% larger than high at 0.5 mM and 48% larger at 1 mM.
- I00 low is 34% larger than high at 0.5 mM and 41% larger at 1 mM.

At 0.1 mM the same direction is present but weaker. Chem shows the reverse:
low/2N cell area is approximately 28–32% smaller than high/4N across the three
conditions. If fuse and chem are expected to share the same 2N starting
population, the opposing distributions cannot be dismissed as ordinary
differences between unrelated cell lines.

The mixture assay complicates rather than resolves the interpretation.
GFP-defined high cells are larger at 0.1 mM, but the classified low cells are
slightly larger at 0.5 and 1 mM. This could reflect genuine mixture-context
biology, but it could also reflect classification error—particularly for cells
near the GFP decision boundary. The mixture assay is therefore the closest
available identity reference, not a perfect gold standard.

### Nuclear-area distributions (Figure 3b)

Nuclear area is less supportive of a complete swap. Chem and mixture high cells
have larger nuclei across the shared conditions. I00 generally shows the same
direction, although more weakly. C00 is high-greater-than-low at 0.1 mM, nearly
equal at 0.5 mM, and slightly low-greater-than-high at 1 mM.

This matters because a clean label swap would ideally reverse both cell and
nuclear size. The nuclear distributions do not provide that clean reversal.
They weaken the strongest version of the swap claim, while still leaving an
unresolved difference between the nominally identical fuse and chem 2N
populations.

### Nominally identical 2N populations (Figure 3c–d)

Figure 3c–d makes the fuse-versus-chem 2N discrepancy explicit by dividing each
low/2N median by the chem 2N median at the same glucose concentration.

- Fuse 2N cell area is 1.50–1.83× chem 2N (Figure 3c).
- Fuse 2N nuclear area is 1.06–1.47× chem 2N (Figure 3d).
- Classified mixture low cells are also larger than chem 2N and generally lie
  closer to fuse than to chem.

Thus the area data do not uniquely identify fuse as the erroneous experiment.
They demonstrate that the putatively same 2N material is not behaving as the
same population. Possible explanations include a fuse label swap, a chem sample
or annotation problem, passage or clonal drift, selection history, or an
experiment-specific change in morphology. The mixture comparison makes a
fuse-only explanation less secure, but classification uncertainty prevents it
from exonerating fuse.

**Conclusion from Figure 3:** cell area provides substantial evidence
consistent with a swap, but nuclear area and the imperfect mixture reference
make the result suggestive rather than definitive. What is firmly established
is a large and unexplained inconsistency between the nominally identical fuse
and chem 2N populations.

## Figure 4 — Cytoplasmic signal supplies a second, largely independent reversal

The original motivation for measuring cytoplasmic red was to test whether GFP
from the 4N reporter leaked into the mCherry channel. Figure 4a–b does not
support that explanation. Fuse low cells and classified mixture low cells have
more cytoplasmic red than their high counterparts; if GFP bleed-through were
dominant, the GFP-positive mixture high cells should instead be brighter.

The more important observation is the lineage reversal. Both fuse batches have
negative high-minus-low contrasts, whereas chem has a positive contrast
(Figure 4b):

- C00 median high-minus-low contrast: −0.226 robust background units.
- I00: −0.281.
- Mixture: −0.210.
- Chem: +0.187.

This is difficult to reconcile with the statement that fuse and chem contain
the same 2N cells. Whatever physical mechanism produces the extranuclear
signal—protein redistribution, membrane permeability, optical background, or
another property—it distinguishes low from high in opposite directions in the
two lineages.

Cell size was an important alternative explanation because fuse low cells are
larger. Figure 4c therefore repeats the matched comparison after dividing cells
into assay-specific cell-area deciles and giving the shared size strata equal
weight. The reversal becomes, if anything, clearer:

- C00 median size-standardized contrast: −0.326; 0/6 pairs positive.
- I00: −0.334; 0/8 positive.
- Mixture: −0.235; 0/14 positive.
- Chem: +0.378; 7/7 positive.

Size stratification does not establish the biological mechanism, but it removes
the simplest explanation that the red reversal is merely a by-product of the
cell-area reversal. That promotes cytoplasmic red from an interesting
correlate to a substantially more independent piece of evidence for a
fuse-versus-chem identity inconsistency.

Figure 4d shows the contributing image context. It displays larger fields,
every segmented cell outline, and phase, mCherry/red, NIR, GFP where acquired,
and a composite. These images are intended to expose segmentation, density, and
channel-quality differences rather than to serve as quantitative replicates.

The mixture result again needs a classification caveat. A consistent negative
contrast after size stratification supports the fuse direction, but systematic
GFP-classification errors could influence which cells enter each group.

**Conclusion from Figure 4:** simple GFP-channel leakage is unsupported. The
key result is an opposite low/high cytoplasmic-red relationship in fuse and
chem, despite their putatively shared 2N population. Because this reversal
survives size stratification, it is independent evidence consistent with a
label or sample identity problem.

## Overall interpretation

The revised evidence chain is:

1. The full time courses contain a reproducible swap-like cell-area pattern in
   C00 and I00, but not an equally clean nuclear-area reversal (Figure 1).
2. The pattern persists across field cell count and total segmented area, so it
   is not explained solely by confluence (Figure 2).
3. In the comparable 24–48 h subset, fuse and chem have opposite high/low
   cell-area ordering, and their nominally identical 2N populations differ
   markedly in both cell and nuclear size (Figure 3).
4. Fuse and chem also have opposite high/low cytoplasmic-red ordering, and that
   reversal persists after cell-size stratification (Figure 4).

Together, the area and cytoplasmic measurements provide two substantially
independent reasons to suspect a label or sample identity problem. The evidence
is stronger than a morphology-only argument, but it still does not determine
whether fuse was swapped, chem was mislabeled, or the populations diverged
before imaging. The most defensible claim is:

> The putatively identical SUM-159-fuse and SUM-159-chem 2N populations are not
> phenotypically consistent. Both morphology and size-standardized
> cytoplasmic-red signal show lineage-specific reversals that are compatible
> with a label swap. Direct sample provenance or an independent ploidy assay is
> required to identify the affected lineage and resolve the swap hypothesis.
