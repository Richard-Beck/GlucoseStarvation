## Manuscript Methods

### Study design, experimental units, and analysis populations

The study combined longitudinal microscopy, extracellular-glucose
measurements, image segmentation and object classification, empirical
trajectory summaries, and mechanistic modeling to compare baseline- and
elevated-ploidy cell populations across starting-glucose conditions. The
microscopy analysis contained 33,600 acquisition identities linked to recorded
acquisition experiment, culture-well coordinate, within-well field position,
cell line, experimental ploidy label, starting glucose, and acquisition time.
An acquisition experiment denotes the recorded acquisition-batch label, a
physical well denotes one culture-well coordinate within that batch, a field
denotes one imaged position within a physical well, and an acquisition identity
denotes one field at one acquisition time with resolved phase, alive-channel,
and dead-channel inputs. These identifiers are operational: the available
metadata do not establish that acquisition experiments or physical wells are
independent biological experiments.

The data were reduced at different levels for different analyses. Hard alive
and dead object counts, including zeros for fields with no detected objects,
were averaged across field positions within each cell-line, acquisition-
experiment, physical-well, ploidy, starting-glucose, and time group. One
field-mean live count for one physical well-time contributed one count-
likelihood term. The modeled trajectory was defined by cell line, acquisition
experiment, ploidy, and starting glucose and did not distinguish physical wells.
Counts from multiple physical wells and times assigned to the same combination
therefore contributed separate conditional Student-\(t\) terms evaluated
against the same ordinary differential equation (ODE) trajectory. They shared
line-specific kinetic parameters, a ploidy effect, a fitted initial-glucose
group, and a line-by-ploidy initial-cell group. No physical-well random effect,
residual-correlation model, or replicate-level variance term was included.
Recorded dead-cell counts did not enter the mechanistic-model likelihood.

For empirical trajectory features, available field-mean well-time counts were
reduced to condition medians at each hour before trajectory summaries were
calculated. For glucose, every nonmissing entry in source-table columns R1--R4
was treated as one luminescence observation. Those column labels do not
identify whether the entries were repeated reader measurements, aliquots,
culture wells, or independent samples. Multiple glucose entries assigned to
the same trajectory and time contributed separate likelihood terms. A source
row without a ploidy designation was matched to both modeled ploidy
trajectories. Each nonmissing R1--R4 entry from such a row was consequently
copied once to each trajectory and contributed two lognormal or left-censoring
likelihood terms. The two copies were not represented as one shared
observation: there was no common observation identifier, covariance term, or
likelihood weight linking them. Initial glucose was shared across ploidies
within each cell-line, calibration-batch, and nominal-starting-glucose group,
but this shared latent state did not remove the two likelihood contributions.

These reductions do not provide a sampling model for biological runs drawn
from a larger population. Field averaging removes within-well field variation
from the count likelihood, condition medians remove replicate identity from
empirical features, and paired empirical effects do not propagate field-,
well-, or run-level uncertainty. Posterior uncertainty is conditional on the
observed records, mechanistic model, likelihood, and priors. It therefore
supports comparisons within the observed cell-line populations and records,
but does not by itself justify generalization to independently repeated
cultures or broader cell populations.

Five cell-line populations were retained: MCF10A, MDA-MB-231, SNU668,
chemically generated SUM-159, and fusion-derived SUM-159. Their experimental
ploidy labels were MCF10A 2N and 4N, MDA-MB-231 3N and parental, SNU668 low and
high, and both SUM-159 populations 2N and 4N. These were experimental group
labels rather than measurements made in the microscopy manifest.

Analysis scope was defined separately for each analysis family. The five-line
population was the reference population for empirical summaries,
mechanistic-data assembly, and model-structure selection. The five
nondominated structures selected in this population were carried forward
unchanged; structures were not reselected in other populations. A four-line
joint population excluding fusion-derived SUM-159 was an exclusion sensitivity
analysis for parameter and posterior-prediction summaries and was the
population used for the main pooled intervention-support summaries. Five
single-line fits assessed sensitivity to hierarchical information sharing.
Directional transfer was tested for MCF10A from lower to higher ploidy and
SNU668 from higher to lower ploidy in the four-line population, with each
transfer case paired to a matched null design. The seven fitted-data contexts
used for posterior intervention simulations were the five-line joint
population, the four-line exclusion population, and one fit for each of the
five individual lines.

Recorded starting-glucose conditions were 0, 0.1, 0.25, 0.5, 1, 5, and
25 mM, although one fusion-derived SUM-159 plate contained only 0.1, 0.5, 1,
5, and 25 mM. Acquisition times were parsed to hours at minute-level
resolution; mechanistic-model observations were subsequently assigned to an
8-h grid. The available treatment record states that MCF10A cells were seeded
one day before treatment in normal medium, washed on day 0, and placed directly
into glucose-defined medium. The other lines underwent 6 h in starvation
medium, after which an equal-strength 2× glucose medium was added to reach the
target concentration. Plate rows and columns assigned starting glucose and,
where applicable, ploidy.

Cell-line source and authentication, the biological identity and generation of
the elevated-ploidy states, passage-stability criteria, plate format, seeding
density, wells and fields per condition, acquisition duration, and the number
and definition of biological and technical replicates are not established by
the current records and require project-owner confirmation before submission.

Segmentation and alive/dead/junk classification were evaluated in a fixed set
of 90 manually annotated microscopy frames with experimental labels, crop
coordinates, and manual alive-cell coordinates. Classifier-training images
were selected with random seed 42 by sampling 20 images without replacement
from each of the five retained cell-line populations; MCF10A control and HRAS
experiments were excluded from the sampling frame.

The SUM-159 analyses included competition and monoculture imaging. The
competition population comprised 1,200 acquisition identities. Targeted
monoculture measurements used a deterministic selection of complete field time
series described below. A separate multimodal subset contained seven manually
enumerated fields with specified source channels, masks, and crop coordinates.

### Microscopy preprocessing, segmentation, object classification, and repair

The analysis manifest contained one resolved phase, alive-channel, and
dead-channel path for each acquisition identity. Upstream channel loading
nevertheless permitted multiple raw TIFFs to resolve to a channel; when this
occurred, their arrays were summed. Each channel array was independently
scaled to \([0,1]\) using its 1st and 99th intensity percentiles, with values
outside that range clipped. Phase contrast was required. A missing alive or
dead channel was represented by a zero-valued array.

Cell masks were generated with a CellposeModel/CPSAM segmentation model from a
two-channel input containing scaled phase contrast and the sum of the scaled
alive and dead channels. Cellpose channel indices were `[2,1]`, and diameter
was determined automatically. The exact CPSAM weights used at execution time
are not identified by an immutable version or checksum.

Nuclear signal was the scaled sum of the alive and dead channels. Nuclear
segmentation was performed separately within each cell mask. It was skipped
for a cell containing fewer than 10 pixels or having a 95th-minus-5th
percentile nuclear-signal range below 0.03. Otherwise, signal was smoothed with
a Gaussian filter with 1-pixel standard deviation and thresholded at the larger
of the Otsu threshold and the within-cell 60th percentile; the 75th percentile
was used if Otsu thresholding failed. Connected components smaller than
10 pixels were removed, binary closing with a 1-pixel radius and hole filling
were applied, and components smaller than 10 pixels were removed again.

The classifier-training collection contained 100 selected acquisitions and
used matched \(200\times200\)-pixel center crops. The classifier composite
contained dead, alive, and phase channels, whereas the segmentation input
contained phase and summed nuclear signal. Objects were curated interactively
as unlabelled, alive, dead, or junk. Only 57 nonzero curated labels entered
fitting: 28 alive, 21 dead, and 8 junk labels from 23 of the 100 selected
images. The remaining 1,355 persisted object labels were unlabelled and were
excluded from fitting.

For classification, the dead/alive/phase composite was robustly scaled before
two crops were made for each object: a tight bounding-box crop with pixels
outside the object mask set to zero, and an unmasked context crop centered on
the object and extending to approximately three times the bounding-box height
and width, clipped at image boundaries. Each crop was resized to
\(224\times224\) pixels with antialiasing and normalized using ImageNet channel
means \((0.485,0.456,0.406)\) and standard deviations
\((0.229,0.224,0.225)\). A ResNet-18 in evaluation mode, with its final
fully connected layer replaced by the identity mapping, produced a 512-value
embedding for each crop. The two embeddings were concatenated to form
1,024 features; geometry and nuclear measurements were not included. The
backbone used the corresponding torchvision default ImageNet weights and was
not fine-tuned.

The production estimator was a fitted Random Forest with 100 trees, Gini
splitting, unlimited depth, minimum split size 2, minimum leaf size 1,
square-root feature sampling, bootstrap sampling, no class weights, and random
seed 42. All 57 labeled objects were used in one fit. There was no augmentation,
hyperparameter search, class reweighting or resampling, train/validation split,
early stopping, or separate tuning stage. The predicted class was the class
with the greatest mean probability across trees. An exact probability tie was
resolved by the encoded class ordering, selecting the smallest tied class
label. No additional numerical probability threshold was applied.

An area-pass flag denoted a cell-mask area of at least 50 pixels and was stored
separately from the alive/dead/junk prediction. Main per-image alive and dead
counts used hard classifier labels without additionally requiring the area-pass
flag. Fields with no detected objects were retained. Field-level outputs were
combined into image- and object-level measurement tables.

The original processing treated \(A_{\mathrm{nucleus}}\leq A_{\mathrm{cell}}\)
as a strict object-level invariant, with both areas measured in pixels.
Thirty-five fields contained an object with
\(A_{\mathrm{nucleus}}>A_{\mathrm{cell}}\) and were selected for reprocessing.
The full image-status row and all object rows were replaced for each selected
field; rows for unaffected fields were retained. Rerun admission required all
35 requested fields and no others, successful status for every rerun field,
zero rerun errors, and agreement between each field's recorded object count
and emitted object rows. Dataset-level checks required 33,600 unique image
rows matching the acquisition manifest, zero error-status rows, 18,089,147
object rows, consistent schema, order and identifiers, and per-image
object-count agreement.

These checks did not require the replacement measurements to satisfy the
original nuclear-area invariant. Among the 9,706 replacement objects,
35 objects, one in each rerun field, still had
\(A_{\mathrm{nucleus}}>A_{\mathrm{cell}}\), with a maximum nuclear-to-cell area
ratio of 1.0751445. Thus, reprocessing converted the affected fields from
processing failures into structurally accepted measurements but did not
restore the triggering invariant. No numerical comparison with the original
object measurements was required. The affected measurements are retained here
as generated; whether to accept and explicitly document nuclear-mask overshoot
or impose a scientifically justified mask constraint and regenerate them is a
project-owner decision.

### Segmentation detection and alive-classification validation

Manual alive-cell coordinates in each of the 90 validation frames were
translated from crop coordinates to full-image coordinates. Each point was
first assigned to the cell-mask label at
\((\lfloor x\rfloor,\lfloor y\rfloor)\). If that pixel was background, the
nearest labelled pixel within a square radius of 3 pixels was used; otherwise,
the point remained unmatched. Of 1,843 manual alive-cell points, 1,764 mapped
to a mask and 79 did not, giving a point-to-mask recovery proportion of
1,764/1,843 (95.7%). This coordinate-recovery proportion is a detection proxy
for annotated alive cells, not a complete instance-segmentation metric; it
does not assess boundaries, merging, splitting, or unannotated objects.

Mapped points assigned to the same segmented object were deduplicated,
yielding 1,744 unique manual-positive objects. An object touched by at least
one mapped point was manual-positive. All segmented objects whose centroids
fell within the context crop were examined. The classification evaluation set
contained every manual-positive object and every otherwise negative object
whose centroid lay within the annotation target rectangle; other segmented
objects were excluded. This gave 1,744 mapped positives and 783 operational
negatives. The 79 unmatched manual points were reported separately and did not
enter the object-classification false-negative count.

Because only alive cells were manually annotated, negative truth was defined
operationally as a target-rectangle object without a mapped alive point rather
than by independent dead or junk annotations. Hard classifier class 1 was
predicted alive and all other hard classes were predicted negative. With true
positives (TP), false positives (FP), true negatives (TN), and false negatives
(FN) defined against mapped manual-positive status, precision was
\(TP/(TP+FP)\), recall was \(TP/(TP+FN)\), specificity was
\(TN/(TN+FP)\), balanced accuracy was
\((\mathrm{recall}+\mathrm{specificity})/2\), and F1 score was
\(2TP/(2TP+FP+FN)\); a zero denominator produced a missing value.
Receiver-operating-characteristic area under the curve was calculated from
the continuous alive probability using rank-based averaging of ties. These
metrics assess alive versus operational-not-alive classification conditional
on a segmented-object evaluation set, not joint segmentation-and-
classification performance.

None of the 90 validation acquisitions was identical to any of the 100
selected training acquisitions or the 23 acquisitions that supplied fitted
labels, and no training crop was pixel-identical to a validation crop. The
training and validation sets nevertheless shared six source fields when
acquisition time was ignored. One of these fields contributed fitted labels;
in every overlap, the training and validation acquisitions occurred at
different times. The benchmark was therefore held out at the acquisition-frame
level, but not at the broader source-field level. The available records also
cannot exclude unrecorded exploratory use of the validation frames during
informal development.

### Extracellular-glucose processing and mechanistic-data assembly

Extracellular-glucose observations and calibration standards came from four
processed assay batches spanning the modeled line-pair and SUM-159
experiments. Only processed calibration and sample tables are documented. The
assay-kit identity, reader model and settings, sample collection and storage,
reagent preparation, sample and standard volumes, raw luminescence units,
physical units encoded by the source calibration values, and biological
meaning of R1--R4 are not established. Prior descriptions of a Glucose-Glo
assay or Varioskan reader are therefore not asserted without current
confirmation.

For batch \(b\), positive calibration coefficients \(\alpha_b\) and
\(\beta_b\) were estimated by Nelder--Mead optimization:

\[
(\widehat{\alpha}_b,\widehat{\beta}_b)
=\arg\min_{\alpha_b,\beta_b>0}
\sum_i\left[\log L_{bi}-\log(\alpha_bG_{bi}+\beta_b)\right]^2,
\]

where \(G_{bi}\) was the recorded calibration glucose value and \(L_{bi}\)
was luminescence. The batch-specific calibration relation was
\(\mathrm{E}(L\mid G)=\alpha_bG+\beta_b\). The standard deviation of
calibration residuals on the log scale was bounded below by 0.1. A batch with
fewer than two valid calibration points was assigned
\(\alpha_b=10\), \(\beta_b=10\), and log-scale standard deviation 0.5.

Sample luminescence was taken from every available entry in R1--R4. An `LO`
entry was replaced by 20 and treated as left-censored; blank, `NA`, and other
missing entries were removed. Numeric observations were uncensored. A recorded
dilution factor of 1 gave \(d=1\); otherwise
\(d=1000/(\text{recorded dilution factor})\). Expected luminescence was
\(\alpha_bG(t)d+\beta_b\). No upper-censoring or outlier-rejection rule was
identified. Glucose observations were assigned to calibration batch, cell
line, nominal starting glucose, and ploidy when present. The two-term treatment
of ploidy-unspecified entries is described in the unit hierarchy above and is
retained as the implemented statistical rule pending a project-owner decision.

When count input supplied total count \(N\) and dead count \(D\), live count
was \(N-D\); when alive and dead counts were supplied directly, total count was
their sum. Absolute ploidy values used for the continuous covariate were
2N and 4N for both SUM-159 populations and MCF10A, 2.6N and 5.2N for SNU668,
and 3.5N and 4.4N for MDA-MB-231. For ploidy \(p\) and the baseline value
\(p_0\) of that line,

\[
z=\operatorname{round}\{\log_2(p/p_0),3\}.
\]

Observed count times in hours and glucose times recorded as day × 24 were
assigned to index \(\operatorname{round}(\mathrm{hours}/8)+1\) on a grid from
0 h through 8 h beyond the maximum observed count time. No interpolation or
additional off-grid tolerance rule was used. The four-line exclusion,
directional-transfer, matched-null, and single-line populations used the same
definitions, with indices remapped after exclusion. The prior center for each
line-by-ploidy initial live-cell state was the arithmetic mean of its live-cell
observations at time zero.

### Empirical trajectory features and display transformations

For empirical summaries, count observations at each condition and hour were
reduced to medians for live count, dead count, total count
\((\mathrm{live}+\mathrm{dead})\), and dead fraction. Luminescence was
back-transformed as

\[
\widehat{G}=\max\left\{0,\frac{L-\beta_b}{\alpha_bd}\right\}.
\]

If at least one uncensored glucose entry was available, the median of the
uncensored values was used; otherwise, the median derived from censored
measurements was used. These condition medians were descriptive trajectory
points; replicate number was not used as a weight and measurement error was
not propagated.

For any finite, time-ordered trajectory, trapezoidal area under the curve
(AUC) was

\[
\mathrm{AUC}=\sum_i(t_{i+1}-t_i)\frac{y_i+y_{i+1}}{2}
\]

and was missing when fewer than two points were available. Initial, final,
minimum, maximum, net-change, fold-change, and event-time descriptors were
computed from finite ordered observations. Adjacent raw rates were
\((y_{i+1}-y_i)/(t_{i+1}-t_i)\). Smoothed extrema used a smoothing spline with
`spar=0.55` evaluated at 200 equally spaced points; raw differences were used
when fewer than four unique times were available or spline fitting failed.
Other descriptors included live- and dead-cell AUCs, dead-fraction summaries,
glucose initial, final, minimum, AUC and drawdown, maximum glucose-drawdown
rate, and first observed time at or below 0.1 mM. A threshold time was missing
if the threshold was never observed.

Six empirical summaries were defined:

1. Low-glucose growth was computed separately for conditions with starting
   glucose at or below 0.25 mM. A smoothing spline with `spar=0.62` was fitted
   to \(\log(1+\mathrm{live\ cells})\), its first derivative was evaluated at
   260 equally spaced times, the outer 6% of the observed time span was
   excluded, and the 90th-percentile derivative was retained. The
   line-by-ploidy summary was the median across eligible glucose conditions,
   in \(\log(1+\mathrm{cells})\) per hour.
2. High-glucose growth used starting glucose of 5 or 25 mM. Within each
   condition, a linear regression of \(\log(1+\mathrm{live\ cells})\) on hour
   included observations from the start through the first time live count
   reached 70% of its observed maximum. At least three observations and two
   unique times were required; the line-by-ploidy summary was the median slope
   across eligible conditions.
3. Zero-glucose live-cell AUC was the median live-cell AUC at 0 mM over the
   interval from the first through last finite glucose observation, in
   cell-hours.
4. Low-glucose yield intercept was the intercept from a regression of
   live-cell AUC on \(\log(1+G_0)\), including an intercept, across finite
   conditions with starting glucose at or below 1 mM.
5. Glucose use per live-cell AUC was the slope from an intercept-bearing
   regression of glucose drawdown on live-cell AUC over finite conditions with
   starting glucose at or below 1 mM. Drawdown was the first minus last
   calibrated glucose value.
6. Net total-cell yield at 1 mM was the median, across 1-mM conditions, of the
   nonnegative difference between the maximum total-cell count in the glucose
   measurement window and total-cell count at the first observation at or
   after the start of that window.

Feature regressions required at least three finite observations and at least
two distinct predictor values. Exactly two ordered ploidy states were required
for each line. The paired effect was
\(\Delta f=f_{\mathrm{elevated}}-f_{\mathrm{baseline}}\). The per-ploidy
effect was \(\Delta f/\Delta z\), where
\(\Delta z\) was the corresponding difference in continuous log2 ploidy.
These stored empirical effects were not standardized across lines, assigned a
fallback denominator, or clipped. They did not carry bootstrap or
replicate-propagated uncertainty, although regression standard errors and 95%
coefficient intervals were retained for fitted features.

The Figure 2 empirical-effect display applied separate display-only
transformations after excluding fusion-derived SUM-159. For each feature, the
four remaining line-specific per-ploidy effects \(e_i\) were divided by their
across-line sample standard deviation:

\[
e_i^\ast=\frac{e_i}{s^\ast},\qquad
s^\ast=
\begin{cases}
\operatorname{sd}(e_1,\ldots,e_4),&
\operatorname{sd}(e_1,\ldots,e_4)\text{ finite and }>0,\\
1,&\text{otherwise}.
\end{cases}
\]

No across-line mean was subtracted. Displayed line effects were clipped to
\([-3,3]\). The displayed feature-level diamond was the median of the four
unclipped standardized values and was then clipped to the same interval.
These operations changed only the visualization, not the stored paired or
per-ploidy estimators. Supplementary paired-feature and condition-level
feature displays used raw quantities without this standardization, fallback,
or clipping. Live-cell AUC was divided by 1,000 in one supplementary display
only to express its axis in thousands.

### Cell- and nuclear-area summaries and morphology covariates

Cell-mask and nuclear-mask areas were measured in pixels; no physical
pixel-size calibration was available. Descriptive area summaries were
restricted to objects hard-classified as alive. Medians and interquartile
ranges summarized area distributions across the applicable fields, physical
wells, conditions, and times. Interquartile ranges were descriptive
distribution summaries, not confidence or credible intervals. Where reported,
elevated-versus-baseline contrasts were log2 area ratios.

Morphology covariates were constructed in two stages. First, within each
cell-line, acquisition-experiment, physical-well, ploidy, starting-glucose, and
time group, alive objects from all imaged fields were pooled. Cell area was the
median segmented-mask area across all pooled alive objects. Nuclear area was
the median nuclear-mask area among pooled alive objects with detected nuclear
area greater than zero. The 50-pixel area-pass flag was not an eligibility
criterion for these covariates. Only finite group summaries from 24 through
48 h inclusive entered the second stage, and all starting-glucose conditions
were retained.

For each cell line and ploidy state, a scalar area value was the unweighted
median of all finite group medians across acquisition experiments, physical
wells, starting-glucose conditions, and eligible times. Each group median had
equal weight; neither object count nor field count was used as a weight. For
cell or nuclear area \(a\), the model covariate was the within-line log2 ratio
to the mapped baseline state,

\[
x_{\mathrm{area}}=
\log_2\left(\frac{a_{\mathrm{state}}}{a_{\mathrm{baseline}}}\right).
\]

The baseline covariate was therefore zero and the elevated-state covariate was
the realized elevated-versus-baseline log2 area ratio. Cell- and nuclear-area
datasets replaced the continuous ploidy covariate with the corresponding area
covariate for each modeled trajectory.

Nonfinite object values were excluded, and nuclear-area values additionally
had to be greater than zero. Construction failed if a mapped
cell-line/ploidy state had no finite cell- or nuclear-area summary, if its
baseline scalar was nonpositive or nonfinite, or if a modeled trajectory could
not be mapped. There was no positive minimum-cell, minimum-field, or
minimum-replicate threshold beyond requiring at least one finite group summary
per mapped state. Object counts, detected-nucleus counts, and detection
fractions were retained as descriptive diagnostics but did not weight or
filter the covariate. Pooling fields before the group median and then
collapsing experiment, physical well, glucose, and time by an unweighted
median removed the corresponding repeated-measure identities. Measurement
uncertainty was not propagated into GPATH, and the scalar covariate had no
sampling-error interval. Nuclear-mask overshoots described above were retained
without clipping.

### SUM-159 competition and targeted monoculture measurements

Competition images were segmented and classified from phase, green, alive,
and dead channels. The available records identify the relevant fluorescence
channels operationally as green/GFP and red, but do not establish the
biological constructs, fluorophores, promoters, selection reagents, or cell
populations producing either signal.

Green-mixture calibration used all competition wells having at least one
eligible object at the global minimum acquisition time, 0 h; a later
earliest-per-well time was not substituted. Eligible calibration objects were
hard-classified alive, passed the 50-pixel cell-area gate, contained at least
20 cytoplasmic pixels, and had finite locally background-corrected green
signal. Cytoplasm was the cell interior after removing a 1-pixel inner
boundary and excluding the nuclear mask dilated by 2 pixels. Local background
excluded all cell masks dilated by 3 pixels and was estimated by a normalized
Gaussian-weighted smoother with a 75-pixel standard deviation. The mixture
score was
\(\operatorname{asinh}(\mathrm{median\ background\mbox{-}corrected\
cytoplasmic\ green}/1)\).

A two-component Gaussian mixture shared component means and variances across
the 0-h calibration wells but estimated a separate mixture weight for each
well. Expectation-maximization fitting used six starts, comprising three
quantile/partition starts and three seeded random starts, at most
500 iterations, and relative tolerance \(10^{-8}\). A component was rejected
if its effective size was below the larger of 25 objects and 0.5% of
observations; the maximum-likelihood noncollapsed fit was retained. The
decision boundary ignored fitted well-specific mixture weights and was the
unique point between the component means at which their Gaussian densities
were equal. The procedure failed rather than applying an alternative rule if
exactly one central root was unavailable. Values below the boundary were
assigned to the lower-green component and values at or above it to the
higher-green component.

These components were oriented only by the mean of the transformed,
background-corrected green signal. No independent ploidy measurement,
pure-population control, reporter-genotype label, or orthogonal validation
linked lower green to low ploidy or higher green to high ploidy. Competition
groups are therefore described operationally as lower-green and higher-green.
The red measurement is likewise described as background-standardized
red-channel signal rather than assigned an unsupported biological target.

Targeted monoculture candidates were restricted to the two fusion-derived
SUM-159 acquisition contexts and one chemically generated SUM-159 context,
with experimental 2N and 4N rows eligible. Within each context, every observed
ploidy-by-starting-glucose condition was formed and the acquisition times
shared by all such conditions were identified. Five shared times nearest five
equally spaced targets spanning the shared minimum to maximum were selected,
with equal-distance ties resolved toward the earlier time. The chosen times
were 8, 40, 72, 112, and 144 h in one fusion-derived context and 0, 32, 56,
88, and 120 h in the other fusion-derived and chemically generated contexts.

A field series had to contain exactly one acquisition at every selected time.
The target for each condition was two complete series. When at least two
physical wells supplied complete series, sorted wells were chosen at rounded
25th- and 75th-percentile indices. Within the first chosen well, the complete
series nearest field position 2 was selected; within the second, the series
nearest position 3 was selected. Ties were resolved by smaller numeric
position and then lexical series identity. If fewer than two physical wells
supplied complete series, complete series were sorted lexically and the first
two were used. No shortage or missing phase/alive/dead channel group occurred
in the retained selection. The deterministic procedure used no random
sampling or seed and yielded 340 acquisitions: two complete field series for
each observed ploidy-by-glucose condition. Preliminary alive/dead counts were
available as candidate metadata but did not enter the selection ranking after
eligibility filtering.

For the targeted red-fluorescence comparison, retained monoculture fields had
experimental 2N or 4N labels, starting glucose of 0.1, 0.25, 0.5, or 1 mM, and
times from 24 through 48 h inclusive. All targeted fields meeting those
criteria were included. For the competition comparison, two fields per
matched starting-glucose-by-hour stratum were sampled with seed 20260729.
The separate multimodal subset contained seven manually enumerated records:
two from each of the three monoculture contexts and one competition field.
The scientific rationale, representativeness criterion, and crop-selection
rule for those seven records are not recoverable, so they are treated as a
fixed display selection rather than a probability sample.

Targeted monoculture objects were matched one-to-one to objects from the full
classification run by minimizing total squared centroid distance. A match was
retained when centroid distance was at most 5 pixels and relative area
difference

\[
\frac{|A_{\mathrm{target}}-A_{\mathrm{full}}|}
{\max(A_{\mathrm{target}},A_{\mathrm{full}},1)}
\]

was at most 0.25. Monoculture eligibility required a successful match, hard
alive classification, cell area of at least 50 pixels, and targeted nuclear
area greater than zero. Competition eligibility required hard alive
classification, the same area gate, nuclear area greater than zero, and a
lower- or higher-green mixture call.

Red cytoplasm was defined after a 1-pixel erosion of the cell mask and a
2-pixel dilation of the nuclear mask and required at least 20 remaining
pixels. Local background excluded cell masks dilated by 3 pixels. A normalized
Gaussian smoother with a 75-pixel standard deviation was applied on a
fourfold-downsampled image and interpolated to the original resolution.
Background-noise standard deviation was estimated as
\(1.4826\times\mathrm{MAD}\) of noncell residuals, with
\(\mathrm{IQR}/1.349\) and then ordinary standard deviation as fallbacks. The
standardized red signal for one cell was

\[
\frac{\operatorname{median}
(\mathrm{raw\ red}-\mathrm{local\ background})_{\mathrm{cytoplasm}}}
{\mathrm{background\mbox{-}noise\ standard\ deviation}}.
\]

Exact cell area for size-standardized comparisons was the number of pixels
carrying that cell's label in the selected mask. The image-level red summary
was the median standardized per-cell signal. Monoculture images were paired by
experimental context, starting glucose, hour, and field position; lower- and
higher-green cells in a competition image shared the image key. Monoculture
red-signal contrasts were experimental 4N minus 2N differences, whereas
competition contrasts were higher-green minus lower-green differences between
paired image-level medians. Cell- and nuclear-area distributions and
24--48-h trajectories were summarized by context, glucose, time, and
operational group. Area-versus-density, confluence, focused area-distribution,
and same-2N comparisons were supplementary quality-control analyses.

### GPATH mechanistic model

GPATH v1 was an ODE model with live cells \(N\), \(R\) extracellular resource
states \(X_r\), \(P\) pathway fluxes, and \(W\) optional waste states \(W_j\).
Negative state values arising during numerical integration were truncated to
zero when rates were evaluated. Resource uptake was

\[
u_r=\frac{a_{e,r}X_r}{a_{h,r}+X_r},
\]

where \(a_{e,r}\) was maximum uptake and \(a_{h,r}\) was the half-saturation
constant. With resource-yield vector \(Y_R\), allocation matrix \(A\), and
elementwise multiplication \(\odot\), pathway flux was

\[
\boldsymbol{\Phi}=A(Y_R\odot\boldsymbol{u}).
\]

Potential growth was a smooth minimum across pathways,

\[
\mu=-\frac{1}{10}\log\left\{\sum_{p=1}^{P}\exp(-10\Phi_p)\right\}.
\]

For waste-mechanism indicator \(M_j\), with 0 denoting multiplicative
inhibition and 1 denoting additive death,

\[
I=\frac{1}{1+\sum_j(1-M_j)W_j},
\qquad
d=\sum_jM_jW_j.
\]

The net cell rate was \(g=\mu I-m-d\), where \(m\) was the
maintenance/death term. State equations were

\[
\frac{dN}{dt}=gN,\qquad
\frac{dX_r}{dt}=-u_rN,\qquad
\frac{dW_j}{dt}=N\sum_rK_{jr}u_r.
\]

Thus, the waste indicator changed how waste affected net cell rate but not its
production equation. Each column of \(A\) was normalized to sum to one and had
one raw anchor fixed at 1. Non-anchor terms in the first column were
transformed as \(x/(1+x)\) to break symmetry; other non-anchor terms were
positive. In strict-allocation variants, off-diagonal terms were fixed to
zero.

For each positive physical parameter \(\theta_l\), the line- and
ploidy-specific value was

\[
\theta_l(i,z)=
\exp\left[\log c_l+s_l\{\eta_{li}+z\beta_lq_l\}\right],
\]

where \(i\) indexed cell line, \(z\) was the continuous log2-ploidy or
morphology covariate, and \(q_l\) indicated whether the structure enabled a
covariate effect on parameter \(l\). Line effects \(\eta_{li}\) and global
covariate effects \(\beta_l\) had independent standard-normal priors. Prior
centers \(c_l\) and log scales \(s_l\) were \(2.4\times10^{-4}\) and 0.5 for
maximum uptake, 1 and 0.5 for half-saturation, 625 and 0.5 for resource
yields, 1 and 0.5 for allocation terms, 0.01 and 1.5 for waste-production
terms, and 0.05 and 0.2 for maintenance/death.

Initial live cells for line-by-ploidy group \(i\) were

\[
N_{0i}=\exp\{\log(\overline{N}_{0i})+\zeta_i\},
\qquad \zeta_i\sim\mathcal{N}(0,1),
\]

where \(\overline{N}_{0i}\) was the empirical time-zero prior center. Initial
glucose was fitted for each cell-line, calibration-batch, and nominal-glucose
group and shared across ploidies. It had a normal prior centered on nominal
\(G_0\), with standard deviation
\(\sigma_{\mathrm{base}}+\sigma_{\mathrm{rel}}G_0\), truncated between zero and
twice the maximum nominal starting glucose. Scale parameters were
\(\sigma_{\mathrm{base}}=\exp\{\log(0.01)+\epsilon_{\mathrm{base}}\}\) and
\(\sigma_{\mathrm{rel}}=\exp\{\log(0.05)+0.2\epsilon_{\mathrm{rel}}\}\), with
standard-normal priors on both \(\epsilon\) terms. Other resource states began
at 1 and all waste states at zero.

For observed live-cell count \(N_{\mathrm{obs}}\),

\[
\log(N_{\mathrm{obs}}+1)
\sim t_{\nu_N}\{\log(\widehat{N}+1),\sigma_N\},
\]

with the observation-transform Jacobian included in the fitted density.
\(\sigma_N=\exp\{\log(0.2)+0.5\epsilon_N\}\), with
\(\epsilon_N\sim\mathcal{N}(0,1)\). Degrees of freedom were constrained to be
at least 2 and assigned a Gamma prior with shape 2 and rate 0.1.

For an uncensored luminescence observation \(L\),

\[
L\sim\operatorname{LogNormal}
\left[\log\{\alpha_b\widehat{G}(t)d+\beta_b\},\sigma_{L,b}\right],
\]

where calibration coefficients and log-scale standard deviation were fixed at
their batch-specific estimates. A left-censored observation contributed the
corresponding lognormal cumulative probability at its censoring value.

During fitting, the ODEs were integrated with the Stan BDF solver using
relative and absolute tolerances of \(10^{-6}\) and at most 50,000 steps.
Nominal time zero was evaluated at \(10^{-8}\) h. In directional holdout
datasets, only time-zero count and glucose observations from held-out wells
contributed to fitting; their complete trajectories contributed only to
generated holdout log likelihood.

### Candidate model structures, MAP optimization, and model assessment

The count-model candidate universe contained the following 12 structures:

| Resources | Pathways | Waste states | Allocation | Waste effect |
|---:|---:|---:|---|---|
| 1 | 1 | 0 | Full | None |
| 1 | 1 | 1 | Full | Multiplicative inhibition |
| 1 | 1 | 1 | Full | Additive death |
| 2 | 1 | 0 | Full | None |
| 2 | 1 | 1 | Full | Multiplicative inhibition |
| 2 | 1 | 1 | Full | Additive death |
| 2 | 2 | 0 | Full | None |
| 2 | 2 | 0 | Strict diagonal | None |
| 2 | 2 | 1 | Full | Multiplicative inhibition |
| 2 | 2 | 1 | Full | Additive death |
| 2 | 2 | 1 | Strict diagonal | Multiplicative inhibition |
| 2 | 2 | 1 | Strict diagonal | Additive death |

All 12 structures were evaluated in each of 47 data populations, giving
564 MAP specifications, each with 500 starts. The population mapping was:

| Population family | Data populations | MAP specifications | Posterior fits |
|---|---:|---:|---:|
| Five-line joint | 1 | 12 | 5 |
| Leave-one-line-out | 5 | 60 | 5, all for fusion-derived SUM-159 exclusion |
| Single-line | 5 | 60 | 25 |
| Five-line directional transfer | 10 | 120 | 0 |
| Matched five-line null | 10 | 120 | 0 |
| Four-line directional transfer | 8 | 96 | 0 |
| Matched four-line null | 8 | 96 | 0 |
| **Total** | **47** | **564** | **35** |

Each MAP specification used 500 independent CmdStan BFGS starts with the
parameter-transform Jacobian enabled. Without an external initializer,
CmdStan random initialization used radius 2. No optimization seed or BFGS
stopping-control override was set, so those settings depended on the
execution-time CmdStan defaults.

For each start, the optimized log posterior objective \(lp\), return code, and
parameter and generated-quantity values were retained. A start was eligible
whenever \(lp\) was finite, irrespective of return code. The selected start
was the first occurrence of the maximum finite \(lp\), resolving an exact tie
by the smallest start index. Optimization stability was described by the
numbers of finite starts lying strictly within 1 and 5 log-posterior units of
the selected start, and return-code-zero counts were recorded. These were
descriptive only: no minimum near-best-start count, return-code requirement,
gradient criterion, or other failure rule prevented selection of an
unfinished BFGS run. The procedure can therefore select a nonzero-return-code
start and did not apply a prespecified optimizer-convergence gate.

The implemented information criterion was

\[
\mathrm{IC}_{\mathrm{post}}=-2lp_{\max}+2k.
\]

This is a posterior-objective information criterion, not conventional
likelihood AIC, because \(lp_{\max}\) contains prior contributions. The
effective parameter count was

\[
k=L_{\mathrm{eff}}(n_{\mathrm{lines}}+1),
\]

where \(L_{\mathrm{eff}}=3R+WR+1\), with an additional \((P-1)R\)
allocation parameters in unconstrained models. It represented line effects
and one shared ploidy-effect vector but omitted observation and
initial-condition nuisance parameters. Criterion differences were calculated
relative to the minimum within each dataset.

A structure was nondominated within a dataset if no other finite structure had
both no greater \(k\) and no greater posterior-objective deviance
\((-2lp_{\max})\), with at least one strict improvement. Equal-complexity and
equal-deviance ties were not otherwise broken. The exact five finite
nondominated structures selected in the five-line population were:

1. one resource, one pathway, no waste;
2. one resource, one pathway, one multiplicatively inhibitory waste;
3. two resources, one pathway, no waste;
4. two resources, two pathways, full allocation, no waste; and
5. two resources, two pathways, full allocation, one multiplicatively
   inhibitory waste.

These five structures were carried unchanged into posterior inference,
four-line and single-line comparisons, and intervention simulations; they
were not reselected or probability-weighted in those settings.

For a directional holdout case, the transfer fit and matched null had the same
parent population, held-out line, direction, and structure. The transfer score
was

\[
\Delta\ell_{\mathrm{holdout}}=
\ell_{\mathrm{holdout}}^{\mathrm{transfer}}-
\ell_{\mathrm{holdout}}^{\mathrm{null}}.
\]

A positive value favored transfer. No normalization, minimum-effect threshold,
or formal tie rule was applied.

### Morphology-covariate model fits

The morphology analysis fit GPATH v1 to alive-cell area and nuclear area in
24 prespecified dataset-by-model specifications spanning resource, pathway,
waste, allocation-constraint, and maintenance variants. The 24--48-h area
covariates described above were substituted for continuous ploidy and compared
with corresponding continuous-ploidy fits.

Each morphology specification used 500 Jacobian-adjusted BFGS starts.
Initial parameter vectors were rank-matched to count-model solutions within a
radius of 0.2. Per-start results were assessed through finite-start stability,
optimized objective, the posterior-objective information criterion defined
above, and nondominance. Matched morphology- and ploidy-covariate fits were
compared using condition-level log-likelihood and information-criterion
differences. For selected fits, morphology effects were summarized per
covariate unit and at the realized baseline-to-elevated covariate change for
each cell line. The absence of an optimizer acceptance rule also applied to
these fits unless a separate rule is defined and applied.

### Posterior inference, diagnostics, and predictions

Posterior inference was run for 35 count-model fits: the five selected
structures in each of seven fitted-data contexts comprising the five-line
population, the four-line fusion-derived-SUM-159 exclusion, and five
single-line populations. Directional-transfer and matched-null analyses were
MAP-only and were not part of these 35 fits.

Each posterior fit used four independently launched NUTS chains, with
1,000 warmup iterations and 1,500 retained sampling iterations per chain,
target acceptance probability 0.995, maximum tree depth 14, and a dense
Euclidean metric. Warmup was saved by CmdStan, whereas downstream posterior
objects retained post-warmup draws. Chain \(c\) was initialized from the
\(c\)-indexed ranked finite MAP start for the fit, cycling only if fewer
eligible starts were available. Named parameter values were converted to Stan
initial values and constrained values were clamped to valid bounds. Sampling
seeds were unique across the 140 chains and ranged sequentially from 8,200,001
through 8,200,140.

Posterior quality control combined the four chains and calculated divergent
transitions, rank-normalized \(\widehat{R}\), bulk effective sample size, tail
effective sample size, maximum observed tree depth, and energy Bayesian
fraction of missing information. Fit summaries retained total divergences,
maximum \(\widehat{R}\), and minimum bulk and tail effective sample sizes.
No diagnostic acceptance thresholds or prescribed actions for a failed chain
or fit were encoded. The 35 fits therefore cannot be described as having
passed a prespecified posterior-diagnostic gate. The project owner must either
define and apply such a gate or retain this as a limitation on posterior
claims.

When a capped posterior sample was required, draws were allocated as evenly as
possible across chains and selected at evenly spaced post-warmup iteration
indices within each chain. Without a cap, all 1,500 retained iterations per
chain were used. Chain-balanced draws reconstructed line- and ploidy-specific
kinetic parameters for the five-line, four-line exclusion, and single-line
analyses.

Posterior predictions were generated separately for each selected structure
in the five-line and four-line exclusion populations. Each
structure-by-population prediction used 400 post-warmup draws, comprising
100 evenly spaced retained iterations from each of four chains. For each draw
and modeled trajectory, initial live cells and initial glucose were
reconstructed from fitted initial-condition parameters; additional resources
began at 1 and waste states at zero.

State trajectories were evaluated at 8-h output times from 0 through 8 h
beyond the maximum observed count time. They were integrated with LSODA.
Relative and absolute tolerances and maximum internal step were not specified,
so execution-time solver defaults applied; the 8-h times were requested
outputs, not fixed integration steps. Expected glucose luminescence at each
observed glucose index was
\(\alpha_b\widehat{G}(t)d+\beta_b\), with one expected value retained per draw
and glucose observation. Current trajectory displays showed the pointwise
posterior median separately for each structure and trajectory, with no
posterior interval; structures were overlaid without probability weighting or
pooling.

Instantaneous growth surfaces used 100 evenly spaced members of each
400-draw prediction set. Glucose was evaluated at 40 log-spaced values from
\(10^{-4}\) to 25 mM, and the second-resource axis at 40 linearly spaced
values from 0 to 1. Any additional resource was fixed at 1. At each point, the
stored quantity was the instantaneous elevated-minus-baseline net growth rate
for the same line and draw. Displays showed, separately for each structure and
line, the pointwise median of the 100 draw-specific differences and the
fraction of those differences greater than zero. A model-agreement map was the
fraction of the five structure-specific posterior medians that were positive
at each line and grid point; it was not a pooled posterior probability.
Phase-trajectory displays used pointwise posterior-median first-resource versus
second-resource trajectories for the three selected two-resource structures.
No prediction display probability-weighted or pooled structures or
populations. The manuscript-visible prediction summaries used the five-line
collection; the four-line collection was generated with the same protocol but
was not pooled into those displays.

### Glucose-management simulations and support summaries

Posterior intervention simulations crossed the seven fitted-data contexts, the
same five selected structures, three initial seeding densities, and
100 chain-balanced, evenly selected posterior draws per structure and context.
A single-line context simulated its one line; the five-line and four-line
contexts simulated every retained line. Each simulation began with equal
fractions of baseline- and elevated-ploidy cells. For a draw, total seeding was
the mean of its fitted baseline- and elevated-ploidy initial counts multiplied
by 0.5, 1, or 2. Schedule-specific glucose initialized the first resource;
additional resources began at 1 and waste states at zero.

Each of three consecutive 48-h segments was integrated with LSODA. The
0.25-h interval specified the requested output sequence from 0 to 48 h within
each segment; it was not a fixed numerical integration step. Relative and
absolute tolerances, maximum internal step, and other solver controls were not
specified, so execution-time solver defaults governed internal integration.
At a boundary, the last state returned at exactly 48 h was the pre-action
state. The scheduled action was applied to that state, and the resulting
post-action state initialized the next 48-h solve. The boundary state was
represented once in the concatenated 0--144-h trajectory, while pre- and
post-action values were retained for transition summaries.

The schedule grid contained 21 schedules. For each
\(X\in\{0,0.1,0.25,0.5,1,5,25\}\) mM, one schedule refreshed glucose to \(X\)
at both 48- and 96-h boundaries, one carried the state unchanged at both
boundaries, and one added \(X\) mM glucose at both boundaries. Carry changed no
state. Glucose addition added \(X\) to the current glucose state and preserved
all other components. Refresh preserved baseline- and elevated-ploidy cell
numbers and their fraction, set glucose to \(X\), reset each additional
resource to 1, and reset waste states to zero. Zero-glucose addition was
state-identical to zero-glucose carry and was omitted from displayed support
summaries, leaving 20 distinct displayed schedules.

At an endpoint, the numerical ploidy contrast was

\[
\log\left\{
\frac{\max(N_{\mathrm{elevated}},10^{-12})}
{\max(N_{\mathrm{baseline}},10^{-12})}
\right\}.
\]

A draw supported elevated ploidy exactly when
\(N_{\mathrm{elevated}}>N_{\mathrm{baseline}}\); equality did not count as
support. Within each fitted-context, simulated-population, structure, line,
schedule, and day combination, support was the fraction of the 100 draws with
a positive log ratio, and outcome was summarized by the median log ratio.
Context summaries used 1× seeding at days 4 and 6.

For the four-line exclusion population, pooled support gave equal contribution
to each retained line and structure rather than weighting by model
probability. At day 6, each density-by-schedule summary pooled 2,000
evaluations (100 draws × 5 structures × 4 lines); corresponding per-structure
and per-line summaries contained 400 and 500 evaluations. Day-specific
summaries at 1× density used the same 2,000-evaluation pool.

Support frequencies report the numerical frequency of an endpoint advantage
for elevated ploidy. They were not adjusted for model weights and did not
apply a viability constraint, minimum effect size, multiplicity correction,
or formal decision threshold. No formal schedule-ranking analysis,
optimal-schedule selection, or outcome-based tie rule was performed. Schedule
order was fixed by starting glucose and intervention family rather than by
outcome.

### Sensitivity and robustness analyses

The four-line population excluding fusion-derived SUM-159 tested sensitivity
of parameter and prediction summaries to that line's exclusion and provided
the main pooled intervention-support population. Its posterior parameter
summaries were compared with the five-line summaries, with divergences,
\(\widehat{R}\), and effective-sample-size diagnostics retained descriptively.
Independent fits to each of the five lines assessed posterior summaries
without across-line hierarchical information sharing.

Directional-transfer sensitivity used MCF10A lower-to-higher and SNU668
higher-to-lower holdout designs and matched nulls. Matched transfer-minus-null
holdout scores were compared within structure. Intervention support was
summarized separately for five-line, four-line, and single-line fitted-data
contexts. For morphology fits, radius-0.2 rank matching to count-model
solutions was retained as an explicit numerical initialization assumption.

Supplementary morphology quality control compared area with density and
confluence and examined focused area distributions and same-2N comparisons.
Supplementary SUM-159 analyses retained the object, field, and operational
green-component rules stated above.

### Code, data, software, and reproducibility limitations

A reproducible release should archive image- and object-level measurements,
aggregate counts, analysis-ready modeling arrays, analysis-population
definitions, optimization plans, posterior reconstruction specifications,
posterior parameter objects, morphology data and fit definitions,
intervention configurations, and the targeted SUM-159 mask manifest. The
accompanying source should include GPATH equations and the Stan
implementation; image preprocessing, segmentation, classifier training and
inference code; empirical-feature and morphology-summary code; optimization
and posterior-prediction utilities; and intervention-simulation utilities.

Several components needed for computational reproduction remain unpinned. The
execution-time preprocessing and embedding implementation for the
alive/dead/junk classifier, PyTorch and torchvision versions, and external
image-classification implementation are not immutably recorded. The fitted
classifier records scikit-learn 1.7.0, but end-to-end runtime equivalence also
depends on those other components. CPSAM segmentation weights lack an
immutable identifier. The default Cellpose model used for preliminary
targeted-field measurements lacks a recorded model identifier, weights path,
and package version. Exact CmdStan and posterior-analysis versions and the
BFGS defaults used at execution are also unrecorded. Posterior prediction and
intervention integrations relied on execution-time LSODA defaults rather than
explicit tolerances.

The absence of optimizer and posterior-diagnostic acceptance rules is a
methodological limitation rather than solely a software-version limitation.
The implemented two-term treatment of each ploidy-unspecified glucose entry
and the retention of nuclear-mask overshoots likewise require explicit
scientific decisions. No formal intervention schedule ranking was performed.

The final code and data repository, accession or persistent location, license,
computational environment, immutable revision, and member-level identities for
raw image collections remain to be designated. These items must be versioned
in the release or disclosed explicitly; end-to-end computational
reproducibility should not be claimed until they are resolved.
