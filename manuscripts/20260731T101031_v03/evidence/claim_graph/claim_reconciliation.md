# Claim reconciliation

## Outcome

This graph was rebuilt from scratch from 23 validated fresh-context figure
audits. It does not migrate any non-authoritative claim, evidence node, or edge
from the prior graph.

- Worker-proposed observation–claim rows: 103
- Accepted into grouped evidence: 81
- Rejected during reconciliation: 22
- Final evidence nodes: 40
- Reconciliation umbrellas: 5
- Final relationships: 56
- Qualification relationships: none

Every authoritative claim retains its exact locked ID, text, tier, and status.
The locked status strings containing the word `qualified` are preserved as
fields; they do not create qualification edges.

## Umbrella claims

| ID | Reconciled umbrella | Role |
|---|---|---|
| U1 | High ploidy is associated with higher glucose use or uptake and lower cell-number or per-cell yield across complementary readouts. | Groups glucose-use/yield tradeoff evidence and supports C0, C1, C2, and C11. |
| U2 | High ploidy often improves persistence or relative growth under zero or low glucose. | Groups direct and modeled low-resource advantage evidence and supports C0, C3, and C4. |
| U3 | Modeled high-versus-low ploidy fitness changes direction across resource schedules and biological or model contexts. | Groups selection-relevant resource reversals and supports C0, C4, and C11. |
| U4 | Higher ploidy is associated with larger cell area in multiple non-fusion or SUM-159-chem contexts. | Groups positive size contrasts and supports C7. |
| U5 | Cell-area scaling with ploidy reverses or fails in SUM-159-fuse and several other contexts. | Groups genuine size reversals and strongly undermines C7. |

These five nodes are the only non-authoritative claims. No panel-level claim
nodes were created.

## Main adjudications

- F2a glucose drawdown per live-cell AUC is retained as strong evidence within
  U1. Its uniformly positive four-line direction is not cancelled by visual
  near-overlap of selected F2c points.
- Pure-culture persistence/growth readouts were not treated as direct selection
  measurements. Their selection implication is expressed transparently through
  U2→C4.
- F5e and S12 provide the core resource-dependent fitness evidence. F5e is a
  relative-growth surface; S12 begins from equal high/low fractions. Both remain
  explicitly model-derived.
- S12 line/model counterexamples are retained as moderate undermining evidence
  for C4 and C11 rather than averaged away.
- Attenuation toward zero and absence of a visible difference are not
  `undermine` relationships. This removed the false counteredges proposed for
  S8 and S10 sensitivity panels.
- Total-cell gain was not equated with volumetric yield. C12 is linked only to
  the explicit S10 volume-proxy yield evidence, which is mixed across lines.
- Cell-area evidence was divided into positive-scaling (U4) and
  reversed/failed-scaling (U5) umbrellas. The strong SUM-159-fuse reversals are
  retained; C7 therefore remains materially conflicted even though its locked
  status is unchanged.
- S18 cytoplasmic-red evidence was downgraded to weak for C0 because its
  biological meaning is unresolved.

## Authoritative-claim state

| Claim | Incoming relationships | Undermining relationships | Reconciled reading |
|---|---:|---:|---|
| C0 | 16 | 0 | Strongly supported by empirical and modeled response shifts. |
| C1 | 2 | 1 | Strong support from F2/F4/S10, with a weak all-line raw-feature exception. |
| C2 | 1 | 0 | Strong model/posterior support; total-cell gain is only auxiliary. |
| C3 | 4 | 3 | Strong zero-glucose persistence support with genuine line/readout counterevidence. |
| C4 | 3 | 1 | Moderate model-derived support and moderate line/model counterevidence. |
| C7 | 2 | 1 | Mixed: moderate positive scaling evidence and strong repeated reversals. |
| C11 | 4 | 1 | Supported by metabolic cost and high-resource fitness evidence; S12 contains real counterexamples. |
| C12 | 2 | 1 | Weak positive and moderate negative volume-proxy evidence; remains a working hypothesis. |

## Worker-row disposition

`worker_relationship_dispositions.csv` records every one of the
103 worker rows, its final grouped evidence node when accepted,
and the explicit reason when rejected. The canonical
`observation_claim_relationships.csv` contains only final evidence-to-claim
relationships and agrees with the graph.
