# Preliminary findings — uncorrected occurrence inputs

**A transferable distribution of realisation propensities improves prediction, but neither model adequately recovers all aspects of interaction support.** This is a technical pilot: confirmed presence corrections remain outstanding and can change both n and the pool of zero-interaction pairs. Do not use these numbers as final manuscript results.

## Primary held-out comparison

Each target dataset is predicted from parameters fitted to the other nine. Values below are equal-dataset means; lower is better.

| Measure | M0: shared probability | M1: shared beta distribution |
|---|---:|---:|
| Negative log predictive probability | 1.165 | 0.729 |
| Absolute relative total-support error | 37.4% | 33.5% |
| Absolute relative regional-richness error | 100.2% | 28.6% |
| Regional-retention curve error | 6.31 percentage points | 3.45 percentage points |
| Absolute mean-gap error | 6.33 percentage points | 3.47 percentage points |

M1 improves predictive likelihood in all ten datasets. Its largest improvement is in regional link richness, not total local support. Both models underestimate the mean gap in every held-out dataset. M1's gap underestimates range from approximately 1.1 to 7.7 percentage points; Garraf PP2 is the largest. Thus heterogeneous propensities recover part of the repeatedly supported tail, but substantial mismatch remains.

The grouped holdout preserves this conclusion. M1 improves predictive likelihood in all ten datasets; its mean absolute errors are 38.8% for total support, 33.4% for regional richness and 3.64 percentage points for the gap. These are sensitivity checks with conservative blocks, not ten independent new validation samples.

## Transfer is not the only problem

Within-dataset M1 fits reduce regional-richness error to 4.0% on average, but still underestimate total support in every dataset (25.0% mean absolute error) and the mean gap in every dataset (3.70 percentage points). Optimising likelihood does not guarantee improved prediction of every derived metric: within-dataset gap error is slightly worse on average than LODO.

Consequently, a simple explanation that only the shared parameters fail to transfer is insufficient. Even local fits do not jointly reproduce the number of realised links and their repeated support. Estimated within-dataset mu spans 0.071–0.256 and phi 0.040–2.841, indicating descriptive variation as well. These estimates are not corrected for occurrence ascertainment or observation effort.

## Reading the ecological implications cautiously

The pilot supports the usefulness of modelling heterogeneity, but not the claim that co-occurrence plus two shared parameters completely explains support. It also cannot establish that residual differences are ecological rather than sampling-related. Interaction-derived occupancy conditions the opportunity data on observed interactions; restoring independent presence records changes the problem being fitted.

The falling conditional K/n curves at low n arise even under a constant probability when only K>0 pairs are shown. Therefore that pattern alone does not demonstrate greater specialist consistency. Compare figures 06 and 07 rather than interpreting the conditional curve in isolation.

Start with **03_retention_curves**, then **05_error_diagnostics** to see why a moderately close normalised curve can coexist with substantial errors in absolute support and richness. **02_support_distributions** locates the missing repeated-support tail. The remaining figures show aggregate prediction and the conditional-selection diagnostic.

The next justified step is to correct and audit occurrences, then rerun these same models before proposing additional mechanisms or models. The current analysis does not by itself resolve the earlier sampling-effort concern.
