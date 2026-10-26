# Fixed-total-support comparison

Existing entire-dataset-holdout parameters were kept fixed. Each model receives the target's opportunity matrix, known regional links, site count, and observed total support M. This evaluates support allocation conditional on M; it is not the earlier prediction task without target support information.

| Model | Mean gap MAE (pp) | Median MAE (pp) | Mean absolute singleton error (pp) | Earlier MAE without fixed M (pp) |
|---|---:|---:|---:|---:|
| Common | 2.502 | 2.495 | 8.844 | 4.561 |
| Pair | 0.976 | 0.982 | 3.421 | 2.881 |
| Site | 2.496 | 2.487 | 8.787 | 3.418 |

Pair improves gap MAE over the fixed-M Common/allocation null in 9 of ten datasets, reducing equal-dataset mean error by 61.0%. It also better describes singleton abundance on average. The Site model produces gaps and singleton fractions close to the allocation null; its transferred site heterogeneity supplies little improvement in this test. This concerns the existing fitted Site specification, not all possible site effects.

| Dataset | Common MAE (pp) | Pair MAE (pp) | Site MAE (pp) |
|---|---:|---:|---:|
| Quercus | 3.325 | 0.979 | 3.330 |
| Nahuel | 1.092 | 1.120 | 1.084 |
| Salix Galpar | 2.418 | 0.985 | 2.408 |
| Gottin HP | 5.500 | 2.006 | 5.452 |
| Gottin PP | 3.275 | 1.137 | 3.207 |
| Garraf HP | 1.431 | 0.659 | 1.451 |
| Garraf PP | 1.554 | 0.282 | 1.599 |
| Garraf PP2 | 1.223 | 0.562 | 1.219 |
| Olot | 2.572 | 1.417 | 2.567 |
| Montseny | 2.634 | 0.610 | 2.648 |

Nahuel is the exception to Pair's improvement over Common. Pair's errors are not uniformly signed across systems: some gaps remain too large and others too small. Its aggregate advantage is therefore not a universal fit or evidence of a uniquely identified mechanism. Most systems use the same transferred Pair parameters (a=0, h=2); Garraf PP2 uses a=−1, h=2. No target-specific parameter adjustment was made.

Fixing M improves Pair relative to its earlier unconditioned curve in 8 datasets. The remaining advantage over Common is evidence that persistent pair differences provide a useful statistical description of allocation beyond opportunity counts and total support. Supplying M is additional target information, so the reduction relative to the earlier exercise should not be described as a like-for-like validation gain.

## Verification and limits

Pair and Common means use exact finite dynamic-programming marginals. Pair intervals use 2,500 independent exact conditional draws; Common intervals reuse 5,000 existing draws. Site intervals and means use four MCMC chains. A small enumerated system checks the samplers, and all 30 dataset/model combinations pass fixed-total, curve and metric consistency checks.

Site diagnostic maximum split R-hat: 1.026; minimum effective sample size: 139.2. Maximum mean-gap Monte Carlo SE: 0.012 pp. Gottin PP and Montseny were extended after the initial diagnostic screen. These diagnostics do not prove convergence.

Remaining sampling flags: Montseny: high_site_fraction, R-hat=1.026, ESS=139.2. Treat those Site summaries as provisional.

Intervals describe conditional allocation variation, not confidence intervals. Gap errors average integer removal levels through 80%, including zero, with equal dataset weighting for aggregate performance. Shared-site/collection dependence, interaction-derived occurrences, excluded uncertain restorations and known-link conditioning remain unchanged; see README and the copied metadata files. No new ecological mechanism is established.
