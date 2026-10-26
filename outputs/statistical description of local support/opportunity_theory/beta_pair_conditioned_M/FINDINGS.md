# Continuous beta versus two-class pair heterogeneity

Both approaches use known regional links, positive support, the same confirmed opportunities and exact observed M. Beta parameters are trained on the other nine datasets using the same four-summary synthetic-likelihood criterion as the saved two-class fit. No target-specific parameter tuning was performed.

| Model | Equal-dataset mean gap MAE (pp) | Median MAE (pp) | Mean absolute singleton error (pp) |
|---|---:|---:|---:|
| Common | 2.502 | 2.495 | 8.844 |
| Two-class pair | 0.976 | 0.982 | 3.421 |
| Beta pair | 0.884 | 0.514 | 4.272 |

Beta improves on Common in 9/10 datasets and on the two-class model in 5/10. These comparisons give each dataset equal weight.

The average gap error is slightly lower for Beta, but the dataset wins are evenly split and its mean absolute singleton error is higher. Both heterogeneity descriptions improve substantially on Common. The defensible conclusion is that the gap improvement survives replacing two classes with continuous pair heterogeneity; this does not establish actual ecological classes or a uniquely correct distribution. Beta's largest gap improvements occur in the two Gottin datasets, which share sites, so they are not independent replications of its advantage.

| Dataset | Two-class gap MAE | Beta gap MAE | Beta mu | Beta concentration |
|---|---:|---:|---:|---:|
| Quercus | 0.979 | 1.927 | 0.500 | 0.5 |
| Nahuel | 1.120 | 1.317 | 0.500 | 0.5 |
| Salix Galpar | 0.985 | 1.483 | 0.500 | 0.5 |
| Gottin HP | 2.006 | 0.517 | 0.500 | 0.5 |
| Gottin PP | 1.137 | 0.249 | 0.500 | 0.5 |
| Garraf HP | 0.659 | 0.510 | 0.500 | 0.5 |
| Garraf PP | 0.282 | 0.373 | 0.500 | 0.5 |
| Garraf PP2 | 0.562 | 0.300 | 0.378 | 2.0 |
| Olot | 1.417 | 1.852 | 0.500 | 0.5 |
| Montseny | 0.610 | 0.312 | 0.500 | 0.5 |

Boundary selections: 0/10. Fits use a prespecified finite grid, not unconstrained maximum likelihood. Parameters minimise training summary discrepancy, not held-out gap errors. The library has 64 replicates per dataset/grid point; simulation and grid resolution can influence close rankings.

Nine folds select mu=0.5 and phi=0.5, corresponding to Beta(0.25,0.25), a continuous U-shaped distribution concentrated toward low and high probabilities. Garraf PP2 selects mu≈0.378 and phi=2. Thus the continuous fit still favours substantial probability heterogeneity; it does not imply that probabilities are concentrated around one common value. These are latent probabilities before the positivity and total-support conditioning.

The beta distribution is placed among known pairs, matching the existing two-class convention. The support distribution integrates a positive-truncated binomial over that beta density. This is a beta mixture of positive binomials, not the conventional beta-binomial truncated after mixing. The distinction is deliberate: changing it would change selection as well as the shape of pair heterogeneity. See README for the equation.

The comparison can show whether the gap result survives replacing two classes with a continuous distribution, and whether one better describes these outputs. It cannot establish two genuine ecological classes, a true beta distribution, or a unique cause of heterogeneity. There are only ten related datasets and the test conditions on observed M.

All input fingerprints, opportunity/support constraints, fixed-M marginals and curve/metric identities were checked. Gaussian quadrature was checked against an analytic example and doubled in resolution for each selected fit. Exact finite enumeration validates the conditional dynamic program. All conditional means use dynamic programming; 2,500 independent draws per dataset supply central 95% allocation intervals, not confidence intervals.

Confirmed-restoration flags, interaction-derived occupancy, excluded uncertain records, shared Gottin sites and unresolved shared-collection overlap remain unchanged. This is exploratory conditional prediction for known links, not independent landscape validation or discovery of unobserved interactions.
