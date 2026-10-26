# Predicting total support from all-pair co-occurrences

The earlier persistent-pair parameters were learned on nine datasets, but its target total support M was observed. Here both the M regression and the persistent-pair parameters exclude the target. Known regional interaction identities, N and occurrence records are still supplied.

The co-occurrence proxy does not preserve the earlier accuracy. The primary proportional estimate increases equal-dataset Pair MAE from 0.976 to 5.739 percentage points; the intercept regression gives 8.410 pp. The proportional approach remains close for Montseny, Olot and Nahuel, but is substantially worse for Quercus, Gottin HP/PP, Salix Galpar and Garraf PP2. This identifies total-support estimation as a major source of error in this particular prediction pipeline, not a failure of co-occurrence information generally.

Observed all-pair realisation fractions range from 8.3% to 34.2%. The regression overestimates support in several low-realisation systems and underestimates it in Salix Galpar and Garraf PP2. A positive across-dataset association therefore does not make one percentage sufficiently accurate for every target. Pair still improves mean gap error over Common when both receive the proportional M estimate (5.739 versus 7.146 pp), but neither matches the observed-M Pair reference.

All-data Proportional regression: M = 0.000 + 0.188 T; slope 18.845%. This is descriptive only and is not used for held-out curves.

All-data Intercept regression: M = 150.663 + 0.175 T; slope 17.525%. This is descriptive only and is not used for held-out curves.

T sums pair-site co-occurrences for all pairs, including K=0. The existing totalCOVstotalM_acrossDatasets.jl filters to K>0, so its slope is not the requested all-pair percentage. The primary model is proportional (zero intercept); the intercept regression is a prespecified sensitivity check, not selected according to validation performance.

| Support estimate | Allocation | Equal-dataset mean MAE (pp) | Median (pp) | Signed error (pp) | Observed-M mean MAE (pp) | Datasets improved |
|---|---|---:|---:|---:|---:|---:|
| Proportional | Pair | 5.739 | 4.795 | 1.815 | 0.976 | 3/10 |
| Proportional | Common | 7.146 | 5.366 | 3.903 | 2.502 | 1/10 |
| Intercept | Pair | 8.410 | 9.144 | 4.764 | 0.976 | 0/10 |
| Intercept | Common | 8.821 | 9.226 | 5.956 | 2.502 | 1/10 |

| Dataset | Observed M | Proportional M estimate, raw | Pair MAE: observed M | Pair MAE: proportional estimate | Pair MAE: intercept estimate |
|---|---:|---:|---:|---:|---:|
| Quercus | 430 | 992.432 | 0.979 | 11.597 | 11.597 |
| Nahuel | 334 | 375.615 | 1.120 | 0.919 | 7.841 |
| Salix Galpar | 3104 | 1755.978 | 0.985 | 6.543 | 6.461 |
| Gottin HP | 357 | 709.025 | 2.006 | 8.426 | 11.505 |
| Gottin PP | 1987 | 5358.615 | 1.137 | 13.280 | 15.956 |
| Garraf HP | 159 | 203.408 | 0.659 | 3.047 | 9.724 |
| Garraf PP | 747 | 548.979 | 0.282 | 2.093 | 0.727 |
| Garraf PP2 | 4288 | 1858.917 | 0.562 | 10.906 | 11.039 |
| Olot | 253 | 284.847 | 1.417 | 0.315 | 8.564 |
| Montseny | 1451 | 1404.689 | 0.610 | 0.269 | 0.683 |

Primary fold slopes range from 14.809% to 27.532%. No target M or individual K calibrates these slopes. OLS weights datasets equally in its squared-error objective, although high-T datasets have greater leverage. This is not the arithmetic mean of dataset percentages.

Predictions are rounded to integer M, then projected onto the known-link feasible range L≤M≤sum(C). These constraints use known identities and occurrences, not held-out K. There are 5 projected dataset/scenario cases. Raw, rounded and used predictions are all retained.

- Quercus, Proportional: raw 992.432, used 963; Above known-link capacity; projected to sum(C).
- Quercus, Intercept: raw 1129.792, used 963; Above known-link capacity; projected to sum(C).
- Nahuel, Intercept: raw 533.127, used 512; Above known-link capacity; projected to sum(C).
- Gottin_PP, Intercept: raw 6133.058, used 5862; Above known-link capacity; projected to sum(C).
- Garraf_HP, Intercept: raw 384.069, used 269; Above known-link capacity; projected to sum(C).

Mean curves are exact conditional dynamic-programming expectations, with no allocation or removal simulation. Gap errors average all integer removal levels m/N≤0.8, including zero. Aggregate results weight each dataset equally. The local-retention fraction is still 1−m/N; the predicted M changes the support allocation and hence regional retention.

This is a plug-in prediction: regression and parameter uncertainty are not propagated. The primary comparison isolates replacing observed M with estimated M; Common is an allocation control. The site and beta models are not refitted or extended. Known links exclude wholly undetected interactions from the generated network, even though zero-K pairs contribute to the total-co-occurrence predictor.

Shared Gottin sites and related collections prevent interpreting leave-one-dataset-out as independent landscape validation. Much occurrence information is interaction-derived; unresolved occupancy corrections and the existing interaction flags remain. See README and copied metadata. Similar curves do not demonstrate an ecological mechanism.
