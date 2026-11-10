# Degree family transitions under removal

Consumer degrees, same four bounded MLE models and fixed initial resource-count bound as the previous analysis. Same seed and 100 nested removal orders; saved original 0/40/80 AICs agree within 1e-6. Grid: nearest integer site counts to 0:5:95%, with duplicates and complete removal excluded. Actual fractions are used throughout.

Fits require at least three connected species and two distinct degrees; other states are grey and cannot establish a transition. No new KS bootstrap is run: this analysis concerns relative AIC support, not absolute goodness of fit. All fitted parameters, KS D, optimiser residuals and sample sizes are retained.

A sustained departure is the first new family winning three successive evaluated levels. It can later reverse; returns and all switches are reported. Confirmation spans differ by site count. Threshold summaries condition on detected transitions; absence is not assigned 95%. Labels alone do not establish ecological regimes. Ties choose the first listed family; AIC weights expose uncertainty.

| Dataset | Sustained departure /100 | Median onset % | IQR % | Median switches | Later returns /100 |
|---|---:|---:|---|---:|---:|
| Quercus | 55 | 4.5 | 4.5–15.9 | 2.0 | 55 |
| Nahuel | 72 | 50.0 | 37.5–50.0 | 2.0 | 0 |
| Salix_Galpar | 85 | 55.1 | 29.9–70.1 | 3.0 | 40 |
| Gottin_HP | 74 | 62.5 | 50.0–75.0 | 2.0 | 22 |
| Gottin_PP | 11 | 81.2 | 65.6–82.8 | 1.0 | 4 |
| Garraf_HP | 21 | 60.0 | 48.0–76.0 | 2.0 | 18 |
| Garraf_PP | 91 | 65.0 | 55.0–75.0 | 1.0 | 7 |
| Garraf_PP2 | 98 | 71.4 | 61.9–76.2 | 1.0 | 0 |
| Olot | 95 | 42.9 | 21.4–57.1 | 3.0 | 25 |
| Montseny | 65 | 61.1 | 50.0–72.2 | 2.5 | 17 |

Mean weights average valid sequences equally. Raster rows retain replicate order; pale versions of family colours indicate a runner-up within 2 AIC. Species bands show central 95% removal-sequence ranges, not confidence intervals. Small samples can produce unstable selections. Existing sampling, site-universe, occurrence and shared-collection limitations remain. Analysis is conditional on observed networks, not independent landscape validation.

## Interpretation
Mean model support generally changes gradually rather than at a shared sharp threshold. Garraf PP and Garraf PP2 show relatively consistent movement from truncated power law toward exponential; Nahuel also moves in this direction. Salix mainly shifts from power law toward truncated power law. Gottin PP largely retains truncated-power-law support. Quercus and Garraf HP remain predominantly exponential on average despite individual switches. Olot shifts from lognormal toward exponential with substantial uncertainty.

Individual sequences are less smooth than their averages. In Quercus all 55 sequences meeting the three-level departure rule subsequently return to the baseline family; its early onset is not an irreversible transition. Garraf HP similarly returns in 18 of 21 detected departures. Garraf PP and PP2 have 91 and 98 detected departures, with 7 and 0 later returns respectively. Threshold medians describe detected sequences only and should always be read alongside counts, AIC uncertainty and remaining species. The operational three-level definition is a descriptive persistence criterion, not a biological change-point estimator.

Figures: 01 shows mean AIC weights; 02 shows each sequence (pale colours mean runner-up delta AIC < 2; grey means uninformative state); 03 shows remaining connected consumers (mean and central 95% sequence interval); 04 shows onset distributions (grey points), median (black), interquartile range (blue), and number detected out of 100.
