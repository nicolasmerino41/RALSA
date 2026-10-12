# Does singleton frequency plus one tail parameter suffice?

**Separating singleton frequency improves the description in every dataset, but a single geometric decay still misses repeated-support structure in several systems.** AIC favours the two-parameter model in all ten datasets, with improvements of 5.3–204.3 AIC units. The one-parameter model underestimates singleton frequency throughout. This is evidence for additional statistical flexibility, not for a distinct ecological mechanism generating singletons.

## Fits and adequacy

The table shows the two-parameter estimates. Tail mean is E[K|K>=2]. D is the maximum absolute discrepancy in the full CDF; bootstrap p uses 999 simulations with parameter refitting.

| Dataset | s | r | Tail mean | D | Bootstrap p |
|---|---:|---:|---:|---:|---:|
| Quercus | 0.496 | 0.780 | 5.43 | 0.033 | 0.372 |
| Nahuel | 0.552 | 0.611 | 3.34 | 0.010 | 0.847 |
| Salix–Galpar | 0.455 | 0.764 | 5.23 | 0.040 | 0.001 |
| Gottin HP | 0.491 | 0.780 | 5.53 | 0.075 | 0.012 |
| Gottin PP | 0.584 | 0.731 | 4.71 | 0.051 | 0.001 |
| Garraf HP | 0.708 | 0.629 | 3.69 | 0.038 | 0.139 |
| Garraf PP | 0.637 | 0.720 | 4.58 | 0.083 | 0.001 |
| Garraf PP2 | 0.394 | 0.842 | 6.67 | 0.055 | 0.001 |
| Olot | 0.500 | 0.732 | 4.50 | 0.028 | 0.620 |
| Montseny | 0.670 | 0.690 | 4.19 | 0.033 | 0.001 |

All estimates are interior, and non-singleton sample sizes range from 26 to 586. Singleton fraction, overall mean and non-singleton mean are recovered exactly by the two-parameter MLE here. Their agreement is therefore not an independent adequacy check.

- **Nahuel:** the closest overall fit, with only 1.0 percentage point maximum CDF discrepancy. Both frequencies and CCDF are compatible with this description at the available resolution.
- **Quercus and Olot:** reasonable coarse descriptions (maximum discrepancies 3.3 and 2.8 percentage points). Some local and extreme-tail deviations remain; lack of detected misfit is not proof of a geometric law.
- **Garraf HP:** tentative compatibility, with limited information from just 26 repeated links. The concentration around K=6 and the isolated K=14 observation are not well represented by a smooth geometric decay, despite no clear global bootstrap rejection.
- **Salix–Galpar, Gottin HP, Gottin PP, Garraf PP and Montseny:** the high-support CCDF generally falls more slowly than the fitted geometric tail. A single slope cannot capture the relatively rapid initial decline and the persistent rare high-K links. Garraf PP has the largest full-CDF discrepancy (8.3 percentage points); among its non-singletons the maximum discrepancy reaches 22.9 points. Montseny's full-distribution discrepancy is only 3.3 points, but reaches 10.1 points conditional on repeated links and is detectable with its larger sample.
- **Garraf PP2:** excess mass near low repeated support and near the upper end remains relative to the fitted smooth tail. The finite bound at N=21 is already included, so this residual is not caused by fitting an unbounded geometric distribution.

The conditional non-singleton diagnostic flags the same six datasets. It is not independent corroboration of the whole-distribution test: once singleton mass is fitted exactly, the remaining discrepancy lies entirely within the repeated component. Bootstrap p-values are exploratory, unadjusted, and assume independent links; links sharing species/sites may violate that assumption. Values of 0.001 are the simulation-resolution floor.

## Interpretation

The two parameters provide a useful compact summary, but not an adequate common description of every support distribution. Matching the mean and singleton fraction can coexist with consequential differences in the rare repeated-support tail. The observed/fitted p90 differs by at most one site for this model, yet the logarithmic CCDF still exposes discrepancies beyond that percentile.

Different N, sampling extent, effort and detection can contribute to the observed shapes and their apparent fit. Nothing here identifies their causes. The current validated interactions and recorded site universe are unchanged; outstanding occurrence corrections affect co-occurrence quantities, not these fits unless K or N changes. Existing Garraf HP/Olot one-link quality flags are retained. No removal-gap predictions or further distribution families were fitted.
