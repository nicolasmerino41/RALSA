# Do these co-occurrence predictors improve parameter transfer?

**The specified regressions provide a modest, uneven improvement over the common-parameter baseline.** They lower held-out gap MAE in seven of ten datasets, but worsen it substantially in Gottin HP, Gottin PP and Garraf PP2. They do not provide a consistently accurate dataset-specific prediction.

## Equal-dataset performance

| Measure | Common baseline | Co-occurrence regressions | Own fit, in-sample |
|---|---:|---:|---:|
| Mean MAE (pp) | 3.29 | 2.97 | 0.56 |
| Median MAE (pp) | 2.68 | 2.71 | 0.52 |
| Mean maximum error (pp) | 5.72 | 5.27 | 1.23 |
| Mean signed error (pp) | +0.38 | +0.61 | +0.56 |
| Mean relative error | 29.4% | 25.1% | 4.7% |
| Median relative error | 17.4% | 21.9% | 4.5% |

The mean improvement is 0.31 pp, about 9.5% of baseline mean MAE. Median absolute and relative errors do not improve. Means and medians therefore tell a mixed story, rather than supporting an unqualified success claim. Relative errors are dataset-level MAE divided by mean empirical gap over the same integer levels, not pointwise percentage errors.

## Which datasets benefit?

Negative MAE change means improvement relative to the common baseline.

| Dataset | Regression MAE (pp) | MAE change (pp) | Relative error |
|---|---:|---:|---:|
| Quercus | 0.95 | -1.40 | 6.2% |
| Nahuel | 0.53 | -0.53 | 4.4% |
| Salix–Galpar | 1.89 | -1.11 | 11.9% |
| Gottin HP | 3.61 | +2.12 | 24.9% |
| Gottin PP | 4.11 | +2.27 | 35.6% |
| Garraf HP | 4.87 | -1.59 | 64.0% |
| Garraf PP | 1.78 | -2.42 | 18.8% |
| Garraf PP2 | 7.64 | +1.57 | 40.9% |
| Olot | 0.81 | -0.80 | 5.5% |
| Montseny | 3.52 | -1.24 | 38.9% |

The strongest absolute improvement is Garraf PP. Quercus, Nahuel and Olot reach MAE below one percentage point. Improving the baseline does not necessarily make a prediction accurate: Garraf HP and Montseny retain relative errors of 64% and 39%.

Gottin HP and PP overestimate the gap, with maximum errors of 7.85 and 7.25 pp. Garraf PP2 underestimates it, reaching a maximum error of 14.12 pp at 16/21 sites removed. Its predicted s=0.615 is much higher than fitted s=0.394, while predicted r=0.697 is lower than fitted r=0.842. These two predictor summaries miss its combination of fewer singletons and slower tail decline. These residuals do not identify the ecological cause.

Both regression slopes are positive in every fold, but this does not ensure accurate held-out parameters. Gottin HP requires extrapolation beyond the training range for both predictors; Garraf HP extrapolates for mean log(C), and Garraf PP for the proportion C=1. No alternative specification was selected in response to these failures. Good gap agreement can also involve compensation between parameter errors; it need not imply both parameters are accurately predicted.

## Occupancy checks and limits of validation

C was computed directly from the intersection of species occurrence sets, using the confirmed restored occupancy version. All 4,304 known regional links satisfy 1<=K<=C<=N. No C was replaced by K, and no problematic links were dropped. The confirmed version includes 629 restored species-site presences; 479 uncertain additions remain excluded. Inputs and every link-level quality check are saved.

Held-out K and parameter estimates were not used to fit regressions or construct predictors, apart from the permitted K-based quality check and the assumed known regional identities. However, many occurrence records themselves originated from interaction observations. This is therefore not a fully independent occurrence-based prediction test. It remains conditional on these recorded species presences, known regional links and N.

Gottin HP and PP share 32 sites; the Garraf/Olot/Montseny collection has unresolved within-collection sampling dependencies. Ordinary leave-one-dataset-out results are exploratory, not independent landscape validation. Ten datasets and two estimated response parameters per dataset also provide little basis for a general mechanism claim. Sampling effort, detection and site count may influence opportunity predictors; their effects are not isolated here.

The result is limited to **these two prespecified regressions**: they add some held-out predictive information on average, with appreciable failures and little improvement in typical-dataset error. Failure in several datasets does not establish that co-occurrence generally lacks predictive information. No additional predictors, distribution families or averaging methods were tried. Existing Garraf HP/Olot one-link flags and outstanding uncertain occupancy records remain documented.

Independent checks reproduced all occurrence intersections, ordinary least-squares coefficients and held-out predictions, as well as exact-combination gap curves and performance metrics. See `README.md` for methods and provenance.
