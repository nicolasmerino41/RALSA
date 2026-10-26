# Known-link opportunities: amount versus allocation

All ten datasets are retained. Both total-support regressions exclude the target, as do the unchanged persistent-pair parameters. Known regional interaction identities and occurrence-derived opportunities remain supplied.

All-ten descriptive Proportional fit: M=0.000 + 0.556 T_known; slope 55.619%, Pearson r=0.861. These all-ten coefficients do not generate held-out predictions.

All-ten descriptive Intercept fit: M=9.120 + 0.554 T_known; slope 55.402%, Pearson r=0.861. These all-ten coefficients do not generate held-out predictions.

| Pair model input | Mean absolute M error (%) | Mean distribution TV (%) | Mean KS (pp) | Mean gap MAE (pp) |
|---|---:|---:|---:|---:|
| Observed_M | 0.000 | 9.615 | 4.819 | 0.976 |
| Known_Proportional | 28.775 | 14.601 | 11.928 | 3.252 |
| All_Proportional | 58.184 | 20.580 | 18.415 | 5.739 |
| Known_Intercept | 31.426 | 15.067 | 12.599 | 3.453 |
| All_Intercept | 90.848 | 27.820 | 25.912 | 8.410 |

Known-link CO changes mean Pair gap MAE from 5.739 to 3.252 pp, improving 7/10 datasets. Distribution TV changes from 20.580% to 14.601%, improving 7/10. With observed M, the allocation-only reference still has TV=9.615%, despite gap MAE=0.976 pp. A good gap curve therefore does not imply an equally close support-frequency distribution.

The known-link predictor substantially reduces total-support error (mean absolute relative error 28.8%, versus 58.2% for all-pair CO). Salix Galpar and Garraf PP now have almost correct amounts and approach the observed-M reference. Major remaining amount errors occur in Gottin HP (+75.6%), Gottin PP (+98.4%) and Garraf PP2 (-51.1%). Those systems also show large distribution and gap errors. This points to poor amount estimation as an important remaining limitation there. Known-link prediction is not uniformly superior: Nahuel, Olot and Montseny have larger gap errors than under all-pair prediction.

Even with observed M, distribution TV ranges from 4.3% to 15.6%; the support-frequency distribution is not recovered exactly. Some small-frequency tail differences remain visible in CCDFs. TV includes empirical sampling fluctuations, so this is descriptive mismatch rather than a test rejecting the allocation model. With known-link-estimated M, Pair still improves mean gap MAE over Common (3.252 versus 5.056 pp), and mean distribution TV (14.601% versus 19.080%). Allocation heterogeneity therefore retains descriptive value alongside the unresolved amount prediction.

| Dataset | M error: known CO (%) | TV: observed M (%) | TV: known CO (%) | TV: all CO (%) | Gap MAE: known CO (pp) | Gap MAE: all CO (pp) |
|---|---:|---:|---:|---:|---:|---:|
| Quercus | 24.794 | 14.671 | 17.740 | 39.098 | 1.306 | 11.597 |
| Nahuel | -14.779 | 9.665 | 13.077 | 8.275 | 3.576 | 0.919 |
| Salix Galpar | -1.578 | 7.006 | 7.072 | 19.436 | 1.135 | 6.543 |
| Gottin HP | 75.610 | 15.625 | 25.929 | 30.398 | 7.065 | 8.426 |
| Gottin PP | 98.381 | 8.115 | 20.094 | 37.659 | 7.277 | 13.280 |
| Garraf HP | -5.907 | 8.535 | 7.470 | 14.066 | 0.255 | 3.047 |
| Garraf PP | 0.829 | 8.134 | 8.135 | 9.876 | 0.307 | 2.093 |
| Garraf PP2 | -51.056 | 7.381 | 29.470 | 32.310 | 9.787 | 10.906 |
| Olot | 7.739 | 12.681 | 11.486 | 10.803 | 0.471 | 0.315 |
| Montseny | 7.078 | 4.343 | 5.535 | 3.877 | 1.344 | 0.269 |

## Separation

Observed-M predictions isolate the remaining allocation mismatch within this model: their mean support equals observed M/L exactly. Estimated-M predictions add a change in amount while holding pair parameters, opportunities, known links and N fixed. The exact signed decomposition is predicted minus empirical = (observed-M model minus empirical) + (estimated-M model minus observed-M model). It is saved for every frequency bin and removal level. Absolute distances are not additive: the components can reinforce or cancel. An estimated-M curve can occasionally improve through compensation, not necessarily better allocation. This is a model-conditional diagnostic, not a causal or variance decomposition.

TV is half the sum of absolute frequency differences: the fraction of probability mass that must move among K bins. KS is maximum CDF separation; Wasserstein is summed absolute CDF separation in supporting-site units. Singleton error and empirical/model 90th percentiles are also exported. Metrics use expected model distributions, without tests of sampling significance.

Known-link regressions have 0 infeasible rounded predictions. All raw estimates are retained; predictions are projected to [L,sum(C)] only where necessary, as in the earlier all-pair analysis.

OLS with an intercept is retained as a sensitivity check, not selected based on held-out performance. Gap errors average integer removal levels through 80%, including zero. Across-dataset summaries give equal weight to each dataset. Shared sites, collection dependence, interaction-derived occupancy, unresolved corrections, known-link conditioning and unpropagated parameter/regression uncertainty remain; see README. No datasets were removed, and no allocation parameters were refitted.
