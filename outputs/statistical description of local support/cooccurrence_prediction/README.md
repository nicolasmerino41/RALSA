# Prespecified co-occurrence predictor test

Run from RALSA:

```
julia --startup-file=no --project=. scripts/cooccurrence_prediction.jl
```

Outputs are restricted to `outputs/statistical description of local support/cooccurrence_prediction`. Existing analyses, interaction records and fitted s,r estimates are preserved. This is one exploratory specification; no alternative predictors, regression forms, averaging methods or distribution families are searched.

## Occupancy decision and provenance

Primary occupancy uses the **confirmed** restorations already prepared in `outputs/support_prediction_zero/confirmed/inputs`. Those restore 629 species-site presences: Quercus 156, Gottin HP 308, Garraf HP 115 and Olot 50. Other datasets retain their existing occupancy. The separate inclusive sensitivity's 474 Gottin PP consumer-site records with missing plant identities, four Garraf PP2 plant records and one Montseny plant record remain excluded because their interpretation is uncertain. No uncertain records are promoted based on predictive performance.

The confirmed restoration remains a separate version rather than a change to primary project exports. Its recorded hashes are checked before use. `input_versions.csv` fingerprints occupancy files, the descriptive link and site tables, existing fitted parameters, baseline curves/metrics and correction provenance. `occupancy_audit.csv` records additions. See `support_prediction_zero/presence_record_provenance.csv` and its source checksums for source rows and evidence.

These records are not a fully independent species-occurrence census. Most consumer occupancy and much resource occupancy originally derive from interaction observations. Restoring confirmed non-interacting presences improves the opportunity definition but does not remove that observation dependence. Accordingly this is an exploratory test conditional on these recorded occurrences, known regional interaction identities and N. It is not clean validation against independently measured occupancy throughout.

## C, predictors and quality control

For each known regional consumer–resource link, C is the size of the intersection of its two guild-specific species site sets in the confirmed occupancy table. It is not the number of sites where that pair interacted. The predictor pool includes every known regional link, equally weighted, with no threshold based on its K. Regional identities are explicitly assumed known.

Within each dataset x1=mean(C==1) and x2=mean(log(C)), using natural logarithms. N is the unchanged descriptive site universe. Held-out K is accessed only for the quality-control condition 1<=K<=C<=N and empirical evaluation. The script records every link in `link_cooccurrence_qc.csv`; violations are saved in `unresolved_inconsistencies.csv` and stop the run. It never replaces C by K or drops inconsistent links. A zero-row inconsistencies CSV means no such violations, not that occupancy provenance is complete.

## Regressions and holdouts

Both existing fitted parameters are verified strictly interior. For each target, two ordinary least-squares regressions with intercepts are fitted to the other nine datasets:

logit(s)=a+b*x1; logit(r)=c+d*x2.

Each training dataset is one equally weighted row. Coefficients use only training parameter estimates; target predictors enter inverse-logit predictions. Target s,r and empirical gaps are accessed only after prediction, for evaluation and the in-sample reference. `regressions.csv` records all coefficients, exact training identities, training ranges and whether the target predictor extrapolates. Training design matrices must have full rank.

Predicted s,r define the same singleton-plus-geometric distribution on K=1,...,N, with the tail normalised using target N. The expected random-removal curve uses the existing exact without-replacement formula, evaluated using stable log-products and -expm1. There is no simulation or local recalibration.

The common baseline is the existing nine-dataset arithmetic-mean parameter transfer. Training identities and means are checked against its saved results. Own-fit reconstruction remains an explicitly in-sample reference. Comparison curves and levels are loaded unchanged from the baseline analysis.

## Metrics and figures

All integer removal levels 0,...,N are saved; primary metrics include exactly m/N<=0.8, including zero. Mean absolute error, maximum absolute error and signed mean error are in percentage points. Positive signed error means overestimating the gap. Change in MAE is co-occurrence MAE minus common-baseline MAE; negative means improvement. Relative error is **100 times dataset MAE divided by that dataset's mean empirical gap**, both in percentage points. It is not an average of pointwise ratios. The mean gap is calculated over the same integer levels.

`across_dataset_summary.csv` weights datasets equally, reporting mean and median MAE and relative error, signed bias, mean maximum error, mean MAE change and number improved. It averages dataset-level relative errors, not a ratio pooled across datasets. Different N imply different main-range grids and last included fractions.

- `predictors.csv`: dataset-level x1,x2, N and known link counts.
- `parameter_predictions.csv`: observed fitted, held-out predicted and common-baseline s,r.
- `performance.csv`: all three comparisons and requested error measures.
- `curves.csv`: all integer levels and all four gap curves; gap columns are fractions, errors in pp.
- `01_parameter_predictions.png`: observed versus held-out s,r; circles are regression predictions, crosses the common baseline; dashed diagonal is equality.
- `02_gap_predictions.png`: empirical, regression-predicted, common-baseline and own-fit gaps on common axes.
- `03_performance_comparison.png`: absolute and relative errors, coloured regression circles versus grey baseline diamonds.

## Retained limitations

Gottin HP and PP share the same 32 documented sites. The Garraf HP/PP/PP2, Olot and Montseny datasets share a provider collection; exact campaign and site overlaps remain unresolved. Related datasets remain in training in these ordinary holdouts. See the preserved metadata audit in `parameter_transfer/metadata_relationships.csv` and its source fingerprints. These are not ten independent landscape validations.

Existing interaction-quality flags remain: the fictitious Gottin HP V1 consumer is removed; Garraf HP and Olot retain unresolved one-link discrepancies. N covers recorded sites, not an independent census of wholly unrecorded sampling locations. Sampling effort and detection remain uncorrected. Neither success nor failure identifies mechanisms, estimates wholly undetected interactions, or establishes whether co-occurrence generally has predictive information.
