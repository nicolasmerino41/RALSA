# Co-occurrence-based total support for gap prediction

This separate extension replaces the observed target total support M in the persistent-pair allocation analysis with a leave-one-dataset-out estimate from total co-occurrences T. Original scripts and outputs are unchanged.

Run from RALSA:

```text
julia --startup-file=no --project=. scripts/cooccurrence_estimated_M/predict_M.jl
```

Optional arguments: repository root, output directory, and `plot` to reuse calculated results. Default output is `outputs/statistical description of local support/opportunity_theory/cooccurrence_estimated_M`.

## Design

T sums the number of co-occurring sites for all consumer-resource pairs, including pairs with no observed interaction (K=0). It is not the number of unique co-occurring pairs and not the opportunity total limited to known links. The all-pair table is verified against confirmed occurrence records, including its total against the sum of consumer richness times resource richness across sites. Positive-pair K/C, site counts and M are checked against the earlier conditioned-M analysis. Hashes verify the original opportunity-theory input versions.

For each target, fit M=bT on the other nine datasets by ordinary least squares through the origin. As an explicitly separate sensitivity analysis, fit M=a+bT by ordinary least squares with an intercept, also excluding the target. Each dataset supplies one row with equal residual weight. Large T can nevertheless give a dataset greater leverage. These slopes are not averages of individual realisation fractions. All-ten regressions are reported only descriptively; their coefficients never generate held-out curves. The existing script totalCOVstotalM_acrossDatasets.jl currently filters to K>0, whereas this analysis includes all K=0 pairs.

Use the existing leave-one-dataset-out persistent-pair a,h parameters, unchanged. Their recorded training identities must equal those of the corresponding M regression. The original parameters were chosen from the original simulation-based summary-score grid, not newly estimated here or selected against this task's results. Predicting M does not change those parameters.

The known regional links remain fixed and positive. Round predicted M to an integer. If it is outside the feasible interval [L,sum(C)] over known links, project to the nearest endpoint and explicitly retain raw/rounded/used estimates and flags. This is a feasibility rule based only on known identities and occurrence opportunities; it is not calibration against target K. Such cases are failed unconstrained M predictions, even though a feasible curve is shown. At either endpoint the support vector is deterministic.

For interior M, use the unchanged exact conditional-support helper from scripts/beta_pair_conditioned_M/conditional_support.jl. Persistent-pair support weights are the existing equal mixture of two separately zero-truncated binomials with probabilities logistic(a-h) and logistic(a+h), conditioned on sum(K)=estimated M. Common is a control with no pair heterogeneity; its probability cancels when conditioning on M. Exact dynamic-programming marginals give expected curves without allocation sampling. The site model is not rerun: the comparison focuses on the requested persistent-pair result, with Common as a control.

For removal of exactly m sites, extinction probability is zero for m<K, otherwise product over j=0,...,K-1 of (m-j)/(N-j). Expected regional retention averages survival over the predicted support distribution. Subtract 1-m/N to obtain the gap. The local-support fraction remains 1-m/N whether M is observed or estimated. Primary errors average integer m/N<=0.8, including zero; all integer levels through 100% are saved. Aggregate means and medians give datasets equal weight. MAEs and signed errors are in percentage points; relative error is 100*MAE/mean empirical gap, not pointwise ratios.

The original Pair observed-M curves are a conditional reference, not a prediction that lacked target support information. Both the new regressions and pair parameters exclude the target. No target K is supplied beyond the explicitly assumed known identities; K supplies training M for other datasets, input quality checks and empirical evaluation. This remains conditional on known regional links and supplied occurrence data. Zero-K pairs enter the T predictor but are not allowed to become new links in the allocated network.

## Outputs and verification

- 01_support_estimation.png: descriptive all-ten regressions and held-out raw M estimates against observations.
- 02_prediction_errors.png: Pair gap errors with observed versus estimated M.
- 03_gap_curves.png: empirical, observed-M Pair, estimated-M Pair and proportional-estimate Common curves.
- dataset_totals.csv, full_data_regressions.csv, regression_folds.csv: totals, slopes and predictions.
- summary.csv, across_dataset_summary.csv, curves.csv, support_frequencies.csv: performance and underlying values.
- result_checks.csv: occurrence completeness, compatibility, exact-M and probability checks.
- input_versions.csv, source_versions.csv, code_versions.csv, parameter_estimates_used.csv: provenance.

## Retained limits

These are plug-in estimates without propagated regression or pair-parameter uncertainty. Ten datasets provide limited transfer assessment. Gottin HP/PP share sites; related collections have unresolved overlap, so folds are not independent landscape validation. Many occurrence records derive from interactions. Confirmed restorations are retained. Outstanding uncertain occupancy issues remain excluded: Gottin PP 474 missing-identity records, four Garraf PP2 zero columns, and one Montseny zero row. Prior fictitious Gottin HP V1 removal and unresolved one-link differences in Garraf HP/Olot remain. N does not represent wholly undocumented sites. Corrections to occurrence records can change T and C even when they do not change K. This test does not estimate wholly undetected interactions or establish ecological mechanisms.
