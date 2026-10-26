# Known-link co-occurrences: amount and allocation of support

Separate extension of the cooccurrence_estimated_M analysis. All ten datasets and existing analyses are retained. No pair, site or beta allocation parameters are refitted.

Run from RALSA:

```text
julia --startup-file=no --project=. scripts/known_link_estimated_M/predict_M.jl
```

Optional arguments: repository root, output directory, and `plot` to regenerate graphics and findings from saved comparisons. Default output is `outputs/statistical description of local support/opportunity_theory/known_link_estimated_M`. The script uses the existing exact conditional-support helper in scripts/beta_pair_conditioned_M/conditional_support.jl.

## Held-out prediction

For each dataset, T_known=sum(C_i) over supplied known regional interactions. C_i is calculated from confirmed species occurrences; individual target K is not used to construct the predictor. Knowing which pairs interact regionally is an explicit assumption and does use regional interaction identities. Unlike the earlier T_all, T_known excludes pairs with no recorded regional interaction.

Fit M=b*T_known by ordinary least squares through the origin on the other nine datasets, one equally weighted dataset row each. Retain M=a+b*T_known as a separate sensitivity check, matching the earlier all-pair exercise; do not select between these specifications using held-out accuracy. Descriptive all-ten fits are not used to predict held-out datasets. Large totals retain greater regression leverage despite equal residual weights.

Round predicted M to an integer and, if necessary, project onto [L,sum(C_i)], retaining raw, rounded and used values with explicit feasibility flags. This matches the previous rule, uses no target K, and does not make an infeasible raw prediction correct.

Use unchanged Common and persistent-Pair parameters from the matching nine-dataset fit. Existing training identities are checked. The Pair model has equal prior low/high probabilities logistic(a-h), logistic(a+h), with separately zero-truncated binomial components. Conditional weights are restricted to positive known-link supports summing to predicted M. Exact dynamic-programming marginals give expected support distributions and random-removal curves; no new simulation is needed. Common is an allocation control. Site and beta models are not extended here.

## Amount versus allocation

Compare the same transferred model under observed M, known-link-estimated M and all-pair-estimated M. Observed-M predictions are an information-rich allocation reference, not fully held-out predictions of amount. They fix mean support to observed M/L while allowing remaining distributional error.

For every K bin and removal level, export the exact signed identity:

```text
estimated-M prediction - empirical
  = (observed-M prediction - empirical)
  + (estimated-M prediction - observed-M prediction)
```

The first term is allocation mismatch conditional on observed amount; the second is the model's response to changing the amount. Pair parameters, known links, opportunities and site counts are held fixed. Changing M necessarily changes the implied distribution, so this is a model-conditional counterfactual comparison, not a unique causal decomposition. Absolute errors and distances do not add, and amount errors can partly cancel allocation errors.

Distribution measures:

- Total variation (TV): half the sum of absolute frequency differences, reported as a percentage of probability mass. This is not the percentage of individually misidentified interactions.
- KS: maximum CDF difference, in percentage points.
- Wasserstein: sum of absolute CDF differences, in units of supporting sites; not normalized by N.
- Mean-support error: (used predicted M - observed M)/L, exactly.
- Signed singleton-frequency error and observed/predicted 90th percentiles.

These compare expected predicted distributions against empirical frequencies. They do not establish statistical rejection or a universal shape family; sampling and regression uncertainty are not propagated. Gap MAE, maximum error and signed error use integer m/N<=0.8, including zero. Relative error divides MAE by mean empirical gap over those levels. Aggregate means and medians weight datasets equally. Full curves through 100% removal are exported.

## Outputs

- 01_support_frequencies.png: empirical frequencies and the three Pair predictions.
- 02_support_ccdf.png: direct comparison of the complementary cumulative distributions on log-log axes.
- 03_gap_curves.png: expected empirical and predicted gaps through 80% removal.
- 04_amount_and_allocation_errors.png: M errors, distribution TV and gap errors by dataset.
- distribution_metrics.csv, comparison_summary.csv: dataset-level and aggregate performance, including both regressions and Common control.
- comparison_distributions.csv, comparison_curves.csv, signed_error_decomposition.csv: complete underlying distributions, curves and signed components.
- regression_folds.csv, full_data_regressions.csv, dataset_totals.csv: held-out and descriptive regressions and totals.
- input_versions.csv, source_versions.csv, comparison_source_versions.csv, code_versions.csv: input hashes and provenance.

Frequency/CCDF panels display through the larger of the empirical maximum and each model's 99.5th percentile (up to N); full support is retained in CSVs and all calculations. CCDFs below half an empirical observation are outside the visible vertical range. Curve interpolation is visual only.

## Verification and retained limitations

The source input hashes, previous all-pair inputs and original conditional-M inputs must match current files. Known-link C/K are cross-checked against the all-pair table. Pair counts, N, opportunities and M must agree across comparisons. Opportunity counts are checked against restored occupancy; the all-pair total is checked against site-level consumer richness times resource richness. Exact conditional distributions must sum to one and have mean used M/L. Reconstructed empirical and observed-M gap curves must reproduce the saved reference. Signed decomposition and its mean-support identity are checked.

Known regional identities, supplied occurrences and site counts remain conditioning information. Gottin HP/PP share sites; related provider collections have unresolved overlap, so folds are not independent landscape validation. Much occupancy derives from interactions. Confirmed restorations are retained; uncertain excluded cases remain Gottin PP's 474 missing-identity records, four Garraf PP2 zero columns and one Montseny zero row. Prior fictitious Gottin HP V1 removal and unresolved one-link differences in Garraf HP/Olot remain flagged. N does not include wholly undocumented sites. Occurrence corrections can change C even if K is unchanged. No datasets are excluded, and no new mechanism or wholly undetected interaction is inferred.
