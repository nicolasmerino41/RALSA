# Equal-dataset parameter transfer

Run from RALSA:

```
julia --startup-file=no --project=. scripts/parameter_transfer.jl
```

Outputs are confined to `outputs/statistical description of local support/parameter_transfer`. Existing analyses and parameter estimates are preserved.

## Procedure

Each dataset is held out once. Its transferred s and r are the unweighted arithmetic means of the other nine datasets' existing two-parameter estimates. Each training dataset has weight 1/9, irrespective of site or link count. Parameters are not re-estimated, pooled over links, transformed before averaging, or calibrated using the target. `training_parameters.csv` explicitly records all 90 training contributions.

The transferred probability mass is P(K=1)=sbar and P(K=k)=(1-sbar)*rbar^(k-2)/sum(rbar^j,j=0,...,N-2) for k=2,...,N, normalised with the held-out dataset's known N. All existing estimates are finite and identifiable. The target's observed K and own-fit parameters enter only evaluation and the separately labelled in-sample reference.

For every integer m=0,...,N, expected regional retention is the support-probability-weighted mean of 1-C(m,K)/C(N,K). Extinction is zero for K>m. The ratio is evaluated through log-products, with -expm1 used for survival to avoid cancellation. Local retention is 1-m/N; the gap subtracts that common baseline. No simulations are used.

Primary errors are computed over integer levels satisfying m/N<=0.8, including zero. Metrics are arithmetic means or maxima across those levels. No interpolation to exactly 80% is added; the actual final fraction is saved for each dataset. All error columns ending in `_pp` are percentage points. Curve retention/gap columns are fractions. Positive signed error means overestimating the gap.

The own-fit reference is independently recalculated and checked against every integer level of the preceding `gap_reconstruction`, as well as its MAE. Equal-dataset summaries average the ten dataset-level metrics, not their links or removal-grid points. `mean_dataset_maximum_error_pp` averages ten maxima; `largest_dataset_maximum_error_pp` is the worst maximum. Mean signed bias can cancel across datasets, so mean absolute dataset bias is also reported. Deterioration is transferred MAE minus own-fit MAE, in percentage points; negative values mean an incidental improvement on this gap metric, not a better distributional fit.

## Interpretation and related datasets

This tests transfer of a common parameter baseline to held-out datasets, conditional on their known site counts. It does not test prediction from co-occurrence or regional network structure, identify ecological mechanisms, estimate wholly undetected interactions, or establish a best averaging method. Own fits use the target's observed support and are explicitly in-sample.

Original metadata were checked before interpretation:

- Gottin HP and Gottin PP use the same `Gottin/metadata.csv` in `Gottin/data_processing_occurrence.R`. Their validated occupancy tables share all 32 site IDs, matching the 32 metadata sites. In either fold the other interaction type remains in training. These are **not independent landscape holdouts**.
- Garraf HP, Garraf PP, Garraf PP2, Olot and Montseny are processed together in `Garraf-Montseny-Olot/data-processing_occurrence.R`, attributed to Jordi Bosch. This establishes a related data collection, not identical sampling campaigns or sites. Exact spatial/temporal independence, especially among the Garraf datasets, is unresolved in the checked documentation. Shared plot-number strings alone do not establish site identity. Other collection members remain in training in these folds.
- No shared-site relationship was established here for Quercus, Nahuel or Salix–Galpar with the other datasets. This absence of evidence is not an exhaustive independence audit.

`metadata_relationships.csv` and `metadata_source_versions.csv` record the checked relationships, sources and hashes. LODO results must not be described as ten independent landscape validations. No unrequested alternative grouping or averaging scheme is selected using held-out performance.

## Provenance and retained flags

The script uses the ten validated interaction datasets and N values from the descriptive analysis, checking all K values against interaction triples. It verifies the preceding fit/reconstruction input fingerprints and saves `input_versions.csv` and `transfer_input_versions.csv`. Site counts refer to the recorded project universe; wholly unrecorded sampled sites cannot be recovered from the exports.

The fictitious Gottin HP V1 consumer remains removed. Garraf HP and Olot retain unresolved one-link discrepancies with earlier degree tables. Confirmed occurrence additions and uncertain records were considered in the separate zero-presence sensitivity analysis but have not replaced primary occupancy exports. They change co-occurrence-based quantities but not this analysis unless interactions or the site universe change. Sampling effort and detection remain potential influences on fitted and transferred support.

## Outputs

- `summary.csv`: transferred and own-fit s,r; training identities; MAE, maximum absolute error, signed mean error and removal fraction at maximum error for both curves; MAE deterioration.
- `curves.csv`: all integer removal levels through 100%, empirical/transferred/own-fit regional retention and gaps, common local retention, signed errors and main-range flag.
- `across_dataset_summary.csv`: equal-weight summaries for transfer and in-sample reference.
- `support_probabilities.csv`: full target-normalised transferred, own-fit and empirical distributions.
- `01_parameter_transfer_gaps.png`: dataset-coloured empirical curves, dashed black transfer, dotted grey own-fit reference, on common axes through the main range.
- `02_transfer_errors.png`: dataset-level absolute error and signed bias, comparing transfer circles with in-sample diamonds.
- `FINDINGS.md`: accuracy, bias, deterioration and validation limitations.
