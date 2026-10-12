# Gap reconstruction from the existing support fit

Run from RALSA:

```
julia --startup-file=no --project=. scripts/gap_reconstruction.jl
```

Outputs are confined to `outputs/statistical description of local support/gap_reconstruction`. Existing analyses are preserved. No parameters are refitted, and no new models, predictors or simulations are used.

## Inputs and interpretation

The analysis uses the descriptive `link_support.csv` and `site_universe.csv`, plus the existing `two_parameter_fit/fit_summary.csv` and `fitted_distributions.csv`. Only the singleton-plus-geometric model is used. The saved fitted probabilities are checked against the saved s and r values; checking this identity is not refitting.

Input fingerprints are recorded in `input_versions.csv` and `reconstruction_input_versions.csv`. The script verifies that source exports still match the fit's provenance, reconstructs every K from the validated interaction triples, and checks the site count N against the descriptive site universe. N is the recorded project universe; sampled sites absent from all exports cannot be inferred.

This is **reconstruction from parameters estimated using the observed support distribution**. It is not independent prediction from co-occurrence, regional network structure or new systems. The fitted distribution already matches the empirical singleton fraction and mean support for these interior fits. Good reconstruction would show preservation of this particular summary curve, not general distributional adequacy.

## Exact calculation

For every integer m=0,...,N, a link with support K survives with probability 1-C(m,K)/C(N,K); extinction is zero for K>m. Extinction probabilities are evaluated via the log of the product over j=0,...,K-1 of (m-j)/(N-j). Survival uses -expm1(log extinction) to avoid cancellation near complete removal. No binomial coefficients are evaluated as floating-point factorials.

Empirical regional retention is the empirical frequency-weighted average survival. Fitted regional retention uses the saved fitted probability mass over every K=1,...,N, including supports not observed in that sample. Each regional retention curve is normalised by its own initial distribution of links. Local-support retention is 1-m/N for both. The gap is regional minus local retention, and its signed reconstruction error is fitted minus empirical.

## Main range and metrics

The main range includes **all and only integer levels with m/N<=0.8**, including zero removal. Thus the last included removal fraction can be below 0.8 (e.g. Nahuel ends at 6/8=0.75). There is no interpolation or simulation. All three requested errors are arithmetic averages or maxima over those integer levels, not integrated areas:

- mean_absolute_error_pp = mean(abs(fitted gap - empirical gap))*100;
- maximum_absolute_error_pp = max(abs(fitted gap - empirical gap))*100;
- signed_mean_error_pp = mean(fitted gap - empirical gap)*100.

The levels within each dataset are equally weighted. Different N give different grids and slightly different endpoints; this should be considered when comparing small errors across datasets. No dataset is pooled into another's metric.

## Outputs

- `curves.csv`: every integer removal level from 0 to N, regional and local retentions, empirical and fitted gaps, signed error in percentage points, and the main-range flag. Retentions/gaps are stored as fractions; error_pp is in percentage points.
- `summary.csv`: the three requested errors in percentage points; the actual maximum removal fraction; the level and signed error at maximum absolute discrepancy; the final main-range error; the empirical peak gap within the main range; and the unchanged s,r estimates.
- `01_gap_reconstruction.png`: all datasets on common 0–80% x limits and a common gap scale. Dataset colour denotes empirical support; dashed charcoal denotes reconstruction. Lines connect integer-level expectations only and stop at each actual main-range endpoint.
- `02_gap_errors.png`: signed errors on common symmetric axes. Positive values indicate overestimated regional retention and gap. Points identify the evaluated integer levels.

The script checks probability normalisation, the initial and complete-removal endpoints, monotonic regional retention, the regional-retention lower bound given by local retention, and the identities at one site removed and one site retained.

## Retained quality flags

The current validated interaction records retain the removal of the fictitious Gottin HP V1 consumer. Unresolved one-link discrepancies for Garraf HP and Olot remain documented rather than silently altered. Confirmed host-presence restorations and additional uncertain records were analysed separately in `support_prediction_zero`; they have not replaced the primary occupancy exports. Such corrections change co-occurrence-based quantities but not this calculation unless the interaction records or site universe change. This exercise does not test sensitivity to sampling effort or imperfect detection.
