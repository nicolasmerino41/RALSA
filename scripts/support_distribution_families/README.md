# Statistical families of local interaction support

Separate extension of `scripts/statistical_description_local_support.jl` and `scripts/two_parameter_support_fit.jl`. Original analyses are preserved. K counts distinct sites supporting an observed regional interaction; N counts sampled sites. This is not a species degree-distribution analysis.

## Run

From RALSA:

```text
julia --startup-file=no --project=. scripts/support_distribution_families/fit_distributions.jl
```

Optional arguments: repository root, output directory, and `plot` or `report` to reuse saved fits. Full runs perform 499 bootstrap replicates per comparison. Outputs are PNG figures, CSV tables and FINDINGS.md.

## Models and assessment

All probabilities are discrete and normalized over every integer 1,...,N, including unobserved K values:

- Exponential: proportional to exp(-lambda K), lambda >= 0; equivalent to the earlier bounded geometric.
- Power law: proportional to K^(-alpha), alpha >= 0.
- Cutoff power law: proportional to K^(-alpha) exp(-lambda K), both parameters >= 0. This exponential cutoff is additional to the hard upper bound N shared by every model.
- Discrete lognormal: proportional to K^(-1) exp(-(log K-mu)^2/(2 sigma^2)). This uses the density kernel evaluated at integers, not probabilities integrated over rounding bins. The convex natural parameterization is exp(-a log K-b log(K)^2), b >= 0. Its b=0 boundary is a power-law shape; sigma is then infinite and mu undefined.
- Singleton plus geometric: P(K=1)=s, with P(K=k|K>=2) proportional to r^(k-2), 0<=r<=1. This reproduces the existing two-parameter fit.

Maximum likelihood uses moment equations for one-parameter fits and projected convex optimization for two-parameter fits. Boundary solutions are retained and flagged. AICc uses the number of regional links as sample size and nominal parameter counts. Existing geometric and singleton-geometric likelihoods are checked against saved results. Convergence residuals are exported, including the worst bootstrap residual.

Goodness of fit uses the discrete Kolmogorov-Smirnov distance and a parametric bootstrap with refitting. The reported probability is (1 + exceedances)/(1 + replicates), with Monte Carlo standard error. Low probability indicates mismatch; a high value does not prove the model. AICc only ranks candidates and cannot establish adequacy. Nested boundary solutions and similar AICc values warrant caution when selecting a family.

The optional tail analysis chooses Kmin>=2 by minimizing power-law KS, requiring at least 50 links and five distinct observed support values. Its semiparametric bootstrap resamples the empirical body, generates the fitted tail and reselects Kmin. Samples without an eligible tail are excluded, with valid counts reported. Alternatives use the same selected tail; tail and full-distribution AICc are not comparable. The observed Kmax/Kmin documents the limited range.

These diagnostics treat links as independent observations, although links share species and sites. They are descriptive working-model checks, not dependence-robust tests of ecological universality. Ten datasets are not ten fully independent landscapes. Differences in sampling effort, site count and spatial coverage remain; rescaling K by N cannot remove them automatically. No gap models, co-occurrence predictors or mechanisms are fitted here.

## Provenance and retained flags

Inputs are the descriptive analysis's link_support.csv, site_universe.csv and input_versions.csv, with existing two_parameter_fit/fit_summary.csv used for reproducibility checks. Source hashes are verified; K is checked against distinct-site interaction records and N against the saved site universe. input_versions.csv, source_versions.csv and code_versions.csv record versions and hashes.

The current records derive from the validated 2026-09-29 exports. Confirmed occurrence restorations (Quercus 156, Gottin HP 308, Garraf HP 115, Olot 50) do not change K because interaction records and the site universe are unchanged. Unresolved occurrence issues remain: Gottin PP 474 records missing plant identity, four Garraf PP2 plant zero columns, and one Montseny plant zero row. Much occupancy remains interaction-derived. The fictitious Gottin HP V1 site was previously removed; unresolved one-link differences for Garraf HP/Olot remain flagged. N cannot account for wholly unrecorded sampling sites.

Gottin HP and PP share the same 32 sites. Garraf HP/PP/PP2, Olot and Montseny share provider/collection provenance, with further overlap unresolved. None of these flags is silently repaired in this analysis.

## Reading the outputs

01: frequencies; 02: complementary cumulative probabilities on semi-log axes; 03: the same on log-log axes; 04: relative fit and bootstrap adequacy; 05: tail-only comparisons. Plotting a straight segment does not establish a power law. See FINDINGS.md for dataset-specific conclusions, fit_summary.csv for all parameter estimates and tail_selection.csv for tail eligibility.

Methodological background: [Clauset et al. 2009](https://doi.org/10.1137/070710111) and [Broido & Clauset 2019](https://doi.org/10.1038/s41467-019-08746-5). This implementation adapts distribution comparisons to finite site support.
