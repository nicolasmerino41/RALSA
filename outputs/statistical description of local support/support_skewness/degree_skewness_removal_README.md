# Degree skewness across site removal

Separate consumer and resource figures, ten datasets each. Each species contributes its retained regional partner count once; only species with positive degree are included, matching the earlier degree-distribution analyses. These are degree distributions, not species occupancy or pair-support distributions.

Uniform random site permutations, 2,000 independent sequences per dataset; nested 20/40/60/80% removals within a sequence. Zero removal is calculated once. Site counts are rounded to the nearest integer, matching previous degree analyses. Plot positions label requested percentages; exact realised fractions are recorded in the CSV.

Skewness is the population third central moment divided by variance to the power 3/2, matching support_skewness. Calculate it separately for each retained network, then average; no pooling of degrees across replicates. States with fewer than three connected species or zero variance are undefined and excluded explicitly. Means therefore condition on defined skewness; counts are reported. No uncertainty bands are plotted, as requested. Monte Carlo standard errors are saved.

Confirmed occurrence/interaction inputs are hash-checked; site counts and link supports are verified. Existing sampling, detection, interaction-derived occurrence and shared-campaign limitations remain. No original files are replaced.

| Guild | Datasets with lower mean skewness at 80% than baseline | Undefined draws | Largest Monte Carlo SE |
|---|---:|---:|---:|
| consumer | 9 / 10 | 0 | 0.0146 |
| resource | 6 / 10 | 0 | 0.0133 |
