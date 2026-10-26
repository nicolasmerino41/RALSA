# Empirical local support versus co-occurrence

Five descriptive figures use the confirmed-restoration occurrence records, with no fitted models or new simulations.

1. `01_absolute_K_including_zeros.png`: K versus n for all pairs with n≥1, including K=0.
2. `02_absolute_K_positive_only.png`: K versus n for known regional interactions, K>0.
3. `03_fraction_including_zeros.png`: K/n versus n, including K=0.
4. `04_fraction_positive_only.png`: K/n versus n for K>0.
5. `05_all_datasets_four_cases.png`: all ten datasets together. Columns include/exclude K=0; rows show K and K/n.

The first four figures have 2×5 dataset panels with linear axes. Their horizontal ranges vary by dataset; K/n shares the same vertical range. The combined figure uses a logarithmic **horizontal** axis to keep small and large co-occurrence counts visible together. Its response axes remain linear and the underlying n values are unchanged.

Points are empirical (n,K) or (n,K/n) combinations. Coincident pairs overlap and are displayed once, so point darkness is not a calibrated measure of frequency. Lines join the arithmetic empirical mean across individual pairs at each exact n, separately for each dataset and inclusion rule; these are not fitted or smoothed relationships. For the fraction plots, all pairs at a given n have the same denominator, so the mean equals mean(K)/n. No pairs are pooled across datasets to calculate a combined trend.

The excluded-zero figures condition on observed regional interaction identities. The included-zero figures contain all recorded consumer–resource pairs that co-occur at least once, including those never observed interacting in the dataset. Pairs with n=0 are not included; their ratio would be undefined. These figures describe recorded interactions and opportunities, not verified ecological absence or the probability that any conceivable species pair interacts.

## Reproduce

```powershell
julia --startup-file=no --project=. scripts/empirical_support_vs_cooccurrence/plot_empirical.jl
```

Optional arguments specify the repository root and output directory. Existing analyses are unchanged.

`empirical_pairs.csv` retains all pair-level values; `plotted_cells.csv` records multiplicity at each plotted point; `empirical_means.csv` contains the exact empirical line values. `dataset_counts_and_checks.csv` gives total, positive and zero-support pair counts. Input fingerprints and the script hash are saved.

Every pair's n and K is verified against distinct sites in the confirmed co-occurrence and interaction records. Positive pairs are additionally matched to the opportunity-theory quality-control table, including site counts. The input manifest is verified before plotting.

## Retained data limitations

Confirmed restored occurrences are included; uncertain restorations remain excluded (Gottin PP 474 missing-identity cases, Garraf PP2 four zero columns and Montseny one zero row). Much occupancy is interaction-derived. The prior Gottin HP V1 correction and unresolved one-link differences in Garraf HP/Olot remain unchanged. K=0 therefore means no recorded interaction, not demonstrated incompatibility. Shared Gottin sites and unresolved provider-collection overlap also remain; dataset curves are descriptive, not independent landscape replicates.
