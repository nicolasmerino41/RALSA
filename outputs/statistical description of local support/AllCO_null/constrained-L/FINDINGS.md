# All CO null with constrained L

All positive co-occurring pairs are eligible. Observed M and L are fixed; realised pair identities vary. The target is uniform over opportunity-cell subsets conditional on both M and L, giving support vectors weights proportional to product binomial(C,K). This differs from the original known-link null, where pair identities are fixed.

Four heat-bath Markov chains, 500 burn-in sweeps and 1,250 retained allocations per chain, separated by five sweeps. Each sweep attempts P pair-block updates. These are correlated MCMC samples, not independent exact draws. A finite-state enumeration checks detailed balance; split Rhat, ESS and Monte Carlo standard errors assess mixing for singleton percentage, gap, selected opportunity counts and overlap with observed links. Interpret intervals only alongside diagnostics.csv.

| Dataset | Singleton % empirical / null [95%] | Gap pp empirical / null [95%] |
|---|---:|---:|
| Quercus | 49.62 / 19.99 [15.04, 25.56] | 15.40 / 23.00 [21.92, 24.02] |
| Nahuel | 55.21 / 42.11 [37.42, 47.24] | 12.18 / 14.74 [13.99, 15.46] |
| Salix_Galpar | 45.47 / 29.83 [27.69, 31.84] | 15.88 / 20.21 [19.78, 20.64] |
| Gottin_HP | 49.07 / 26.12 [20.37, 32.41] | 14.51 / 21.47 [20.09, 22.74] |
| Gottin_PP | 58.39 / 40.66 [38.28, 43.02] | 11.54 / 16.24 [15.77, 16.68] |
| Garraf_HP | 70.79 / 62.55 [56.18, 68.54] | 7.62 / 9.37 [8.26, 10.47] |
| Garraf_PP | 63.69 / 56.85 [53.54, 60.00] | 9.45 / 11.76 [11.11, 12.41] |
| Garraf_PP2 | 39.40 / 14.76 [12.93, 16.55] | 18.66 / 26.42 [26.01, 26.83] |
| Olot | 50.00 / 29.51 [22.83, 35.87] | 14.89 / 19.76 [18.53, 20.93] |
| Montseny | 67.04 / 53.45 [51.20, 55.73] | 9.04 / 12.23 [11.83, 12.64] |

singleton_pct: empirical below / inside / above central 95% allocation interval: 0 / 0 / 10 datasets.

gap_pp: empirical below / inside / above central 95% allocation interval: 10 / 0 / 0 datasets.

Diagnostics: maximum split Rhat = 1.0012732447053552; minimum ESS = 2924.6179650745808; failing Rhat<1.01 and ESS≥1000: 0 / 40.

Link richness, zero-pair percentage and mean positive support M/L match observations by construction. Remaining differences concern allocation and identities, not total support or link richness. Gap is the arithmetic mean of exact expected regional retention minus 1-m/N over integer m/N≤0.8, including zero. Singleton percentages and distribution curves use positive-support links. Intervals describe conditional allocations, not confidence intervals or removal-sequence variation.

The same confirmed inputs are hash-checked and all C and K counts reconstructed as in free-L. Sampling/detection limitations, interaction-derived occurrence records and shared collections remain; see copied occurrence_audit.csv and metadata_relationships.csv. No mechanism is established by rejecting equal allocation. Existing analyses are preserved.

Run with Julia using the RALSA project and this run.jl; append plot to redraw saved outputs.
