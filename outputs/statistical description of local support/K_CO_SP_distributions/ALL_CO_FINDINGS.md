# All-pair co-occurrence distribution fits

Every pair with C≥1, including K=0, from confirmed/pairs.csv. Same five finite-support families, maximum-likelihood methods and 499 refitted parametric-bootstrap KS replicates as the known-pair analysis. Original figures and tables are unchanged; 08 adds All CO to the empirical comparison.

| Dataset | Best AIC model | AIC | Bootstrap p | Any candidate p≥0.05 |
|---|---|---:|---:|---|
| Quercus | Cutoff_power_law | 5382.8 | 0.004 | false |
| Nahuel | Singleton_geometric | 2676.7 | 0.368 | true |
| Salix_Galpar | Cutoff_power_law | 12045.2 | 0.002 | false |
| Gottin_HP | Cutoff_power_law | 3859.2 | 0.22 | true |
| Gottin_PP | Cutoff_power_law | 22650.7 | 0.002 | false |
| Garraf_HP | Discrete_lognormal | 1385.0 | 0.98 | true |
| Garraf_PP | Power_law | 3433.4 | 0.138 | true |
| Garraf_PP2 | Cutoff_power_law | 14493.7 | 0.162 | true |
| Olot | Power_law | 1905.6 | 0.52 | true |
| Montseny | Cutoff_power_law | 9485.9 | 0.004 | false |

AIC only ranks models for this particular sample. Do not compare absolute AIC with known-pair AIC, because sample sizes and populations differ. Goodness-of-fit p<0.05 flags mismatch under the independent-pair working model; p≥0.05 does not prove a true family. Dependence through shared species/sites, larger all-pair sample sizes, and sampling affect these diagnostics. Near-threshold p-values have Monte Carlo uncertainty. All current records and original known K/CO counts were verified against confirmed occurrences. Existing occurrence corrections and exclusions are unchanged; see the original README.

K and known-pair CO refer only to observed regional links; All CO includes co-occurring pairs with no observed interaction; SP pools consumers/resources once per species. Each curve has its own denominator. Positive-C pairs only: C=0 combinations are not included. No original fits were rerun.

Bootstrap convergence flags: 0 of 50 fits.