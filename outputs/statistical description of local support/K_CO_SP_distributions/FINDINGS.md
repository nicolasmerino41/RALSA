# Distribution-family results

Fits describe the complete positive distribution on integers 1,...,N. CO uses known interacting pairs, matching the earlier K/CO plots. SP pools consumers and resources with equal species weight, including confirmed species never observed interacting. K is included empirically only, not refitted here.

| Dataset | Quantity | Lowest AIC | AIC | Bootstrap p | Other models within 2 AIC |
|---|---|---|---:|---:|---|
| Quercus | CO | Power_law | 752.6 | 0.118 | Cutoff_power_law, Discrete_lognormal |
| Quercus | SP | Power_law | 477.0 | 0.074 | Cutoff_power_law, Discrete_lognormal |
| Nahuel | CO | Power_law | 614.0 | 0.53 | Cutoff_power_law, Discrete_lognormal, Singleton_geometric |
| Nahuel | SP | Singleton_geometric | 365.2 | 0.404 | Power_law |
| Salix_Galpar | CO | Cutoff_power_law | 5031.8 | 0.174 | Singleton_geometric |
| Salix_Galpar | SP | Cutoff_power_law | 1504.9 | 0.434 | Discrete_lognormal |
| Gottin_HP | CO | Exponential | 700.8 | 0.566 | Singleton_geometric |
| Gottin_HP | SP | Power_law | 496.5 | 0.304 | Cutoff_power_law, Discrete_lognormal |
| Gottin_PP | CO | Cutoff_power_law | 4562.2 | 0.004 |  |
| Gottin_PP | SP | Power_law | 1247.5 | 0.78 |  |
| Garraf_HP | CO | Discrete_lognormal | 339.9 | 0.758 | Cutoff_power_law |
| Garraf_HP | SP | Power_law | 362.0 | 0.784 | Cutoff_power_law, Discrete_lognormal |
| Garraf_PP | CO | Cutoff_power_law | 1414.9 | 0.912 | Discrete_lognormal |
| Garraf_PP | SP | Power_law | 730.5 | 0.602 | Cutoff_power_law, Discrete_lognormal |
| Garraf_PP2 | CO | Power_law | 4576.6 | 0.006 | Cutoff_power_law, Discrete_lognormal |
| Garraf_PP2 | SP | Power_law | 1002.8 | 0.002 | Cutoff_power_law, Discrete_lognormal |
| Olot | CO | Power_law | 454.8 | 0.146 | Cutoff_power_law, Discrete_lognormal |
| Olot | SP | Power_law | 280.1 | 0.376 | Cutoff_power_law, Discrete_lognormal |
| Montseny | CO | Singleton_geometric | 3109.6 | 0.052 | Cutoff_power_law |
| Montseny | SP | Power_law | 1012.3 | 0.444 | Cutoff_power_law, Discrete_lognormal |

AIC ranks alternatives, while refitted parametric-bootstrap KS goodness-of-fit tests assess detectable mismatch. A best-ranked model may still fit poorly. p≥0.05 does not establish the true family; compare curves, AIC differences, boundaries and Monte Carlo errors. AICc is also provided because species samples can be small. AIC values cannot be compared between SP and CO or between datasets with different samples.

499 simulations per fit, each refitted. Candidate definitions match the previous K analysis: decreasing bounded exponential/geometric, decreasing bounded power law, power law with exponential cutoff, discrete lognormal density kernel, and singleton-plus-geometric. All are normalised over 1,...,N. This lognormal evaluates a density kernel at integers rather than integrating rounding bins. A zero-curvature lognormal boundary represents its limiting power-law kernel, not a finite lognormal. There is no tail cutoff selection or unbounded scale-free inference here.

Observations share species/sites and sites have unequal sampling; independent-observation likelihoods and bootstrap tests are working diagnostics, not dependence-robust ecological tests. Shared campaigns further limit independence. Positive K/CO and all-species SP refer to distinct observational units. Similar winning families do not imply equal parameters, common mechanisms or identical shapes.

All source hashes and reconstructed C were verified. Confirmed occurrence restorations and existing data flags remain unchanged: uncertain Gottin PP missing-identity records, Garraf PP2 zero columns and Montseny zero row are excluded as before; many presences are interaction-derived. Full N includes the established site universe but not undocumented sites. Existing analyses are preserved.

Bootstrap convergence flags: 0 of 100 fits have at least one replicate KKT residual ≥1e-6; consult fit_summary.csv before interpreting borderline p-values.