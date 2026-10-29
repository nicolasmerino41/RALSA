# Distribution model selection across network states

Consumers; ten communities. At 40% and 80%, cells report median parameters, KS D and KS p for the selected family across 100 removals. Median p is a summary of separate tests, not a combined test.

| Community | 0% (Baseline) | 40% Removal | 80% Removal |
| :--- | :--- | :--- | :--- |
| **C1 — Quercus** | Exponential<br>Param: λ = 0.547<br>KS D: 0.079<br>KS p: 0.165 | Exponential<br>Param: λ = 0.630<br>KS D: 0.086<br>KS p: 0.120 | Exponential<br>Param: λ = 0.874<br>KS D: 0.060<br>KS p: 0.405 |
| **C2 — Nahuel** | Truncated power law<br>Param: α = 1.23, λ = 0.314<br>KS D: 0.024<br>KS p: 0.505 | Truncated power law<br>Param: α = 1.06, λ = 0.399<br>KS D: 0.019<br>KS p: 0.722 | Exponential<br>Param: λ = 1.033<br>KS D: 0.036<br>KS p: 0.532 |
| **C3 — Salix Galpar** | Power law<br>Param: α = 1.45<br>KS D: 0.029<br>KS p: 0.860 | Power law<br>Param: α = 1.46<br>KS D: 0.038<br>KS p: 0.643 | Truncated power law<br>Param: α = 1.27, λ = 0.032<br>KS D: 0.050<br>KS p: 0.445 |
| **C4 — Gottin HP** | Power law<br>Param: α = 1.39<br>KS D: 0.073<br>KS p: 0.790 | Power law<br>Param: α = 1.45<br>KS D: 0.085<br>KS p: 0.637 | Exponential<br>Param: λ = 0.432<br>KS D: 0.128<br>KS p: 0.372 |
| **C5 — Gottin PP** | Truncated power law<br>Param: α = 1.20, λ = 0.040<br>KS D: 0.026<br>KS p: 0.775 | Truncated power law<br>Param: α = 1.22, λ = 0.049<br>KS D: 0.031<br>KS p: 0.653 | Truncated power law<br>Param: α = 1.20, λ = 0.082<br>KS D: 0.040<br>KS p: 0.472 |
| **C6 — Garraf HP** | Exponential<br>Param: λ = 0.557<br>KS D: 0.073<br>KS p: 0.315 | Exponential<br>Param: λ = 0.647<br>KS D: 0.061<br>KS p: 0.545 | Exponential<br>Param: λ = 0.860<br>KS D: 0.065<br>KS p: 0.530 |
| **C7 — Garraf PP** | Truncated power law<br>Param: α = 1.29, λ = 0.271<br>KS D: 0.016<br>KS p: 0.590 | Truncated power law<br>Param: α = 1.11, λ = 0.367<br>KS D: 0.021<br>KS p: 0.427 | Exponential<br>Param: λ = 1.003<br>KS D: 0.039<br>KS p: 0.375 |
| **C8 — Garraf PP2** | Truncated power law<br>Param: α = 0.62, λ = 0.077<br>KS D: 0.021<br>KS p: 0.965 | Truncated power law<br>Param: α = 0.51, λ = 0.101<br>KS D: 0.027<br>KS p: 0.850 | Exponential<br>Param: λ = 0.215<br>KS D: 0.057<br>KS p: 0.465 |
| **C9 — Olot** | Log-normal<br>Param: μ = 0.98, σ = 1.04<br>KS D: 0.109<br>KS p: 0.220 | Log-normal<br>Param: μ = 0.73, σ = 1.03<br>KS D: 0.094<br>KS p: 0.325 | Exponential<br>Param: λ = 0.445<br>KS D: 0.112<br>KS p: 0.445 |
| **C10 — Montseny** | Truncated power law<br>Param: α = 0.77, λ = 0.159<br>KS D: 0.018<br>KS p: 0.865 | Truncated power law<br>Param: α = 0.65, λ = 0.218<br>KS D: 0.023<br>KS p: 0.667 | Truncated power law<br>Param: α = 0.53, λ = 0.347<br>KS D: 0.031<br>KS p: 0.540 |

## What changes?

- At 40% removal, the best-supported family changes in 0/10 communities.
- At 80% removal, the best-supported family changes in 6/10 communities: Nahuel, Salix Galpar, Gottin HP, Garraf PP, Garraf PP2, Olot. In those communities, the selected new family wins a median 76% of individual removals.

At 40% removal the family is unchanged in all ten communities, while its parameters change. In all seven communities initially described by an exponential or truncated power law, the fitted decline/cutoff rate λ increases. This supports a shorter degree scale within the same family at moderate removal. At 80%, label changes become more common, but do not establish a universal regime shift.

A change in the winning label is evidence of a change in relative model support, not by itself a structural regime shift. Overlapping candidate shapes, sampling variation and model uncertainty remain. The overview gives each community equal weight and assigns it to its highest mean AIC-weight family; it does not display that uncertainty.

## How the fits work

1. For each community, calculate the original positive consumer degrees. Generate 100 uniform, fixed-count site removals at each level (nearest integer number of sites); a link survives whenever at least one supporting site remains. Use nested 40%/80% removals within each random ordering. Fit each resulting degree sequence separately, excluding disconnected species.
2. Fit four discrete distributions by maximum likelihood on degrees 1 through U, where U is the initial number of resource species with an observed link. This same biological upper bound is used at all removal levels. No fitted tail threshold and no pooled replication.
3. Unnormalised probabilities: power law k^(-α); truncated power law k^(-α) exp(-λk); log-normal k^(-1) exp(-(log(k)-μ)^2/(2σ^2)); exponential exp(-λk). Each is normalised by summing over k=1,…,U. α and λ are nonnegative. The log-normal is a discrete density-kernel model, matching the previous implementation. A zero-curvature boundary is a power-law limit, not a finite μ/σ estimate.
4. MLE finds the parameters that maximise Σ log P(observed degree). α is estimated from the individual degree observations, not from the slope of a plotted cumulative line. Larger α means faster decline; λ controls exponential decline or cutoff. μ and σ describe position and spread in log-degree space.
5. AIC = 2 × parameter count − 2 × maximised log likelihood. Power law and exponential use one parameter; the other models use two. Within each sequence, weight = exp(−ΔAIC/2) / sum exp(−ΔAIC/2). At baseline select the highest weight; after removal select the family with the highest arithmetic mean weight across the 100 sequences. All communities have equal weight in the overview.
6. KS D is the largest absolute difference between the empirical and fitted cumulative probabilities. To calculate p, simulate 199 degree samples of the same size from the fitted family and refit each sample. p = (1 + number of simulated D values at least as large as the observed D) / 200. This handles discreteness and estimated parameters. It is a model-conditional diagnostic, not adjusted for choosing the winning family. Species degrees are dependent, so it is not a dependence-robust ecological test. A large p does not establish that a family is true.
7. The removal columns show medians of the selected family's separately fitted parameters and KS results. They are not fits to the averaged degree distribution. The full per-sequence diagnostics, weights, p-value ranges, failure fractions and optimisation checks are in the CSVs. No R² is used for fitting or selection.

Inputs are hash-checked against the preceding analysis. Existing site-universe and occurrence limitations are unchanged; see DoesDegreeDistributionChange/outputs/outstanding_data_issues.csv. The same Göttingen sites and other collection dependencies limit independence between communities.