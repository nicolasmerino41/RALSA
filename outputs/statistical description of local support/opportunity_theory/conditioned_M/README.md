# Pair and site allocation conditional on exact total support

This separate test asks whether the existing persistent-pair and shared-site models describe the allocation of local support when its total M is known. It preserves the original analyses and uses only the existing entire-dataset-holdout parameter estimates. No parameters are refitted.

Run from the RALSA root:

```powershell
julia --startup-file=no --project=. scripts/conditioned_M/conditioned_M.jl
```

Default output: `outputs/statistical description of local support/opportunity_theory/conditioned_M`. Optional arguments specify repository and output directories.

## Conditioning and methods

For every target, the known regional interaction identities, site universe, confirmed-restoration opportunity matrix, observed total support M, and positivity of every link are fixed. Individual link support K is not supplied to the generator: only the sum M is used. Empirical K is used for evaluation. The existing a and h were learned from the other nine datasets and remain unchanged.

This conditions the original positive-row models on M; it does not rescale probabilities to match M in expectation, force extra observations into simulated matrices, or retain uniformly sampled latent states after conditioning.

- **Common:** conditioning on exact M makes all eligible positive-row binary allocations equally likely. Its common probability cancels. This is the existing allocation null. Exact conditional expected support frequencies are computed by dynamic programming; its intervals reuse the 5,000 saved null draws.
- **Pair:** before conditioning on M, each known link has an equally likely low/high latent probability logistic(a−h) or logistic(a+h), and its support follows the corresponding zero-truncated binomial. The mixture is formed from the two separately truncated components, preserving the existing model's definition. Conditional support vectors have weights proportional to the product of these mixture probabilities, restricted to sum K=M. Dynamic programming and backward sampling provide exact means and 2,500 independent conditional draws. A common exponential tilt improves numerical stability and cancels on conditioning; it is not a fitted parameter. Conditioning on M automatically changes the implied mixture of latent pair states.
- **Site:** each site initially has an equally likely low/high latent probability, shared by its eligible pairs. The target joint density for allocation A and site states z is proportional to `2^(-N) × product_i [product_j Bernoulli(A_ij; p_j) / (1 − product_j(1−p_j))]`, restricted to opportunity cells, positive row sums and total M. Four Markov chains alternate symmetric occupied–empty cell swaps (preserving M and positive row sums) with site-state flips using the full density, including the row-positivity normalisers. Site states are therefore reweighted consistently with M. Initial allocations are generated without observed individual K. Two chains start with all low/all high states and two with random states. Each chain uses 1,000 burn-in sweeps followed by 4,000 sweeps, saving every fifth: 3,200 draws in total. Each sweep attempts M swaps and N site-state flips. These are Monte Carlo approximations to the exact conditional target, not independent exact samples.

Small-system exhaustive enumeration validates Pair marginals, Common/null equivalence and Site sampling. Split R-hat, autocorrelation-based effective sample size and Monte Carlo standard error are recorded for Site mean gap, singleton fraction and high-state fraction. The screening thresholds are R-hat <1.01 and ESS ≥400; these are diagnostics rather than proof of convergence. The first pass flagged Gottin PP and Montseny. Their final runs use 5,000 burn-in and 40,000 retained-phase sweeps per chain, saving every twentieth (8,000 draws across four chains). This extension is based on sampling diagnostics, not model performance; fitted parameters are unchanged. Any remaining failures are flagged in the findings.

All gaps use the exact expected survival probability under removal of exactly m sites without replacement. Primary means and errors are arithmetic averages over integer levels m/N≤0.8, including zero. No removal sequences are simulated. Pair/Common mean curves use exact conditional support marginals; Site curves average sampled support frequencies. Intervals are central 95% conditional-allocation intervals, not confidence intervals. Common intervals come from existing independent draws; Site intervals use correlated MCMC draws. Comparing against the old unconditioned MAE shows the effect of supplying M but does not represent the same information-restricted prediction task.

## Outputs

- `01_conditional_allocations.png`: empirical mean gaps and singleton percentages against model means and allocation intervals.
- `02_prediction_errors.png`: gap MAE with M fixed and change relative to the same model without fixing M; negative changes indicate smaller error.
- `03_gap_curves.png`: empirical and model expected gap curves.
- `summary.csv`, `across_dataset_summary.csv`: dataset-level and equal-dataset performance.
- `curves.csv`, `support_frequencies.csv`, `draws.csv`: numerical outputs and sampled summaries.
- `sampler_checks.csv`, `site_diagnostics.csv`: algorithm and sampling diagnostics.
- `input_versions.csv`, `parameter_estimates_used.csv`: provenance and unchanged fitted parameters.

## Retained limitations

Confirmed restored occurrences and the same interaction records/site universe are used. Outstanding uncertain occupancy restorations remain excluded (Gottin PP 474 missing-identity records, Garraf PP2 four zero columns, Montseny one zero row). Much occupancy was derived from interaction records. The prior fictitious Gottin HP V1 correction and unresolved one-link differences in Garraf HP/Olot are retained. Occurrence-only corrections affect opportunities, not K unless interactions or the site universe also change. N does not identify wholly undocumented sites.

Gottin HP/PP share sites, and related provider collections have unresolved overlap; the held-out fits are not independent landscape validation. The models concern known regional links and do not predict wholly unobserved interactions. Their two-state heterogeneity is a statistical specification, not an established ecological mechanism. Equal-dataset summaries are descriptive and no combined independence-based significance claim is made.
