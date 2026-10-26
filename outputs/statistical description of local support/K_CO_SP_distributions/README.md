# K_CO_SP_distributions

Run from RALSA:

```text
julia --startup-file=no --project=. scripts/K_CO_SP_distributions/run_analysis.jl
```

Optional arguments: repository root, output directory, and `plot` to regenerate plots/report from saved fits. Outputs default to `outputs/statistical description of local support/K_CO_SP_distributions`. All five figures and tables are saved together; existing analyses are preserved.

## Populations

- K: support of observed regional interactions, K≥1. Included in the joint empirical figure, not refitted here.
- CO: number of sites where both partners of each known regional interaction occur, C≥1. Pairs never observed interacting are excluded, matching the earlier CO plots.
- SP: number of distinct occupied sites per recorded species, pooling consumer/resource identities with equal species weights. Includes confirmed occurrence-only species. It is not degree-weighted.

All empirical curves are complementary cumulative distributions: proportion with at least n sites. The shared figure compares raw counts, not mean-normalised values. Lines connect occupied discrete values for display; fitting and goodness-of-fit use every integer from 1 to the dataset's full N, including unobserved values.

## Models and fitting

Reuse the likelihood and optimisation routines from the existing local-support family comparison (copied unchanged into distribution_helpers.jl). Fit each dataset and quantity separately by maximum likelihood:

- Bounded exponential: p(k) proportional to exp(-lambda*k), lambda≥0.
- Bounded power law: p(k) proportional to k^(-alpha), alpha≥0.
- Exponential-cutoff power law: p(k) proportional to k^(-alpha)*exp(-lambda*k), both nonnegative.
- Discrete lognormal kernel: p(k) proportional to exp(-alpha*log(k)-beta*log(k)^2), beta≥0. For beta>0 this corresponds to density-kernel parameters sigma=sqrt(1/(2beta)), mu=(1-alpha)/(2beta). This is not a rounded continuous lognormal. beta=0 is a boundary power-law kernel and should not be described as a finite fitted lognormal.
- Singleton plus geometric: p(1)=s and p(k|k≥2) proportional to r^(k-2), 0≤r≤1.

Every distribution is normalised over 1,...,N, not truncated at the observed maximum. Nominal parameter counts are one for exponential/power law and two for the others. Boundary estimates are flagged. AIC and small-sample AICc are reported; rankings are within a dataset and quantity only. No tail-only selection or extra models.

Goodness-of-fit: observed discrete KS distance; simulate 499 independent samples of the original size from each fit, refit each replicate and recompute KS. p=(1+number of simulated KS≥observed KS)/(500). Report p Monte Carlo standard error and optimisation residuals. These are independent-observation working diagnostics: species/pairs share sites and partners, so the tests do not establish independence-robust ecological significance or a universal family. Small SP samples can have low discriminatory power. AIC winners are not automatically adequate fits.

## Outputs

1. 01_CO_fitted_ccdf.png: co-occurrence data and fitted curves, 2×5.
2. 02_SP_fitted_ccdf.png: species occupancy data and fitted curves, 2×5.
3. 03_CO_AIC_goodness_of_fit.png: delta-AIC and refitted bootstrap p for every family/dataset.
4. 04_SP_AIC_goodness_of_fit.png: same for species occupancy.
5. 05_K_CO_SP_empirical_ccdf.png: three empirical curves together per dataset, 2×5.

fit_summary.csv includes absolute AIC/AICc, log likelihood, delta values, KS, bootstrap p, parameter estimates, sample sizes, boundaries and convergence diagnostics. best_models.csv and FINDINGS.md summarise rankings and adequacy. empirical_distributions.csv, fitted_distributions.csv and observations.csv contain the numerical data. Input/source/code hashes and relationship metadata retain provenance.

Confirmed occupancy versions are checked and every C is reconstructed from species occurrences. Existing data-quality flags, incomplete detection, interaction-derived presences, site-count/effort differences and related sampling campaigns remain limitations. See FINDINGS.md. No inference of unseen species or pairs is made.
