# Persistence of interaction propensities

Run from RALSA:

```text
julia --startup-file=no --project=. scripts/propensity_persistence/run_analysis.jl
```

Optional arguments: repository root, output directory, and `plot` to redraw saved results. Default outputs are `outputs/statistical description of local support/opportunity_theory/propensity_persistence`. Only this analysis's outputs are replaced on rerun.

## Controlled experiment

Use confirmed opportunities C for every known regional link, observed total support M, site count N, and positive support for all L known links. No parameters are fitted. For support statistics, opportunity counts suffice: conditional on a row's K, its supporting sites would be selected uniformly among its actual opportunity sites. No additional geographic process is simulated.

Hold the latent mean probability at 0.5 and vary rho over 0, 0.1, 0.3, 0.6, 0.9. This symmetric midpoint is a theoretical control, not an empirical realisation estimate. Conditioning on M does not generally eliminate dependence on the latent mean for heterogeneous models; this exercise does not claim robustness to all choices of that mean.

- Independent local variation: draw a fresh probability for each opportunity. Independent Bernoulli observations marginalise to common probability 0.5 regardless of its distribution. Conditional on positive links and M, this is exactly the existing Common allocation null. There is therefore one baseline, not redundant simulations for each heterogeneity level.
- Persistent two-state: each pair has probability (1−sqrt(rho))/2 or (1+sqrt(rho))/2, equally likely, shared across its opportunity sites.
- Persistent beta: each pair has probability Beta(a,a), a=(1/rho−1)/2, shared across its opportunity sites. At rho=0 use constant probability 0.5.

Both persistent families have mean 0.5 and variance rho/4 before conditioning. For C opportunities, unconditioned variance is C/4 + C(C−1)rho/4; independent local variation has variance C/4. Rho is the pre-conditioning correlation between binary observations for a pair. Neither this variance formula nor the interpretation of rho as a correlation is asserted after conditioning.

First integrate over latent probabilities to obtain binomial-mixture or beta-binomial row weights. Then condition jointly on all rows positive and sum(K)=M. This differs deliberately from the historical fitted Pair/Beta Pair models, which mixed separately zero-truncated binomials. These scenarios test a generative principle rather than rerun or rank the historical fitted models.

Exact dynamic programming produces marginal frequencies and independent backward draws, with no MCMC or rejection approximation. Use 300 draws per scenario for nonlinear statistics and central 95% allocation intervals. Seed is fixed. The baseline must reproduce saved Common frequencies; each draw must obey C and M. A small enumerated example checks the dynamic program against exhaustive probabilities for every family and heterogeneity level. Source hashes must match the authoritative input-version table.

## Outputs

- 01_singletons.png: singleton proportion versus heterogeneity.
- 02_concentration.png: support Gini versus heterogeneity.
- 03_conditional_variation.png: link-weighted within-C population variance of K/C, restricted to C classes with at least three links.
- 04_support_ccdf.png: exact expected CCDFs, with rho=0.6 prespecified as an illustration; no best-fit selection.
- summary.csv: means, empirical statistics, central 95% allocation intervals and Monte Carlo standard errors.
- frequencies.csv: exact expected support frequencies and empirical frequencies at every K.
- draws.csv: all replicate singleton, Gini, top-10% support share and within-C variance statistics.
- checks.csv, source_versions.csv, input_versions.csv and relationship metadata: numerical verification and provenance.
- FINDINGS.md: results and retained data limitations.

The highest 10% support share uses fractional ranking; Gini uses the population formula. Bands are allocation intervals, not confidence intervals, parameter uncertainty or removal-order variation. No gap predictions or new distribution fitting are added.

Persistent differences may arise from abundance, compatibility or environmental context; the experiment cannot distinguish those biological explanations. Fixed known identities and M make this a conditional organisation test. Sampling and opportunity-record limitations remain; see FINDINGS.md.
