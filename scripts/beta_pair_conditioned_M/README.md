# Continuous beta pair heterogeneity, with total support fixed

This separate Julia/CairoMakie test replaces the two equally likely pair probabilities with a continuous beta distribution. It retains known regional link identities, positive support for every link, the same confirmed opportunities, observed total support M, and the same entire-dataset holdouts. The saved two-class results are reused unchanged. No Site model is added.

Run from RALSA:

```powershell
julia --startup-file=no --project=. scripts/beta_pair_conditioned_M/beta_pair.jl
```

Optional arguments specify repository and output directories; a third argument `plot` regenerates figures and the findings note from saved results without refitting. Default output is `outputs/statistical description of local support/opportunity_theory/beta_pair_conditioned_M`.

## Precisely which beta model?

For a known regional pair, p has distribution Beta(mu*phi, (1-mu)*phi). Mu controls the mean latent probability and phi its concentration; lower phi means greater heterogeneity. Given p and the pair's C opportunity sites, K follows Binomial(C,p) conditional on K≥1. Thus

```
w_C(k) = integral Beta(p; mu*phi, (1-mu)*phi)
                  * Binomial(k; C,p) / [1-(1-p)^C] dp,   k=1,...,C.
```

This is the continuous counterpart of the existing selected-pair two-class model: the beta distribution is specified **among known pairs**, and row positivity follows the probability draw. Strictly, it is a beta mixture of zero-truncated binomials, rather than the conventional beta-binomial truncated after mixing. Those operations differ: ordinary post-mixture truncation reweights the beta distribution toward higher p, depending on C. Here the same conditioning order is retained so the comparison changes heterogeneity shape without also changing the selection convention. It should be labelled **Beta pair** in figures, not presented as a standard zero-truncated beta-binomial.

The joint support distribution is then conditioned on the target's exact M:

```
P(K_1,...,K_L | M,O) is proportional to product_i w_Ci(K_i),
subject to 1≤K_i≤C_i and sum_i K_i=M.
```

M is additional target information, but individual K values are not predictors. As before, gap expectations depend on row supports; conditional on K, supporting sites are uniformly chosen among that pair's opportunities.

## Training and evaluation

- Hold out one entire dataset. Fit using the other nine with equal dataset weights. No target K, M, singleton fraction or gap is used to select parameters. Target M enters only the conditional allocation step after training.
- Preserve the earlier four summaries, regularised Gaussian synthetic-likelihood score and 64 simulation replicates per dataset/grid point: log(mean K+0.1), singleton fraction, log(mean K²+0.1), and log(1+site residual clustering). Site residuals use the same opportunity-conditional expectations and variances as the original script.
- Prespecified grid: logit(mu)=−4,−3.5,...,4; phi=0.05,0.1,0.25,0.5,1,2,4,8,16,64,256. Boundary selections are reported, not silently expanded using held-out performance. Both pair families have two parameters, but their grids and numerical simulation libraries differ; close performance rankings should not be overinterpreted.
- Integrate the beta density with 128-node Gauss–Jacobi quadrature. Check each selected fit with 256 nodes. Quadrature moments, an analytic uniform-beta example and probability normalisation are verified. This controls numerical integration error, not model adequacy.
- Reuse the exact-total dynamic-programming method from the fixed-M analysis, copied as `conditional_support.jl` without executing or changing that analysis. Compute conditional marginal support frequencies and mean curves deterministically. Generate 2,500 independent conditional draws per dataset for allocation intervals.
- Use the exact random fixed-count-removal survival law. Compute all integer levels through complete removal, with primary means and errors over m/N≤0.8 including zero. These use the same averaging as the saved two-class results.
- Additional unconditioned beta curves are exported as a diagnostic only. The requested comparisons and all primary figures fix M for both pair models and Common.

The empirical opportunity records are shared across comparisons. The original beta-binomial pilot included all co-occurring pairs (including K=0) and fitted likelihoods of counts. Its parameters are not reused: that would compare a different pair population and fitting objective.

## Files

- `01_allocations.png`: empirical mean gap and singleton percentage versus Common, two-class Pair and Beta Pair means and central 95% conditional allocation intervals.
- `02_two_class_vs_beta.png`: paired errors for the two heterogeneity descriptions, including equal-dataset means. This is the direct comparison requested.
- `03_gap_curves.png`: empirical and predicted gap curves with the same axes.
- `summary.csv`, `across_dataset_summary.csv`, `curves.csv`, `support_frequencies.csv`, `draws.csv`: numerical results.
- `fitted_parameters.csv`, `training_grid_scores.csv`, `simulation_library.csv`: training audit and reusable simulation reference.
- `algorithm_checks.csv`, `quadrature_checks.csv`, `result_checks.csv`: numerical checks.
- `input_versions.csv`, `code_versions.csv`, copied relationship metadata: provenance.
- `FINDINGS.md`: results and interpretation.

Intervals describe possible allocations under fixed constraints, not confidence intervals or removal-sequence variation. Predictive performance can establish whether these summaries are sensitive to the choice of probability distribution; it cannot establish true ecological classes or identify a causal mechanism.

## Retained limitations

The same confirmed restored occurrences, interaction supports and site universe are used. Uncertain restorations remain excluded: Gottin PP 474 missing-identity cases, Garraf PP2 four zero columns and Montseny one zero row. Occupancy is often interaction-derived; the prior Gottin HP V1 removal and unresolved one-link differences in Garraf HP/Olot remain documented. Occurrence-only changes affect opportunities, not K unless interactions or the site universe also change. Wholly undocumented sites and interactions are not estimated.

Gottin HP and PP share sites; provider-collection overlap among Garraf/Olot/Montseny datasets is unresolved. Therefore these holdouts are exploratory rather than independent landscape validation. The Montseny Site sampler flag in the previous analysis is irrelevant to this comparison: only its exact Common and Pair results are reused.
