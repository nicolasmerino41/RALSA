# Support prediction: technical pilot

Run from RALSA: `julia --startup-file=no --project=. scripts/support_prediction.jl`.
Use `--plot-only` to regenerate PNG and vector PDF figures from the saved tables.
No new Julia dependencies are required. Existing manuscript outputs and inputs are untouched.

## Input status

These are the uncorrected 2026-09-29 exports. Confirmed independent host-presence records have not yet been restored for Quercus, Gottin HP, Garraf HP and Olot. Salix–Galpar includes some independent host records. The current consumer presences are interaction-derived; other unresolved raw-record cases require review before corrections. Consequently these results are a technical pilot, not a final ecological test.

`input_versions.csv` records SHA-256 hashes of raw and processed files. All processed records were checked against the original raw exports and matched. `occurrence_audit.csv` records the present input status and occupancy beyond linked species-site records. This records the exports' limitations; it does not establish absence of independent observations in the original studies. After correcting inputs, review this audit and the raw-checksum safeguard in the script before rerunning.

## Models and validation

The observational unit is a consumer–resource pair that co-occurs at least once. n is its number of co-occurring sites; K is its number of interaction sites, including K=0. N is the number of sites in the union of occupancy, co-occurrence and interaction tables. Species identities and guild roles remain separate in `pairs.csv`.

M0 assumes K|n ~ Binomial(n,p). M1 assumes K|n ~ Beta-binomial(n,mu*phi,(1-mu)*phi). Its underlying pair propensities share a beta distribution independent of n. Smaller phi allows more heterogeneity; large phi approaches M0. Neither model estimates a pair-specific parameter from held-out interactions.

The objective averages pair log probabilities within each training dataset, then averages those dataset values equally. M0 has an analytic optimum. M1 uses transformed mu and phi and bounded Nelder–Mead optimisation from eight starting points. Saved convergence and boundary diagnostics accompany the parameters. Logit(mu) bounds are [-12,12]; log(phi) bounds are [-10,log(1e8)].

LODO is the primary leave-one-dataset-out test. Predictions use only the target's n and N; target K is used for scoring and observed curves. Within-dataset fits are descriptive, not predictive validation. Grouped holdouts additionally exclude Gottin HP and PP together, or Garraf HP, Garraf PP, Garraf PP2, Olot and Montseny together. The latter is a conservative shared-provider collection block, not proof of one sampling programme or exact site overlap. The original processing scripts identify Ingo Grass for Gottin and Jordi Bosch for the Garraf–Montseny–Olot collection. Quercus, Nahuel and Salix–Galpar each form a separate block, so their grouped results equal LODO. Summaries weight datasets equally, including in the grouped analysis.

## Predictions and metrics

Expected support counts are H[k] = sum over pairs of P(K=k|n), including zero. Predicted initial link richness is L=sum(H[k], k>=1); total support is T=sum(k*H[k]). Fractions among realised links divide positive H by L.

For m uniformly removed sites, a link with support k survives with probability 1-C(m,k)/C(N,k), with extinction probability zero when m<k. Expected regional retention is sum(H[k]*survival)/L: a **ratio of expected counts**, not the expectation of a random ratio. Observed curves use the empirical support histogram in the same expression. Local retention is 1-h. At h*N between integers, interpolate the floor/ceiling removal-count expectations. Evaluate h=0,0.02,...,0.8; the mean gap is the trapezoidal area between regional and local retention divided by 0.8. Cross-link spatial overlap is unnecessary for these expectations, although it matters for variability.

`metrics.csv` contains per-dataset NLPD (lower is better), total support and link-richness relative errors (prediction/observation - 1), singleton fractions, mean support per realised link, retention-curve mean absolute error and signed mean-gap error. `equal_dataset_summary.csv` averages absolute errors across datasets. Proportional errors and gap values are stored as fractions, not percentages. Within-dataset M0 recovers T by construction; this is not predictive evidence. M1's likelihood fit does not enforce that moment.

The conditional diagnostic compares mean K/n among K>0 pairs with mu/P(K>0|n). Even M0 gives a declining curve because restricting to realised links selects successes, particularly at small n. This diagnostic alone cannot establish a specialist strategy; n is not regional degree.

## Figures

All figures have matching PNG/PDF versions. Orange denotes M0 and blue M1. Black denotes observed values. Figures 01–03 and 06–07 use LODO predictions.

1. **01_transfer_accuracy**: observed versus predicted T, L, mean gap and singleton fraction. Each point is one dataset; dashed lines show equality. Dataset-specific errors are identified in figure 05.
2. **02_support_distributions**: complementary cumulative distribution P(K>=k|K>0), showing the prevalence of repeatedly supported links. Histograms including zeros are in `support_histograms.csv`.
3. **03_retention_curves**: regional retention compared with the grey local-retention baseline. The area between each regional curve and that baseline is the gap.
4. **04_transfer_vs_local_fit**: NLPD under LODO, grouped holdout and within-dataset fitting, on a shared y scale. Within-dataset points use target interactions.
5. **05_error_diagnostics**: signed errors for every dataset and validation scheme. Negative values are underestimates. Count errors use percent; gap errors use percentage points. Colour ranges differ by metric.
6. **06_conditional_realisation**: observed mean K/n among realised pairs at each n and model conditional expectations. Point diameter increases logarithmically with the number of realised pairs contributing to that mean.
7. **07_regional_realisation**: fraction of all co-occurring pairs with K>0 at each n. Point diameter increases logarithmically with the number of co-occurring pairs. Read with 06: matching the selected positive pairs need not reproduce how many pairs realise any interaction.

These lines are analytic expectations, not confidence bands. This pilot estimates neither parameter uncertainty nor sampling/detection uncertainty. Interactions sharing species or sites need not be independent; the likelihood is a working marginal model, not a test of ecological independence.

## Numerical checks

All 50 fits converged; none reached parameter boundaries. The maximum spread in M1 objective values across eight starting points was 2.45e-14. Independent Python log-gamma likelihood calculations agreed within 1.56e-15, and local parameter perturbations did not improve the fitted values. Julia assertions verify known binomial/beta-binomial probabilities, normalisation, the binomial limit, expected total support, retention monotonicity, the local-retention bound and within-dataset M0's exact total. This checks numerical implementation, not model adequacy.
