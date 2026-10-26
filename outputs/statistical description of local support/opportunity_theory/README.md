# Opportunities, support and the retention gap

This separate analysis tests exact constraints and three small generative scenarios. It preserves the existing analyses and does not refit the earlier singleton-plus-geometric distributions. Start with `FINDINGS.md`; mathematical details are in `BOUNDS.md`.

## Run

From RALSA:

```
julia --startup-file=no --project=. scripts/opportunity_theory.jl
julia --startup-file=no --project=. scripts/verify_opportunity_theory.jl
```

The script also accepts repository and output-directory arguments. Outputs default to `outputs/statistical description of local support/opportunity_theory`. Julia packages come from the existing project and manifest. The independent Python verification is optional; it needs NumPy and pandas. Nothing is written into other analyses. Seeds are fixed in the code: site split 731029; library 10000+dataset index; recovery 30000+index; scenarios 40000+index; prediction 50000+100*index+within-indicator; numerical sampler check 917.

The simulation library is cached in this output directory. The same code, inputs, grid, seed and replicate count must be used with this cache. To run a changed specification, use a different output-directory argument. Optional `OPPORTUNITY_STAGE` values `recovery`, `scenarios`, `evaluation`, and `plots` run those stages with the existing library. An ordinary run executes verification, bounds, recovery, scenarios, predictions and figures in that order.

## Inputs and framework

`data/processed/*_interactions.csv`, the descriptive `link_support.csv` and `site_universe.csv`, and the **confirmed restored occupancy** tables in `outputs/support_prediction_zero/confirmed/inputs` define the matrices. Input hashes from the previous co-occurrence analysis are verified before use; additional interaction-file hashes are recorded. Rows are binary distinct pair–site occurrences. Repeated observations are not interaction abundance, strength, rate or a measured functional response.

`METHODS.md`, the support-prediction, co-occurrence and transfer scripts, the manuscript `Articles/Co-occurrences/CO.docx`, and `The story so far.docx` were inspected. No separately named framework file was found in the project. The latter document is an older framework centred on homogeneous-model comparisons. Its assertions are historical context, not validation of this analysis; the current request defines the present tests. Manuscript bundling sections and older modelling outputs are not brought into the present analysis.

For every originally observed regional link, O is the intersection of its partners' recorded site sets. C is its number of opportunities; A is binary recorded interaction presence; K is its row sum; M is the sum of K. All 4,304 links satisfy **A ≤ O and 1 ≤ K ≤ C ≤ N**. No inconsistencies were repaired, links dropped or C replaced by K. Site splits use the saved descriptive site ordering; `site_split.csv` records the indices. The original site table supplies the corresponding names.

The confirmed version restores 629 species-site presences. The 479 uncertain additions remain excluded: 474 Gottin PP records with missing plant identities, four Garraf PP2 plant zero columns and one Montseny plant zero row. Much occupancy was originally interaction-derived, even after restoration. Known regional identities, O and N are supplied for the full target dataset; this is not an independently measured occupancy validation throughout and not prediction of undetected regional links.

Other retained flags: fictitious Gottin HP V1 removed; unresolved one-link discrepancies for Garraf HP and Olot; N describes recorded sites rather than a census of wholly unrecorded sampling locations; detection and sampling effort are not modelled. The input versions and occupancy audit remain linked to the existing provenance tables. Gottin HP and PP share the same 32 sites. Garraf HP/PP/PP2, Olot and Montseny share a provider collection, with exact overlap unresolved. `metadata_relationships.csv` preserves this distinction; no further dependence is assumed. Ordinary dataset holdouts are not independent landscape validation.

## Removal design and integrated quantity

For exactly m sites removed uniformly without replacement, extinction is choose(m,K)/choose(N,K), zero if m<K. Each curve is an **expectation over removal**, with regional retention divided by original L and local occurrence retention divided by original M. Local retention is 1−m/N. The gap is regional minus local retention. Stable products evaluate survival, without removal simulations.

All integer m=0,…,N are exported. Main errors and the integrated quantity use the same discrete measure: **the arithmetic mean over integer levels m/N≤0.8, including zero**. This is an integral with respect to a uniform distribution over these removal counts, not a continuous trapezoidal area. Different N can give different last included fractions. Independent site removal with probability r instead has gap r−mean(r^K); it is used in the mathematical checks, never substituted into the empirical curves.

## Three conditional generative rules

Each rule retains the actual O and the originally observed link identities. For an eligible row with probabilities p_j, sample its independent Bernoulli opportunities **conditional on at least one success**. Its normaliser is Z=1−product(1−p_j). The sequential sampler evaluates this conditioning directly; it never discards a zero row or changes L.

- **Common:** p=logistic(a) for all eligible opportunities.
- **Pair:** draw one sign z_l∈{−1,+1} with equal probability for each known link; p_l=logistic(a+h z_l) persists across its sites.
- **Site:** draw one sign z_j∈{−1,+1} with equal probability for each site; p_j=logistic(a+h z_j) is shared by all eligible links there.

These are finite, two-level random-effect distributions, not fitted singleton components. h is half the difference between the two log odds. The models have one or two parameters. Pair signs can represent persistent compatibility/encounter differences; site signs can represent context or sampling. They identify neither interpretation.

**Ascertainment convention:** latent sign distributions are specified for the selected population of known regional links and their sites; positive-row conditioning follows the sign draw. Thus the full generative law is P(z) times the product of truncated row laws. It is not the law obtained by simulating all possible species pairs and subsequently selecting whichever happened to interact. In particular, equal mixing weights describe this selected population, not a pre-selection community. This explicit conditional family permits a fixed L and avoids claiming anything about absent pairs. Probabilities p are latent Bernoulli propensities before row conditioning, not the final marginal occurrence probability p/Z.

A combined pair-and-site model is not fitted. The tests already reveal limited transfer, grid-boundary estimates and weak-effect identification limits. Adding a combined scenario now would make attribution harder without resolving those issues.

## Controlled scenarios

Common h=0 and Pair/Site h=1,2,3 are evaluated for all systems. For each scenario, a is solved by bisection to match **expected** M to empirical M. The site expectation is computed exactly by summing over the binomial number of high-propensity sites among a row's C opportunities. Activity is not fixed separately in each realisation. There are 96 fresh matrices per scenario; expected M, realised M and variability are exported. Changes after this control refer to equal expected activity, not equal realised activity. This exercise uses target M and is not prediction.

Scenario and predictive normalised positions use bounds calculated with each generated matrix's own M. This ensures a meaningful comparison when generated activity varies. Empirical bounds always use observed M. No site totals are additionally fixed.

## Estimation and recovery

The parameter grid was fixed before examining held-out errors: a=−4,−3.5,…,4; h=0 for Common and h∈{0.5,1,2,3} for Pair/Site. **Estimation uses an approximate Gaussian synthetic likelihood of summaries, not the full matrix likelihood or AIC.** This makes the three candidate rules comparable without introducing many fitted pair or site parameters.

For each grid point and actual opportunity matrix, 64 independent positive matrices supply a reference mean and covariance. Summaries are log(mean K+0.1), singleton fraction, log(mean K²+0.1), and log(1+site residual clustering). For partially observed training matrices, zero-row frequency is also included. Clustering compares each site total with its expectation if each row's observed support were placed uniformly among that row's opportunities. It averages squared residuals divided by marginal null variance plus one. For masked training data, both row support and opportunities in this diagnostic are restricted to training sites. This support-conditioned diagnostic isolates spatial clustering beyond row totals. These K/C values standardise a training-data residual; they are never assigned as predictive interaction probabilities. Covariance is regularised as 0.8V+0.2diag(V)+10⁻⁵I, fixed in advance.

The loss is (Mahalanobis distance squared + log determinant)/number of summaries. Training datasets contribute equal arithmetic weights. Each model chooses its minimum-scoring grid point; no held-out support, M or gap is used for calibration. The grid is finite and some estimates hit its edge; these are explicitly flagged, not claimed to be continuous optimum estimates. Synthetic likelihood normality and the finite simulation library are additional approximations.

Recovery uses independent simulations on Nahuel, Gottin HP and Salix–Galpar opportunity matrices. It evaluates Common and Pair/Site at h=0.5,1,2, with a=−1 and 24 repeats per condition. Model recovery is the lowest synthetic-likelihood score across candidates, without a complexity penalty. This tends to confuse the common rule with weak heterogeneity and is **not a hypothesis test for heterogeneity**. Parameter recovery is also reported within the correct generating family. Good recovery of these deliberately separated two-state scenarios does not establish identification of arbitrary ecological processes.

## Validation and information supplied

**LODO:** fit shared a,h using observations from the other nine datasets. Supply target O, known link identities and N; supply no target interactions, K, M or fitted support parameters. Draw unknown pair/site signs anew. 400 predictive matrices per candidate give expected support, gap and occurrence probabilities. The original arithmetic-mean s,r transfer and own-fit reconstruction are loaded unchanged as external gap comparisons. The latter is explicitly in-sample.

**Within:** one reproducible random approximately 20% site block is hidden per dataset. Fit the same two-parameter family to training summaries, using simulations masked to those training sites. Rows absent from training remain in the selected link population; they are not discarded. Simulations retain the full O and condition on full-dataset positivity, accounting for known regional identities. Target M is not supplied.

Predictions then condition on observed training cells. Common and Pair integrate the two possible pair propensities exactly using the training-cell likelihood. For Site, latent signs are sampled with two Metropolis chains from opposite initial states, using individual, half-site-block and all-site flips; 800 burn-in sweeps and 1,000 retained sweeps, thinning by five, give 400 matrices. Likelihood includes the probability of at least one success in hidden cells when a row has no training success, divided by the full-row positivity normaliser. All hidden site effects are integrated; none is estimated from hidden interactions. Global moves address the initial local-only sampler's diagnosed trapping. Chain summaries are exported. Final local acceptance is about 5% in Gottin PP and higher elsewhere. The two chains give close mean site probabilities; this is a diagnostic, not proof that every latent mode has been explored. A separate four-site test agrees with exhaustive posterior integration to within Monte Carlo error.

Within-system gap curves describe **completion of a landscape whose other approximately 80% of sites are observed**. They are not directly comparable in difficulty with wholly withheld systems. Cell Brier scores and marginal log loss evaluate only hidden eligible cells. This is one spatial split per dataset, not repeated or geographically separated validation.

## Diagnostics and uncertainty

Outputs include the full support frequency/CCDF, singleton frequency, K≥5 tail, 90th percentile, K conditional on C, M, site residual clustering, gap curves and position within envelope B. Brier and marginal log scores evaluate eligible pair–site cells (equal cell weight within a dataset, then equal dataset weight in aggregate); they do not score the full joint matrix law. Monte Carlo probabilities use (successes+0.5)/(draws+1) to avoid numerical zero/one log scores.

Gap MAE, maximum error and bias are in percentage points. Relative error is 100 × MAE / mean empirical gap, using the same removal levels and units. Bands and diagnostic intervals are 2.5–97.5% **generative prediction intervals conditional on fitted parameters**, not confidence intervals, parameter uncertainty, or removal-sequence variability. Every generated matrix has an exact expected removal curve. Batch Monte Carlo standard errors are exported; Z uses a 100-matrix subsample. No formal coverage claim is made from ten dependent systems. Parameter uncertainty and uncertainty in O are not included.

## Files

- `envelope_curves.csv`, `envelope_summary.csv`, `extremal_support_vectors.csv`: both exact envelopes and empirical positions.
- `controlled_scenarios.csv`, `scenario_curves.csv`: activity-matched scenarios and their realised variation.
- `simulation_library.csv`: reproducible synthetic-likelihood reference simulations.
- `recovery.csv`, `recovery_summary.csv`, `parameter_recovery.csv`: model and parameter recovery.
- `fitted_parameters.csv`, `site_split.csv`, `site_sampler_checks.csv`: calibration and conditional sampling audit.
- `prediction_metrics.csv`, `performance_summary.csv`, `baseline_comparison.csv`: held-out performance and prior comparisons.
- `prediction_curves.csv`, `support_predictions.csv`, `conditional_support.csv`, `prediction_diagnostics.csv`: expected curves, support and diagnostic failures.
- `quality_control.csv`, `input_versions.csv`, metadata files: provenance and constraints.
- `algorithm_checks.csv`, `conditional_sampler_verification.csv`, `verification.json`: exact and independent checks.
- Figures 01–03: envelopes, positions and controlled scenarios; 04–05: transfer and within-system gaps; 06–07: diagnostic and support predictions. Colours identify datasets; line styles identify scenarios.
