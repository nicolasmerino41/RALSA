# Opportunity and support organisation

A separate Julia/CairoMakie analysis of whether uneven support is inherited from opportunities or reshaped by interaction realisation. Existing analyses and model parameters are preserved. No additional distribution families or transferable prediction rules are fitted.

## Run

From the RALSA root with its existing Julia environment:

```text
julia --startup-file=no --project=. scripts/opportunity_support_organisation/run_analysis.jl
```

Output folder: `outputs/statistical description of local support/opportunity_theory/opportunity_support_organisation`.

Optional arguments are repository root, output directory and `plot` (reuse saved numerical outputs). Both Julia files must remain together. The existing conditional-support helper at scripts/beta_pair_conditioned_M/conditional_support.jl is reused. The code requires no new dependencies beyond the RALSA project. Running again replaces only this analysis's own outputs.

## Inputs and questions

Use the confirmed-restoration opportunity-theory quality_control.csv, fitted_parameters.csv, site universe, restored occupancy and validated interaction records. K is the number of distinct supporting sites for an observed regional link; C is its number of co-occurring sites. All full-data links satisfy 1<=K<=C<=N. Interactions and opportunities are rebuilt from site records and checked against saved counts. Original input hashes must match; source_versions.csv records every additional input. Known interaction identities are assumed; wholly unobserved interactions are outside scope.

1. Compare C and K frequency distributions, CCDFs and concentration among the same known regional links.
2. Hold N, known links, each C and observed M fixed. Compare the existing Common and persistent-Pair allocation models using exact expected support frequencies and 499 independent exact conditional draws per model/dataset. Means and central 95% allocation intervals describe singleton fraction, support Gini, top-10% support share and within-C realisation spread. No model is refitted.
3. Group links by exact C. Compare empirical mean and population SD of K to allocation distributions calculated for the identical groups. Display only classes with at least three links; retain all classes and group sizes in CSV. The aggregate spread is the square root of link-weighted within-class population variance of K/C, calculated over eligible classes, displayed in percentage points. Model points transform mean variance; they are not the mean of simulated SDs.
4. Draw 300 uniform site permutations per dataset and recompute K/C at 4, 8, 14, 20, 32, 64, 128, 256 sites where available, plus the full N. The primary common comparison uses eight sites, the minimum N across datasets. Summaries condition on K>0 in the subsample; matched C distributions use those same links. Separately retain fractions of original links with K=0 and with K=0,C>0. Full-site samples reproduce inputs. Nahuel has no variation at eight sites because all its sites are retained.

## Allocation models

The Common model has one probability across eligible cells. Conditional on positive support for every known link and exact M, this probability cancels, giving the uniform eligible positive-row binary allocation null. Unequal C alone can create unequal K under this model.

The Pair model uses unchanged entire-dataset-holdout parameters, learned on the other nine datasets. Each known pair initially has equal low/high latent probability logistic(a-h), logistic(a+h). The two separately zero-truncated binomial components form the existing mixture. Conditional allocation weights are products of those mixtures restricted to sum(K)=M. Exact dynamic programming provides expected support frequencies and independent backward draws. A numerical tilt cancels on conditioning and is not a fitted parameter. Small-system enumeration verifies exact marginal calculations; full distributions reproduce the prior conditional-M results.

Intervals describe allocation variation, not confidence intervals, parameter uncertainty or site-removal variation. Empirical positions and interval exceedances are descriptive and not multiplicity-adjusted hypothesis tests. Shared records and sampling affect empirical summaries. Pair improvement supports a statistical description of persistent heterogeneity; it does not prove two biological classes or distinguish ecological heterogeneity from heterogeneous detection.

Concentration uses sorted ascending link support (or C); the top 10% share uses fractional weighting at rank boundaries. Model Gini, top-share and Lorenz summaries are calculated for each allocation before averaging. Opportunity and support concentration curves rank each quantity separately, so they are not a statement about the identities of the highest-ranked links.

## Outputs

| Figure | Content |
|---|---|
| 01_opportunities_and_support.png | Full-data C/K CCDFs and conditional model expectations |
| 02_allocation_features.png | Four empirical features versus allocation means and 95% intervals |
| 03_support_at_equal_opportunity.png | Mean K among links sharing C |
| 04_variation_at_equal_opportunity.png | SD of K among links sharing C |
| 05_concentration.png | Opportunity/support Lorenz curves and allocation intervals |
| 06_sampling_sensitivity.png | Full-data versus eight-site summaries and link disappearance |
| 07_eight_site_comparison.png | Mean C/K CCDFs on the matched eight-site samples |
| 08_support_frequencies.png | Full-data frequencies and conditional expectations |

FINDINGS.md reports the results and limits. empirical_summary.csv, distributions.csv, concentration_curves.csv, allocation_summary.csv, conditional_C.csv, site_subsampling_summary.csv and eight_site_distributions.csv underlie the figures. allocation_draws.csv and site_subsampling_draws.csv retain all simulated summaries. input_checks.csv and algorithm_checks.csv record validation. Source hashes, code hashes and parameters_used.csv document provenance.

## Interpretation limits and retained data flags

This tests recurrence and candidate explanations for distributional features, not universality or a particular tail family. Matching site number does not match within-site effort, spatial coverage, detection, temporal coverage or spatial clustering. The eight-site support range is short. Subsampling selects surviving links and changes the cohort; lost links are reported rather than silently forgotten. Models are not rerun on subsamples, so the sampling check addresses empirical descriptive robustness, not sampling robustness of model preference.

Confirmed occurrence restorations are retained. Uncertain restorations remain excluded: Gottin PP 474 missing-identity records, four Garraf PP2 zero columns and one Montseny zero row. Much occupancy is interaction-derived. Prior fictitious Gottin HP V1 removal and unresolved one-link differences for Garraf HP/Olot remain. N cannot recover wholly undocumented sites. Occurrence-only corrections can change C without changing K. Gottin HP/PP share sites; related collections have unresolved overlap. These are not ten independent landscape validations. Metadata relationship files are copied unchanged. No original output is overwritten.
