# Site accumulation: regional diversity and local support

Separate analysis comparing empirical site accumulation with existing opportunity-constrained Common and persistent-Pair allocations. Previous outputs and parameters are unchanged.

## Run

From RALSA:

```text
julia --startup-file=no --project=. scripts/site_accumulation/run_analysis.jl
```

Alternatively use run.cmd in this folder. Optional arguments specify the repository root, output directory, and `plot` to reuse saved numerical results. Uses the existing RALSA Julia environment and conditional-support helper. Default outputs: `outputs/statistical description of local support/opportunity_theory/site_accumulation`. Re-running replaces this analysis's own outputs only.

## Design

For each of ten datasets, reconstruct the interaction-by-site matrix from validated records and the compatible opportunity sites from confirmed occurrences. Verify row sums K/C and input hashes against the previous opportunity analysis. Retain N, known regional links L, full-region total support M and each link's allowed sites.

Use unchanged entire-dataset-holdout Common and Pair parameters. Their recorded training sets exclude the target. For each of 300 replicates:

1. Draw one random ordering of all sampled sites, shared by empirical and both model trajectories.
2. For each model, draw a full positive support vector conditional on exact M using the existing dynamic program and independent backward sampling. Pair weights are the original equal mixture of separately zero-truncated binomials at logistic(a-h), logistic(a+h). Common has uniform probability on allowed positive-row binary allocations once M is fixed.
3. Given each pair's sampled K, choose that many distinct sites uniformly from its allowed opportunities. This is the original pair/common model's site placement conditional on row support. No shared-site model or extra spatial autocorrelation is introduced.
4. Add sites in the shared order and calculate every integer prefix n=1,...,N: accumulated link richness L_n, accumulated support M_n, mean positive support M_n/L_n, Gini, singleton fraction and support share held by the highest-supported 10% of links.

Empirical matrix stays fixed; model allocations vary between replicates. At every prefix check the mean-support identity. Full prefixes retain every regional link and exactly M. Intermediate M_n and site richness are not constrained to equal empirical values. Their expectations agree (nM/N), but their realised distributions can differ.

Positive-support metrics exclude undiscovered links. An empty prefix has L_n=M_n=0 and undefined mean support/concentration, and is counted explicitly. Gini is the finite empirical population formula. Top-10% share uses fractional rank weighting. Means of nonlinear statistics average replicate statistics, never substitute the statistic of the averaged distribution or the ratio of averages.

## Exact expectations and verification

For a link supported at K of N sites, the probability of being discovered in n uniformly sampled sites is 1-choose(N-n,K)/choose(N,K), with the fraction zero if N-n<K. Stable product calculations reuse the existing random-removal helper. Empirical K or exact conditional model support marginals give analytical expected L_n. Expected M_n=nM/N. The richness figure uses these exact means; replicate summaries supply variability bands. Model full-region support frequencies must reproduce the original fixed-M outputs. Input and endpoint checks are recorded.

Site accumulation is the complementary view of random removal. In particular, expected richness under random order is determined entirely by the full K distribution, so its reproduction is not independent validation of a new mechanism. Joint fluctuations and nonlinear concentration summaries can also reflect placement among sites. This is random accumulation of sampled communities, not chronological ecological assembly, contiguous geographic expansion or a direct beta-diversity decomposition.

## Outputs

- 01_regional_link_accumulation.png: exact expected fraction of regional links discovered and replicate variability.
- 02_mean_support_accumulation.png: mean positive support along accumulation.
- 03_concentration_accumulation.png: support Gini along accumulation.
- 04_diversity_and_concentration.png: concentration versus expected regional link richness; small/medium/large points mark nearest 25%/50%/100% site counts. This is a trajectory of averages, not a fitted relationship between independent samples.
- 05_singleton_accumulation.png: singleton prevalence along accumulation.
- accumulation_draws.csv: every replicate and integer prefix; this is the largest output.
- accumulation_summary.csv: means, central 95% intervals, Monte Carlo SEs and empty-prefix counts.
- exact_expectations.csv: analytical richness and total-support expectations, including n=0.
- paired_differences.csv: empirical minus model differences using shared site orders; intervals include model allocation and site-order variation.
- stage_summary.csv: results near 25%, 50% and 100%, with actual integer counts.
- FINDINGS.md, checks.csv, dataset_inputs.csv and provenance CSVs.

Empirical bands reflect random-order variation only. Model bands reflect allocation plus random-order variation. Neither is a confidence interval or parameter-uncertainty estimate. A different band width is not by itself evidence of model failure. Corresponding prefixes are paired within dataset and replicate; adjacent steps are dependent. Each dataset is reported separately; none is excluded.

## Retained limitations

Known identities and occurrence records are supplied, observed M is imposed, and no amount prediction is attempted. Confirmed occurrence restorations are retained. Uncertain excluded restorations remain Gottin PP's 474 missing-identity records, four Garraf PP2 zero columns and one Montseny zero row. Much occupancy derives from interactions. Prior fictitious Gottin HP V1 removal and unresolved one-link differences in Garraf HP/Olot remain flagged. N does not represent undocumented sites. Shared Gottin sites and related provider collections mean datasets are not fully independent. Sampling effort, observation, spatial coverage and temporal differences remain. No new ecological mechanism is inferred from a better curve fit alone.
