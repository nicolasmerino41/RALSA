# Singleton and bounded-geometric support fits

Run from the RALSA root:

```
julia --startup-file=no --project=. scripts/two_parameter_support_fit.jl
```

This separate analysis writes only to `outputs/statistical description of local support/two_parameter_fit`. It uses Julia, CSV, DataFrames and CairoMakie with the existing project dependencies. No co-occurrence predictors, other distribution families or removal-gap predictions are included.

## Data and versions

The input is the descriptive analysis's `link_support.csv` (all 4,304 observed regional links, K>=1) and `site_universe.csv`. K is the number of distinct supporting sites for a consumer–resource pair. N retains the descriptive analysis's site universe, the union of site IDs in the validated interaction, occupancy and co-occurrence exports. Entirely unrecorded sampled sites cannot be inferred from these files.

The script checks descriptive input hashes against current validated exports and reconstructs every K directly from interaction triples, checking all values and site counts. `input_versions.csv` copies the descriptive source fingerprints; `descriptive_input_versions.csv` fingerprints the actual descriptive inputs. The frozen source exports were imported on 2026-09-29.

Existing flags remain: the fictitious Gottin HP V1 consumer was removed; one-link discrepancies with earlier degree tables for Garraf HP and Olot are unresolved. Confirmed host-occurrence restorations have been tested in the separate zero-presence analysis but have not replaced primary occupancy exports. These affect co-occurrence-based quantities, not K when interaction records and the site universe are fixed. Neither the restored nor original inputs establish complete independent presence sampling for both guilds.

## Models and estimation

Two-parameter model:

P(K=1)=s; P(K=k)=(1-s) r^(k-2) / sum(r^j,j=0,...,N-2), for k=2,...,N.

One-parameter comparison:

P(K=k)=r^(k-1) / sum(r^j,j=0,...,N-1), for k=1,...,N.

Both constrain 0<=r<=1. r=0 puts all geometric-component mass at its lowest support; r=1 gives a uniform component. Larger r means slower decay. The baseline's reported s is its implied singleton probability, not a second fitted parameter. The two-parameter family includes the one-parameter family.

Maximum likelihood gives s=observed singleton fraction. Interior r solves the finite geometric component's expected offset equal to its observed average offset; the script uses monotone bisection. If the observed component mean is at least the uniform-component mean, the constrained optimum is r=1. A zero offset gives r=0. No unconstrained r>1 solution is used.

N=1 has no tail; N=2 has only one possible non-singleton support, so tail r is unidentifiable. With no non-singletons, r and their implied conditional mean are reported NaN (undefined), and two-parameter AIC is suppressed when N>2. Fewer than five non-singletons are explicitly flagged as weakly informative. Identifiable boundary estimates are retained and flagged; nominal AIC at boundaries merits caution. Actual datasets have sufficient non-singletons; see fit_status and counts for each result.

Log likelihood sums the probabilities of individual observed link supports, without a histogram multinomial coefficient. AIC=2p-2 log likelihood, with p=1 or 2 for the regular cases. delta_AIC is relative to the better of these two models within a dataset. Exact recovery of singleton fraction and mean by an interior two-parameter fit is a property of the fit, not evidence of distributional adequacy. The one-parameter interior MLE also matches the observed mean.

## Adequacy diagnostics

The script measures maximum absolute CDF discrepancy over all K, plus the same discrepancy conditional on K>=2. Each model/dataset receives 999 parametric bootstrap replicates, seed 20261012; each replicate samples the same number of links and refits the parameters. The Monte Carlo p-value is (1 + number of simulated discrepancies >= observed)/(B+1). The conditional-tail check is undefined if fewer than five non-singletons are available; its valid replicate count is saved. A minimum reported p of 0.001 is a simulation-resolution limit, not an exact probability.

These checks assess a working independent-link distribution model. Shared species and sites may violate independence; p-values are descriptive diagnostics, not calibrated evidence for ecological mechanisms. They are unadjusted across datasets/models. Non-rejection does not establish a true family, and large datasets can detect small discrepancies. Read effect sizes and the repeated-support tail alongside AIC and bootstrap results.

## Outputs

- `fit_summary.csv`: parameters, conditional non-singleton mean, log likelihood, AIC, boundary/status flags, singleton fraction, mean and p90 (observed/fitted), fraction strictly above the observed p90, CDF discrepancies and bootstrap diagnostics.
- Both observed and fitted p90 use the discrete inverse-CDF definition: the smallest K with cumulative probability >=0.9. This differs from interpolated descriptive quantiles, deliberately using the same definition for both fitted and observed distributions.
- `fitted_distributions.csv`: observed and fitted frequencies and CCDFs for every K=1,...,N, including empty bins.
- `01_frequency_fits.png`: exact empirical frequencies as dataset-coloured points/stems; two-parameter fit in the same colour and one-parameter fit dashed grey.
- `02_ccdf_fits.png`: the same comparisons on a logarithmic probability axis to expose repeated-support tails. Zero probabilities are omitted, not replaced by pseudocounts. Displayed x ranges end up to two sites beyond observed maximum support, bounded by N; fitted probabilities are computed over the full 1,...,N range. Plotting bounds do not truncate the models. Lines connect discrete fitted values.

Figures retain the descriptive dataset colours and contain no explanatory footers. `FINDINGS.md` interprets adequacy and remaining discrepancies. No power-law claim follows from straightness on any log view. The singleton component is a statistical distinction, not evidence for a separate ecological process.
