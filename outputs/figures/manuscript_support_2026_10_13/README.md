# Separate manuscript figure set

This set uses existing confirmed-restoration, exact-envelope, exact-allocation and predictive-validation outputs. It does not refit models, simulate removal, or replace current manuscript figures.

Run from the repository root:

```powershell
julia --startup-file=no --project=. figures/manuscript_support_2026_10_13/build_figures.jl
```

Optional arguments are the repository root and output directory. Default outputs go to `outputs/figures/manuscript_support_2026_10_13`. The script requires the existing project environment and source outputs. It verifies saved input fingerprints, support and opportunity records, exact empirical retention, arithmetic averaging, saved allocation summaries, and predictive error summaries before rendering. CSV files alongside the PNGs expose display transformations and validation results. Captions are in `CAPTIONS.md`; panel sources are listed in `PROVENANCE.csv`. The fixed dataset order follows the current project configuration.

## Input compatibility and limitations

- K and the site universe come from the descriptive local-support analysis; their hashes are recorded by the later opportunity analysis. Confirmed restored occurrences are used to check C against the opportunity analysis. Thus panels combine the same interaction supports, site sets and confirmed opportunity version.
- The confirmed version restores 629 presence records: Quercus 156, Gottin HP 308, Garraf HP 115 and Olot 50. Uncertain restorations remain excluded: 474 Gottin PP records with missing plant identities, four Garraf PP2 zero columns and one Montseny zero row. Much occupancy is interaction-derived rather than independently surveyed. These limitations particularly affect opportunity-constrained and predictive results.
- The fictitious Gottin HP V1 record was removed previously. Outstanding one-link differences in Garraf HP and Olot remain as documented in the source analyses. Occurrence-only corrections do not change K unless interaction records or the site universe also change. N describes recorded sampled sites; it does not recover wholly undocumented sites.
- Gottin HP and PP share 32 sites. Garraf HP, PP, PP2, Olot and Montseny share provider/collection context; exact overlap remains unresolved. Neither the figure nor its caption claims independent landscape validation.
- Mean gap means an arithmetic average across integer removal levels up to 80%, including zero, not a continuous area. Only Figure 2B uses interpolated curves, solely for display.
- The allocation intervals are conditional allocation variability, not statistical confidence intervals. Prediction models and saved draws are used unchanged. The within-system task completes a partly observed full landscape and must not be interpreted as the same information setting as entire-dataset holdout.
- No required panel was unavailable. Scientific outcomes and model choices are unchanged.
