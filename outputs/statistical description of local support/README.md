# Descriptive local-support analysis

Run from RALSA:

```
julia --startup-file=no --project=. scripts/statistical_description_local_support.jl
```

Outputs are confined to `outputs/statistical description of local support`. Existing analyses are preserved. This stage fits no models and uses no inferred or simulated interactions.

## Definitions and provenance

An interaction is an ordered (consumer, resource) pair, preserving guild identity. K counts distinct supporting sites per observed regional interaction; only K>0 interactions enter the distributions. Total support is the sum of K, not interaction abundance. N is the project's current sampled-site universe: the union of site identifiers in the validated interaction, occupancy and co-occurrence exports. `site_universe.csv` enumerates it. Any sampled site absent from all three exports cannot be recovered here; N should not be interpreted as an independent sampling-effort audit.

Inputs are `data/processed`, checked against the frozen normalised `data/raw` exports imported on 2026-09-29. The script verifies raw SHA-256 checksums, processed/raw row equality, unique records, nonempty identifiers and 1<=K<=N. `input_versions.csv` records every path, hash and row count. These are the current validated project records, not newly re-imported field data. `data/README.md` documents removal of the fictitious Gottin HP V1 consumer and unresolved one-link differences in Garraf HP and Olot; neither is silently altered here.

Confirmed independent host-presence additions for Quercus, Gottin HP, Garraf HP and Olot were analysed separately in `outputs/support_prediction_zero`. Ambiguous additional records were retained as an inclusive sensitivity. Their incorporation into the primary occupancy dataset remains outstanding. Since neither version changes interaction triples or N, they leave this descriptive K analysis unchanged. Future interaction or site-universe corrections require rerunning it.

## Files and figure definitions

- `summary.csv`: N, regional link count, total support, mean, median, sample SD (denominator L-1), CV=SD/mean, 25th/75th percentiles, 90th percentile and singleton fraction. Median is the second quartile. Quantiles use Julia's default linear interpolation; noninteger quantiles are valid. Additional descriptors are maximum K, mean K/N, Gini and top-decile support share.
- `link_support.csv`: pair identities and their K, N and K/N.
- `distributions.csv`: exact frequencies and CCDFs at every integer K=1,...,N, including zero-frequency bins. Frequencies are normalised by regional link count, never total support. CCDF(K)=number of links with support >=K divided by regional link count.
- `concentration.csv`: ordered cumulative link fractions and cumulative support shares, including (0,0), ranking links from lowest to highest K. Ties do not affect the curve. The Gini is one minus twice its trapezoidal area; top-10% share uses linear interpolation at the 90th percentile of link rank to avoid rounding dataset-specific link counts.
- `01_support_frequencies.png`: exact frequency stems in dataset facets, with a shared y scale and dataset-specific raw-K ranges. Unoccupied bins have zero height.
- `02_support_ccdf_comparison.png`: all datasets on common raw-K and K/N axes. Raw K uses a logarithmic x axis to separate low-support values while retaining the full tail; probabilities remain linear. Steps represent the discrete distributions. No tail exponent is estimated.
- `03_support_concentration.png`: cumulative concentration curves, plus top-decile shares. Dashed references indicate equal support. Dataset colours are consistent across all figures.

`FINDINGS.md` gives the descriptive interpretation and sampling limitations. No significance tests or uncertainty intervals are implied by differences between curves.
