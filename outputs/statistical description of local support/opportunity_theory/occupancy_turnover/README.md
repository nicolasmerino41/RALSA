# Occupancy, singleton origins and network dissimilarity

Run from RALSA:

```text
julia --startup-file=no --project=. scripts/occupancy_turnover/run_analysis.jl
```

Optional arguments: repository root, output directory, then `plot` to redraw without recalculation. Uses existing Julia/CairoMakie environment; modifies only this analysis's outputs.

## Inputs and comparisons

Use confirmed occurrence tables, processed interactions and the established site universe. Verify original input hashes and reconstruct every K and C to match quality_control.csv. Deduplicate interaction-site records. Species identities retain trophic roles.

1. Compare singleton frequency among all recorded species, species participating in at least one regional interaction, degree-weighted species, known-link co-occurrence C, and realised support K. Species frequencies have a different denominator from pair frequencies. High–low–high means species singleton fraction > C=1 fraction and K=1 fraction > C=1 fraction. This is a descriptive comparison, not a fitted statistical test. Degree-weighted occupancy repeats each species once for each regional partner. Combined weights total 2L. Full occupancy frequency distributions and CCDF plots are included.
2. Classify each singleton interaction: at least one endpoint occupies one site; both occupy more than one site but intersect at only one; or they intersect at more than one site but interact at only one. These mutually exclusive classes exhaust K=1. The first two sum exactly to C=1.
3. Compare every unordered pair of sites using the qualitative framework of Poisot et al. (2012), The dissimilarity of species interaction networks, doi:10.1111/ele.12002. Species presence includes recorded species without an interaction at that site. βS is Sørensen species dissimilarity. Let a be shared interactions, r be differing interactions among species present at both sites, and t be differing interactions requiring at least one species absent from the other site. Compute βWN=(r+t)/(2a+r+t), βOS=r/(2a+r), βST=βWN−βOS.

The βOS denominator restricts to shared-species networks. Thus βST is a subtraction of separately normalised dissimilarities, not a literal fraction of differing interactions. For clarity, also report t/(r+t) pooled over all pairwise differing-link occurrences. The same regional link can contribute to many site comparisons; this is not a percentage of unique regional links. Avoid interpreting these as independent samples or causal fractions.

If 2a+r=0, βOS is undefined. The all-pair figure explicitly uses zero as a convention; dataset_summary.csv also reports the number affected and means excluding these comparisons. In Salix Galpar 89% of site pairs are affected, so its low all-pair mean βOS is not evidence of universally consistent interactions among shared species. Two entirely interaction-empty sites give undefined βWN and are excluded from interaction means. Original N and site identities remain unchanged.

## Files

- 01_singleton_stages.png: recorded species → degree-weighted species → known-link CO → realised support. CO and support curves coincide across roles because they refer to the same regional pairs.
- 02_singleton_origins.png: three structural sources of singleton support.
- 03_network_dissimilarity.png: original Poisot decomposition and the separate direct classification of differing links.
- 04_degree_weighted_occupancy.png: unweighted, linked-only and degree-weighted combined species CCDFs.
- singleton_stages.csv, species.csv, occupancy_distributions.csv: species-weighting comparisons.
- link_categories.csv: every known regional interaction and its category.
- site_pair_dissimilarity.csv: complete pairwise values and undefined flags.
- dataset_summary.csv, FINDINGS.md, source/input versions and collection metadata.

Validation includes small analytically known decomposition cases, βWN=βOS+βST, occupancy/interaction compatibility, K/C agreement with authoritative records, probability normalisation and exhaustive singleton categories. Figures are visually inspected. No new models are fitted.

## Limits

Low occupancy is observed site rarity, not necessarily geographical-range restriction or habitat specialism. Interaction differences may reflect sampling, detectability, timing or ecological differences, not necessarily behavioural rewiring. Confirmed occurrence corrections and all existing data-quality limitations remain; see FINDINGS.md. Spatially dependent pairs and shared campaigns are not independent replication. The analysis distinguishes structural routes to limited support, without identifying abundance or compatibility as a cause.
