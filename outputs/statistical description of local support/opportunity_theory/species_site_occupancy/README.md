# Species site occupancy frequencies

Confirmed species-occurrence records, including recorded presences without an observed interaction. Each species' occupancy is the number of distinct sampled sites where it is recorded. N comes from the established site universe. No undocumented or wholly unrecorded species are inferred. Y is the proportion of species in each group. Absolute and percentage plots show identical discrete frequencies: percentage=100*n/N, with no rebinning. Percentage plots span the full occupancy range; absolute ranges adapt to observed occupancy. Y scales match across datasets and both versions within each group. Combined panels pool species with equal species weight, retaining trophic-role identity (not an equal-weight average of consumer/resource distributions).

These describe site occupancy, not partner specialism. Sampling effort, detectability, geographic extent and interaction-derived occurrence records affect interpretation. Confirmed restorations and unresolved occurrence flags are unchanged: Gottin PP missing-identity records, Garraf PP2 zero columns and Montseny zero row remain excluded as before. Shared sites/collections limit independence.

Reproduce with: julia --startup-file=no --project=. scripts/species_site_occupancy.jl REPO OUTPUT_FOLDER
