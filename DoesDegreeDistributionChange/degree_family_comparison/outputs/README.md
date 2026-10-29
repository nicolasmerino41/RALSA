# Degree family comparison — consumers

Compact table only. Rows are datasets; columns are 0%, 40% and 80% target site removal. Connected species only. Candidates are the previously used bounded exponential, power law, truncated power law, lognormal and singleton-plus-geometric; bounds are 1 through the initial number of connected resource species.

Zero removal is fitted once. At 40% and 80%, the first 100 of the original 2,000 removal outcomes are fitted separately. Winner: lowest mean within-replicate ΔAIC (equivalent to lowest mean AIC), not highest R². A second name appears only for exact numerical equality in that score, with at most two names. R² and KS shown are medians across those fits. KS is maximum CDF discrepancy, smaller is better; it is not a p-value. R² uses CCDF probabilities at all integer degrees in ordinary probability space. No significance tests are claimed for dependent species degrees.

Figures 6–9 are overwritten using all 2,000 original draws: linear x = raw number of links; logarithmic y = P(species degree ≥ x). No mean scaling. Probability is first calculated exactly at integer degrees within each replicate, then averaged with equal replicate weight. Straight segments connect distinct cumulative values; zero probabilities are omitted on the logarithmic axis. Original network, observed removal and shuffled support remain the three curves. Random seed, site sets and fixed-count removal are unchanged. Actual removal rounds to the nearest integer site count.

Input hashes are checked against ../outputs/source_versions.csv; all prior data-quality flags remain applicable. This script reuses the existing fitting helper in RALSA/scripts/K_CO_SP_distributions/distribution_helpers.jl. Run: julia --project=PATH_TO_RALSA run.jl PATH_TO_RALSA.
