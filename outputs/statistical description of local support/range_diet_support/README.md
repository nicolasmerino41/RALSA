# Range, partner breadth and support

Run `julia --startup-file=no --project=. scripts/range_diet_support/run_analysis.jl` from RALSA. Optional arguments are repository root and output folder. Outputs default to `outputs/statistical description of local support/range_diet_support`. Existing analyses are preserved; rerunning overwrites only this folder's results.

1. Species-level degree versus site occupancy, separately for consumers and resources. Points are species; lines are medians at each degree, without fitted regressions. Degree-zero species remain in the exported table but are not treated as diet specialists.
2. Singleton-link frequency by endpoint degree class. Every link occurs once per role; consumer and resource views are not additive.
3. Species share versus link share versus singleton-link share for species with at most two partners. Thresholds one and three are also exported.
4. Exact structural decomposition for K<=h, h=1,2,3: endpoint occupancy<=h; otherwise co-occurrence C<=h; otherwise K<=h despite C>h. Categories partition each rare-link set. They identify restrictions, not causal contributions.
5. Conditional tests of degree–K association and singleton-rate contrast (degree<=2 minus degree>2), shuffling K within exact C classes, 1,999 replicates. All 40 two-sided Monte Carlo comparisons receive BH correction. Intervals are null reference ranges; circles mark adjusted p<0.05, crosses other cases. Movable fractions are exported. These are conditional exchangeability diagnostics; no independent-link or causal inference is asserted.
6. Descriptive rank correlations linking degree, occupancy, mean K and mean C, including partial rank correlations. No statistical significance or causal mediation is attached to these coefficients.

All regional degrees use observed links K>0. Species names retain trophic roles. Resource degree describes visitor/enemy breadth, not diet. Species occupancy is within the sampled site universe, not geographic range area. Current K and C are verified directly from confirmed occurrences and interactions. See FINDINGS.md and provenance tables for data-quality and collection limitations.
