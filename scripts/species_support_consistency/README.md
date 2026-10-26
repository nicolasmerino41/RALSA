# Do species realise similar support across partners?

Run from RALSA:

```text
julia --startup-file=no --project=. scripts/species_support_consistency/run_analysis.jl
```

Optional arguments: repository root, output directory, and `plot` to redraw saved results. Outputs default to `outputs/statistical description of local support/opportunity_theory/species_support_consistency`. Existing analyses are untouched; reruns replace this analysis's outputs only.

For every observed regional link, use validated K, C and N. Test consumers and resources separately. A species needs at least two realised partners to contribute. Species with one partner cannot inform across-partner consistency; they remain in the permutation pool. Analyse both K and q=K/C.

For each species, calculate the average squared difference between its partners' values. This equals twice the sample variance. Average over species with equal species weights, so highly connected species do not dominate by their number of partner pairs. Small values indicate consistency. Species-level means, pooled rates and partner counts are also exported.

Two permutation comparisons (1,999 permutations, fixed seeds):

1. Unrestricted: reassign observed (K,C) tuples to the existing regional links. This measures total grouping by species, including opportunity differences. It is a descriptive shuffle, not a feasible allocation preserving the original pair-specific opportunity matrix.
2. Primary, within C: reassign K only among links with identical C. Preserve regional identities, network degrees, each pair's C, total M, and the entire observed K distribution within every opportunity class. This controls both the mean and variability of q associated with C without assuming a regression or trend. Singletons with C=1 remain q=1. This test is conditional on observed regional links and does not investigate K=0 pairs.

Similarity effect = 100*(1 - observed dispersion / mean shuffled dispersion). A value of 20 means squared differences among a species's partners are 20% smaller than expected under that shuffle. This is not a correlation, proportion of variance explained, or percentage of species with similar partners.

One-sided permutation p=(1+number of null dispersions <= observed)/(2000). Adjust all 40 primary dataset/role/metric tests with Benjamini–Hochberg; adjust the 40 unrestricted tests separately. Non-significance is not proof of no species effect. Exact-C strata with no variation provide no randomisation information: eligible species, eligible links and movable link fractions are reported. The movable fraction ranges approximately 55–83% across these datasets.

Figures show observed effects and central 95% null intervals scaled relative to the same null mean. Intervals are null reference ranges, not confidence intervals around observed effects. Circles identify adjusted p<0.05; crosses identify other results. Both figures were inspected for legibility.

These tests are marginal, not a crossed consumer-resource model: a role-level signal may include partner composition and shared environments or observation effort. No causal abundance/compatibility attribution follows from these tests. Absence of a detectable species-wide signal does not prove a pair-specific mechanism.

Files: 01_opportunity_controlled.png, 02_unadjusted.png, consistency_tests.csv, species_summaries.csv, checks.csv, provenance/relationship tables, and FINDINGS.md. Provenance and pre-existing occurrence limitations are retained in the findings note.
