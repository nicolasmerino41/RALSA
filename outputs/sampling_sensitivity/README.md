# Sampling sensitivity

Run `julia --project=. scripts/sampling_sensitivity.jl` from RALSA. Use
`--plot-only` to redraw existing results. Main manuscript outputs are not modified.

## Design

For each dataset, draw 300 random site permutations (seed 90100 plus dataset index).
Retain nested samples of 5, 8, 14, 20, 30, 50, 100 and 200 sites where available,
plus the full sample. Full samples are evaluated once because their results are exact.
All ten datasets can be compared at 5 or 8 sites; larger sample sizes include only
eligible datasets. Comparisons use the same eligible datasets at full and reduced sizes.

Rebuild each starting network from its sampled sites. Links not observed in that
sample are absent from its initial link inventory. Recompute interaction support,
co-occurrence support and empirical realisation, including all co-occurring pair–site
opportunities in the latter denominator. Further site removal operates within that
sampled network. Its outcomes are calculated analytically, without additional
Monte Carlo removal error.

For target fraction h and n sampled sites, mix exact expectations for floor(h*n)
and ceil(h*n) removed sites with weights giving expected removal h*n. This avoids
different realised removal fractions driving comparisons at small n. Full-sample
results use the same convention. The original manuscript's rounded-removal pink
areas are separately reproduced and checked against its tables before analysis.

Mean gap and mean co-occurrence-only fraction are trapezoidal areas over 0–80%
removal, divided by 0.8. The former is regional interaction retention minus local
support retention. The latter is interaction loss probability minus co-occurrence
loss probability, averaged over initially realised links.

Species contrasts use both refitted groups (the manuscript's split of distinct
degree values within each sampled dataset/guild) and fixed full-network groups.
Only species with interactions in the sampled starting network enter either
comparison. Guilds are averaged equally. A comparison requires both groups in
both guilds; otherwise it is recorded as missing (NaN). Positive differences mean
larger values for generalists. Persistence and absolute link loss are evaluated at
80% expected further removal. Species activity is computed from the union of sites
supporting any of its links, rather than assuming independent partner losses.

## Files

- `summary.csv`: sample-size means and 2.5–97.5% rarefaction intervals.
- `comparisons.csv`: matched-dataset full versus standardised means, ranking
  correlations, between-dataset variation and descriptive p_emp correlations.
- `replicates.csv`: all network metrics per subsample.
- `group_summary.csv`, `group_replicates.csv`: generalist–specialist contrasts.
- `correlation_replicates.csv`: p_emp–pink correlations across rarefaction draws.
- `validation.csv`: verification of original full-sample pink results.
- `sampling_effects.png/pdf`: within-dataset trajectories and common-size comparisons.
- `sampling_species_groups.png/pdf`: group contrasts with refitted classifications.

Intervals describe variation among subsets of the observed sites; they are not
confidence intervals for true landscape values. Nested sample sizes and replicates
from the same dataset are not independent ecological datasets. Full-sample intervals
collapse because all observed sites are included.

## Interpretation

This evaluates sensitivity to the number and selection of sampled sites. It does
not equalise within-site effort, detectability, sampled area, geographical extent,
or environmental coverage. Downsampling changes the observed regional species and
link pool as well as their support. A difference between full and rarefied results
cannot by itself be assigned entirely to observation error or entirely to ecology.
Likewise, a remaining contrast at equal site counts is not proof of an ecological
mechanism. Generalisation classes are relative degree classes, not independent
measurements of niche breadth.
