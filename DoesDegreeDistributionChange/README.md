# DoesDegreeDistributionChange

Run from this folder or any directory:

`julia --startup-file=no --project=PATH_TO_RALSA run_analysis.jl PATH_TO_RALSA`

Optional second argument: output folder. Optional third argument: paired replicate count (default 2000). Outputs are restricted to this folder's `outputs/`; existing analyses are preserved. Read `outputs/START_HERE.md` first.

The primary statistic is change in Gini of positive realised degrees. Exact expected proportional degree retention tests the proposed generalist advantage directly. Mean-scaled cumulative distributions assess shape without confusing it with mean-degree loss. Species that disconnect are reported separately.

The null uniformly permutes complete site-support columns across existing regional links. It preserves site interaction totals, K distribution and the regional graph, while breaking assignment to endpoints. A single shuffled assignment and one random site ordering are drawn per replicate, shared across 40% and 80% removals. Endpoint occupancy feasibility is deliberately not preserved. This isolates support assignment collectively, not a pure causal effect of support–degree correlation. Both guilds use the same link permutation.

Intervals describe removal/allocation variation or Monte Carlo precision, as labelled; they are not biological confidence intervals. No hypothesis testing across dependent links or claims of ten independent landscapes are made. All source hashes and outstanding data issues are retained.
