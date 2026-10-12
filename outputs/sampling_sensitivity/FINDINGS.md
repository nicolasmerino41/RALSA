# Sampling sensitivity: findings

300 nested random subsamples per dataset, with exact expected subsequent removal.
All ten systems compared at 5 and 8 sites. Full-sample outcomes are calculated with
the same removal-fraction convention as reduced samples (see README.md).

| Outcome, averaged over 0–80% removal | Full samples | 8 sites | 5 sites |
| --- | ---: | ---: | ---: |
| Regional–local retention gap | 13.20 percentage points | 9.26 | 7.81 |
| Co-occurrence-only fraction | 9.46% | 6.77% | 5.88% |

At eight sites the mean gap is about 30% smaller and the co-occurrence-only fraction
about 28% smaller than in full samples. These are changes under rarefaction, not
estimates of the fraction caused by observation bias.

## Which claims persist?

**Magnitude and ranking are sampling-sensitive.** The gap's dataset-rank correlation
between full samples and eight-site means is only 0.47 (Spearman). Salix–Galpar shifts
from a 15.82-percentage-point mean gap at 341 sites to 1.29 points at eight sites, a
92% reduction. Its co-occurrence-only fraction shifts from 9.98% to 1.22%. This system
has little repeated support when only a few of its observed sites are included.

**Variation remains at equal site counts.** Eight-site mean gaps range from 1.29 to
16.22 percentage points; mean co-occurrence-only fractions range from 0.86% to 13.72%.
Site count alone does not erase differences in the observed networks, but equal site
counts do not equalise geographical extent, environmental coverage or detection.

**The empirical realisation–pink association persists.** Across the ten eight-site
dataset means Pearson r is -0.943. Across 300 independently paired rarefaction draws,
median r is approximately -0.87 and the central 95% range approximately -0.95 to -0.71.
This remains descriptive: the metrics share information about interaction realisation.

**Generalist contrasts persist in direction, but shrink in magnitude.** At eight sites,
all ten dataset means show higher support per partner, more absolute partner loss,
and greater network participation for generalists when groups are refitted. Fixed
full-network groups also retain positive mean support and participation differences
in all ten datasets. This does not imply every subsample shows the same direction:
Salix–Galpar's support difference is positive in only 69% of eight-site subsamples
with refitted groups, versus 91% or more in the other datasets (Nahuel's full eight
sites form a single deterministic sample).

## Implication for the manuscript

The regional–local gap is a property of the observed support distribution, and its
cross-system magnitude cannot be interpreted independently of sampling scale.
The stronger findings after equalising site counts are the persistent differences
among observed systems, the co-occurrence/realisation relationship, and the mean
generalist–specialist contrasts. Site-number standardisation should accompany
cross-system comparisons. This analysis has not tested unequal within-site effort
or interaction detectability and cannot establish that remaining differences are
entirely ecological.
