# Opportunity, realisation and local support

This analysis holds known links, site count, opportunities and observed total support fixed to compare allocation. Pair parameters are unchanged nine-dataset estimates. Common is the uniform eligible positive-row allocation conditional on M. Results concern whether opportunities alone plus a shared probability suffice, and which features persistent pair heterogeneity improves. They do not identify a biological mechanism.

The strongest recurring pattern is a departure from common-probability allocation: all ten empirical datasets exceed its central 95% intervals for singleton abundance, Gini and support in the highest-supported 10% of links. Opportunity heterogeneity and known total support therefore do not reproduce these features under this null. Persistent-pair heterogeneity moves concentration toward the observations in all ten datasets and singleton abundance in nine, but does not fully reproduce the distributions.

At equal opportunities, the picture is more specific: nine datasets exceed the Common interval for within-C realisation variance, yet the transferred Pair model overshoots the empirical variance (empirical below its interval) in seven. Pair is closer than Common on this scalar in only three datasets. Thus heterogeneous realisation is useful for the marginal distribution, but the existing two-state specification does not correctly distribute all that heterogeneity among opportunity classes. This is an important limitation that gap accuracy alone concealed.

Site count has a substantial effect on apparent shape: reducing Salix Galpar from 341 to eight sites raises the average singleton fraction from 45.5% to 93.8%, while about 92.7% of its original links disappear from the subsample. Other systems remain more broadly supported at eight sites. These results support recurring concentration relative to the Common null, not a sampling-invariant universal support distribution. The null comparison itself has not yet been repeated on subsamples.

## Which features need more than common realisation?

| Feature | Empirical above Common interval | Below Common interval | Inside Common interval | Inside Pair interval | Pair closer than Common |
|---|---:|---:|---:|---:|---:|
| singleton | 10 | 0 | 0 | 6 | 9 / 10 |
| gini | 10 | 0 | 0 | 2 | 10 / 10 |
| top10 | 10 | 0 | 0 | 5 | 10 / 10 |
| within_C_variance | 9 | 0 | 1 | 2 | 3 / 10 |

Intervals are central 95% conditional allocation intervals from 499 independent exact draws, not confidence intervals, removal-sequence variation or multiple-testing-corrected significance tests. Empirical positions are descriptive ranks. Mean Monte Carlo SE is exported.

Opportunity/support CCDFs compare the same known links. Common preserves each link's C and total M: its discrepancies cannot be attributed merely to differing opportunities or total support. Pair is a specific transferred two-state model; improvements support persistent heterogeneity as a useful statistical description, not proof of a two-class ecological mechanism. Observation differences could also contribute.

## Variation at equal opportunity

Each conditional point groups links with exactly the same C; figures show only groups with at least three links. Both empirical and model summaries use identical groups. The within-C scalar is the link-weighted mean squared deviation of K/C from its within-class mean, across eligible classes. The figure displays its square root in percentage points. This prevents variation in C alone from creating the measured within-class spread. Sparse groups and omitted classes are retained in conditional_C.csv, with link counts.

## Eight-site sampling check

| Dataset | Full singleton % | Eight-site mean singleton % | Full Gini % | Eight-site mean Gini % | Mean original links lost % |
|---|---:|---:|---:|---:|---:|
| Quercus | 49.624 | 54.970 | 48.998 | 33.356 | 39.885 |
| Nahuel | 55.215 | 55.215 | 36.237 | 36.237 | 0.000 |
| Salix Galpar | 45.474 | 93.790 | 50.424 | 5.493 | 92.740 |
| Gottin HP | 49.074 | 65.262 | 53.483 | 33.976 | 54.009 |
| Gottin PP | 58.387 | 71.724 | 48.196 | 27.887 | 58.505 |
| Garraf HP | 70.787 | 79.547 | 36.831 | 20.707 | 56.787 |
| Garraf PP | 63.692 | 78.970 | 47.653 | 26.763 | 68.278 |
| Garraf PP2 | 39.400 | 46.689 | 53.100 | 39.453 | 32.332 |
| Olot | 50.000 | 53.817 | 44.368 | 35.005 | 24.373 |
| Montseny | 67.044 | 72.047 | 42.225 | 30.190 | 42.094 |

Each replicate samples sites uniformly without replacement, preserving their interaction and occupancy records. K and C are recomputed. Positive-support summaries use links observed in that subsample; opportunity distributions use those exact same surviving links. Original known links disappearing from the sample are explicitly counted, including those with remaining co-occurrence. This separates apparent positive-distribution shape from loss of link identities. Three hundred permutations per dataset provide intervals at 8 sites and a sequence of other sample sizes; each dataset's full-site result reproduces the input. Nahuel already has eight sites and has no eight-site subsampling variation.

Matching site count tests one source of sampling dependence. It does not equalise within-site effort, detection, spatial extent, spatial clustering or seasonal coverage. Eight sites provide only eight possible support values and cannot establish a tail family. Pair/Common are evaluated on full datasets; they are not refitted or reconditioned in the subsampling exercise. Comparing subsampled positive links changes the cohort, so persistence of a qualitative shape is evidence of descriptive recurrence only.

## Files and limitations

All distributions and concentration curves, exact expected model frequencies, conditional group summaries, allocation draws and site-subsampling draws are exported. Gini and top-10% support share use fractional-rank interpolation; model concentration curves average curves calculated separately for each allocation, not curves of expected individual support. No power-law fitting, new predictive model or ecological mechanism test is added. Shared sites, related collections, interaction-derived occupancy and unresolved data flags are retained in README and metadata. All source versions are hashed and verified.
