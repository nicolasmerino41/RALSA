# What distribution describes local support?

K is the number of distinct sites supporting a known regional link (K≥1). These are fits to the empirical support distribution, not degree distributions, co-occurrence models, or gap predictions. All candidate laws are discrete and normalised on 1,...,N, including unobserved support levels.

A single exponential shows detectable mismatch in all ten datasets. A bounded power law is compatible with eight, but it is not uniquely supported: several alternatives have similar AICc. Salix Galpar clearly favours an exponential cutoff, Quercus favours the singleton-plus-geometric model, and no tested family adequately describes Garraf PP2. These results support recurrent broad distributions with dataset-specific details, rather than one established universal law.

| Dataset | Lowest-AICc family | KS bootstrap p for winner | Other families within 2 AICc units |
|---|---|---:|---|
| Quercus | Singleton geometric | 0.416 | — |
| Nahuel | Power law | 0.826 | Cutoff power law, Discrete lognormal, Singleton geometric |
| Salix Galpar | Cutoff power law | 0.532 | — |
| Gottin HP | Power law | 0.288 | Cutoff power law, Discrete lognormal |
| Gottin PP | Cutoff power law | 0.516 | Discrete lognormal |
| Garraf HP | Power law | 0.370 | Cutoff power law, Discrete lognormal |
| Garraf PP | Power law | 0.284 | Discrete lognormal |
| Garraf PP2 | Power law | 0.022 | — |
| Olot | Singleton geometric | 0.644 | Power law |
| Montseny | Power law | 0.438 | — |

AICc ranks candidates relatively. A low bootstrap p flags systematic distributional mismatch even if a model has the lowest AICc. Conversely, p≥0.05 means this check did not detect mismatch; it does not establish a true family or universality. Values near 0.05 should be read with their Monte Carlo errors.

| Family | Datasets with bootstrap p ≥ 0.05 | Lowest-AICc datasets |
|---|---:|---:|
| Exponential | 0/10 | 0/10 |
| Power law | 8/10 | 6/10 |
| Cutoff power law | 8/10 | 2/10 |
| Discrete lognormal | 8/10 | 0/10 |
| Singleton geometric | 4/10 | 2/10 |

## Tail-only check

The lower cutoff is selected by minimum power-law KS distance among cutoffs Kmin≥2 retaining at least 50 links and five distinct observed support values. The power-law bootstrap regenerates the tail, resamples the empirical body, and reselects the cutoff. Alternatives are compared on that same selected tail. Tail AICc values must not be compared with whole-distribution AICc.

| Dataset | Kmin | Tail links | Observed Kmax/Kmin | Power-law bootstrap p |
|---|---:|---:|---:|---:|
| Quercus | 2 | 67 | 9.5 | 0.010 |
| Nahuel | 2 | 73 | 4.0 | 0.832 |
| Salix Galpar | 8 | 99 | 4.0 | 0.072 |
| Gottin HP | 2 | 55 | 13.0 | 0.435 |
| Gottin PP | 2 | 325 | 15.0 | 0.256 |
| Garraf HP | Insufficient data | — | — | — |
| Garraf PP | 2 | 118 | 19.0 | 0.458 |
| Garraf PP2 | 2 | 586 | 10.5 | 0.006 |
| Olot | Insufficient data | — | — | — |
| Montseny | 2 | 233 | 9.0 | 0.308 |

Any compatible tail is a limited finite-range description. K is bounded by N and observed tails cover a restricted range; these data cannot establish an unbounded, scale-free asymptotic law. A cutoff power law here adds an exponential decay below the hard upper bound N; the hard N bound applies to every model.

Tail bootstrap samples without an eligible tail are excluded from the conditional bootstrap calculation. Gottin HP retains 420 of 499 samples and Quercus 498; other eligible datasets retain all 499. This selection particularly limits interpretation for Gottin HP. Valid replicate counts are recorded in tail_selection.csv.

## Interpretation limits

The tests treat regional-link support counts as independent draws. Links share species and sites, so bootstrap p-values and AICc are working-model diagnostics, not independence-robust ecological tests. Neither counts of non-rejections nor common rank winners establish universality. Site counts, sampling effort, spatial extent and observation can affect shape; K/N alone would not eliminate those effects. Shared campaigns further reduce independent replication.

Existing geometric and singleton-geometric likelihoods are reproduced. Input interaction triples and site counts are checked against the descriptive analysis. Occurrence restorations do not alter these K distributions because neither interactions nor N changed; outstanding occurrence and interaction flags remain in README. No original analysis is overwritten.

Methodological references: Clauset, Shalizi & Newman (2009), https://doi.org/10.1137/070710111; Broido & Clauset (2019), https://doi.org/10.1038/s41467-019-08746-5. Our finite-support, discrete bootstrap is implemented for this bounded site-support problem, not copied as an unbounded power-law test.
