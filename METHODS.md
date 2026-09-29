# Analysis and figure definitions

## Site removal

Sites are removed uniformly without replacement. At nominal fraction `f`, the
retained count is `max(1, round(N*(1-f)))`. Rounding means realised removal fractions
can differ slightly among datasets. Curves use nominal fractions on the horizontal axis.
For support `k`, exact regional loss probability is `choose(m,k)/choose(N,k)` when
`m >= k`, otherwise zero. Species support uses guild-specific occupancy.
This represents sampled-site removal, without geographical clustering, rewiring,
abundance changes or subsequent population dynamics.

## Figure 2: interaction support

(a) Expected proportions of regional species, regional links and local pair–site
interaction occurrences retained. Thin curves show datasets; thick curves give
equal-weight dataset means. (b) Supporting-site counts for realised regional links.
Points are links; boxes show medians and interquartile ranges, with 1.5-IQR whiskers.
The vertical axis is logarithmic.

## Figure 3: co-occurrence and interaction loss

(a) Expected states of originally realised regional links: interaction retained,
co-occurrence retained without a recorded interaction, and co-occurrence lost.
These fractions sum to one. For each link, co-occurrence-only probability is
`loss(K) - loss(n)`, where `K` is interaction support and `n` co-occurrence support.
(b) Cross-dataset relationships between the mean co-occurrence-only fraction and
empirical realisation (left) or mean co-occurrence support (right). The response is
trapezoidal area over 0–80% removal divided by 0.8. Empirical realisation is all
realised pair–site occurrences divided by all co-occurring pair–site opportunities,
including pairs never recorded interacting. Mean co-occurrence support is restricted
to originally realised regional links. Lines are descriptive OLS fits; r is Pearson
correlation across ten datasets. These relationships are not independent causal effects.

## Figure 4: generalists and specialists

Groups are fixed from initial regional degree, separately within each dataset and
guild. Specialists occupy the lower half of **distinct degree values**, generalists
the upper half. This is not a median split of species. Exact group sizes and degree
ranges are saved in `group_definitions.csv`.

(a) Supporting-site count averaged over realised partners, then species within
groups and guilds. Guild means receive equal weight within each dataset, and datasets
receive equal weight overall. Faint points and lines show datasets; large points are
dataset means and bars are across-dataset interquartile ranges.
(b) Mean absolute partners lost per original species, including disconnected species
(upper), and fraction of original species retaining at least one interaction (lower).
Thin lines show dataset means and thick lines equal-weight means across datasets.
Each dataset uses 500 nested random-removal permutations, seeded with 4300 plus its
position in `DATA`. Network participation is not independent demographic survival.
(c) Salix–Galpar consumer complementary cumulative degree distributions among active
consumers at 0%, 40% and 80% removal. Curves are replicate medians of `Pr(degree >= k)`;
absent high degrees contribute zeros. Only positive change points are connected,
matching the established plotting convention. This panel is an illustrative example.

## Scope

Input definitions and unresolved source discrepancies are documented in `data/README.md`.
