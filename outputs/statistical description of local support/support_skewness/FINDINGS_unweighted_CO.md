# Unweighted species occupancy to all co-occurrences to support

All recorded species from both guilds receive equal weight. All CO includes every pair co-occurring at least once, including pairs never observed interacting. Each pair receives equal weight; histogram counts only recover the underlying observations. Local support includes realised regional links (K>0). CO and support therefore describe different pair sets. Skewness is the third central moment divided by variance to the power 3/2, with no sample-size correction. Species and support are verified against the existing summary; all-CO skewness is calculated from its saved frequency counts.

Species to CO increases in 9/10 datasets. CO to support increases in 7/10. Both steps increase in 6/10.

| Dataset | Species unweighted | CO | Support | Both increase |
|---|---:|---:|---:|---|
| Quercus | 0.714 | 1.534 | 2.018 | true |
| Nahuel | 0.916 | 1.506 | 1.667 | true |
| Salix_Galpar | 4.844 | 3.188 | 3.275 | false |
| Gottin_HP | 1.055 | 2.428 | 2.991 | true |
| Gottin_PP | 1.696 | 3.138 | 3.761 | true |
| Garraf_HP | 2.359 | 4.893 | 3.890 | false |
| Garraf_PP | 3.780 | 5.842 | 6.077 | true |
| Garraf_PP2 | 0.855 | 2.237 | 1.751 | false |
| Olot | 0.800 | 1.856 | 1.467 | false |
| Montseny | 1.690 | 3.374 | 3.526 | true |

Descriptive comparisons, not significance tests or a causal assembly test. Species and pairs are different observational units; skewness need not increase when per-pair counts decrease. Existing sampling, detection, occurrence and shared-collection limitations apply. Weighted plots and all previous analyses are preserved.
