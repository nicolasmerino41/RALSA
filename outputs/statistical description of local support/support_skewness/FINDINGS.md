# Skewness across organisational stages

Species occupancy is weighted by initial regional degree, pooling both guilds. CO and K follow the same known regional pairs. Confirmed occurrence-only species have zero degree weight. Unweighted species results are retained separately in the summary CSV.

Skewness = weighted third central moment / weighted variance^(3/2). Weights sum to one; this is descriptive population-moment skewness, without sample-size correction. Positive values indicate a right tail. No prior scaling is required. Zero variance would return NaN.

Species → CO increases: 9/10. CO → support increases: 10/10. Both steps increase: 9/10.

| Dataset | Species weighted | CO | Support | Both increase |
|---|---:|---:|---:|---|
| Quercus | -0.452 | 0.566 | 2.018 | true |
| Nahuel | -0.445 | 0.663 | 1.667 | true |
| Salix_Galpar | 2.251 | 2.053 | 3.275 | false |
| Gottin_HP | -0.045 | 0.894 | 2.991 | true |
| Gottin_PP | 0.129 | 1.199 | 3.761 | true |
| Garraf_HP | 0.744 | 2.967 | 3.890 | true |
| Garraf_PP | 0.528 | 3.337 | 6.077 | true |
| Garraf_PP2 | 0.137 | 1.594 | 1.751 | true |
| Olot | -0.427 | 0.418 | 1.467 | true |
| Montseny | 0.223 | 1.452 | 3.526 | true |

This is a direct descriptive comparison, not a significance test or evidence of causal assembly. Shared species and sites prevent treating observations as independent replicates. Skewness is sensitive to rare high values; singleton fraction and mean remain in the CSV for context. Finite site count, sampling effort and detection affect the distributions. Degree weighting uses interaction information. Existing confirmed-occurrence exclusions and shared-collection limitations remain. Reconstructed C, K and unweighted occupancy distributions must exactly match the saved family-fitting observations before results are written. Other analyses are preserved.
