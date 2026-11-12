# All-co-occurrence null, free regional link richness

Each allocation selects exactly M distinct pair-site opportunity cells uniformly without replacement from all positive co-occurring pairs. Regional link identities and richness are free. No parameter fitting. 5,000 independent draws per dataset, fixed seeds. Cell selection is uniform, not uniform over support-count vectors: vectors receive probability proportional to the product of binomial(C,K).

K=0 fraction uses all positive-C pairs. Singleton fractions and support distributions use only K>0 links within each allocation. Gap uses that allocation's own original realised link count, with empirical and null local retention 1-m/N. Mean gaps average exact expected retention over integer removal levels m/N≤0.8, including zero. Thus gap differences include changed regional richness and allocation, not allocation alone.

| Dataset | Regional links empirical / null | Zero pairs % empirical / null | Singleton % empirical / null | Gap pp empirical / null [95%] |
|---|---:|---:|---:|---:|
| Quercus | 133 / 321.7 | 88.8 / 73.0 | 49.6 / 74.7 | 15.40 / 5.46 [4.56, 6.42] |
| Nahuel | 163 / 277.7 | 82.6 / 70.4 | 55.2 / 83.0 | 12.18 / 3.62 [2.82, 4.53] |
| Salix_Galpar | 939 / 1662.9 | 70.4 / 47.5 | 45.5 / 60.6 | 15.88 / 9.91 [9.46, 10.36] |
| Gottin_HP | 108 / 250.0 | 87.1 / 70.1 | 49.1 / 71.6 | 14.51 / 6.38 [5.35, 7.50] |
| Gottin_PP | 781 / 1492.0 | 87.9 / 76.8 | 58.4 / 77.5 | 11.54 / 5.01 [4.61, 5.43] |
| Garraf_HP | 89 / 137.4 | 85.7 / 77.9 | 70.8 / 88.1 | 7.62 / 2.54 [1.61, 3.59] |
| Garraf_PP | 325 / 503.7 | 76.0 / 62.8 | 63.7 / 77.2 | 9.45 / 5.61 [4.96, 6.28] |
| Garraf_PP2 | 967 / 2148.8 | 73.7 / 41.5 | 39.4 / 56.6 | 18.66 / 11.23 [10.84, 11.63] |
| Olot | 92 / 186.7 | 83.5 / 66.5 | 50.0 / 74.5 | 14.89 / 5.66 [4.49, 6.98] |
| Montseny | 707 / 1143.9 | 81.5 / 70.1 | 67.0 / 81.5 | 9.04 / 4.12 [3.72, 4.53] |

regional_links: empirical below / inside / above allocation interval: 10 / 0 / 0 datasets.

zero_pct: empirical below / inside / above allocation interval: 0 / 0 / 10 datasets.

singleton_pct: empirical below / inside / above allocation interval: 10 / 0 / 0 datasets.

gap_pp: empirical below / inside / above allocation interval: 0 / 0 / 10 datasets.

Across all ten systems, empirical interactions occupy fewer pairs and repeat more often within the realised subset than equal random allocation across all opportunities predicts. Both the lower singleton fraction and higher gap lie outside the null central 95% intervals in every dataset. This differs from a known-link null because regional richness is no longer fixed. At fixed M, fewer realised links necessarily raise mean support M/L; the gap comparison cannot isolate allocation shape from this change in mean support. A future constrained-L test would separate those contributions.


These are conditional allocation intervals, not confidence intervals or removal-sequence variation. They test equal allocation among recorded opportunities, not ecological neutrality generally. Departures do not identify abundance, compatibility or detection mechanisms. No unseen species/pairs are estimated. Restored confirmed occurrences, interaction-derived presences, incomplete detection, sampling differences and shared collections retain their existing limitations (see copied audit and relationship tables). All opportunities were reconstructed from confirmed occupancy; known K was reconstructed from confirmed interaction records. Simulation mean link richness is checked against its exact hypergeometric expectation. Existing analyses are preserved.

Run: julia --startup-file=no --project=RALSA_PATH "PATH_TO_THIS_FOLDER/run.jl".
