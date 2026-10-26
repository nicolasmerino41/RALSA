# Species consistency of local support

Primary test: shuffle observed K only among links with exactly the same C, retaining pair identities and network topology. Each eligible species has at least two partners and receives equal weight. Statistic: average squared difference between two partners' q or K, averaged equally over species. Positive effect = smaller differences than the null; 100 × (1 − observed dispersion / mean null dispersion). 1,999 permutations; one-sided p for similarity; Benjamini–Hochberg correction across all 40 primary tests. Unrestricted tests have a separate 40-test correction.

Plot segments show central 95% null intervals centred relative to the null mean; they are not confidence intervals around observed effects. Crosses denote lack of adjusted evidence for similarity, not evidence of no structure. K=1,C=1 links have q=1 by selection, and are handled by exact-C conditioning. Only observed regional links are tested.

- Within_C, consumer, q: 0/10 datasets show adjusted evidence for similarity. Median effect 6.7%. Significant datasets: none.
- Within_C, resource, q: 1/10 datasets show adjusted evidence for similarity. Median effect -4.0%. Significant datasets: Garraf_PP2.
- Within_C, consumer, K: 0/10 datasets show adjusted evidence for similarity. Median effect -13.5%. Significant datasets: none.
- Within_C, resource, K: 3/10 datasets show adjusted evidence for similarity. Median effect 3.1%. Significant datasets: Garraf_PP2, Olot, Montseny.
- Unrestricted, consumer, q: 6/10 datasets show adjusted evidence for similarity. Median effect 28.0%. Significant datasets: Quercus, Salix_Galpar, Gottin_PP, Garraf_HP, Garraf_PP2, Montseny.
- Unrestricted, resource, q: 5/10 datasets show adjusted evidence for similarity. Median effect 18.1%. Significant datasets: Salix_Galpar, Gottin_PP, Garraf_PP, Garraf_PP2, Montseny.
- Unrestricted, consumer, K: 4/10 datasets show adjusted evidence for similarity. Median effect 39.2%. Significant datasets: Salix_Galpar, Gottin_PP, Garraf_PP2, Montseny.
- Unrestricted, resource, K: 6/10 datasets show adjusted evidence for similarity. Median effect 49.7%. Significant datasets: Salix_Galpar, Gottin_PP, Garraf_PP, Garraf_PP2, Olot, Montseny.

The unrestricted comparison moves (K,C) together and asks whether species group similar observed supports/rates overall. Exact-C shuffling asks whether grouping exceeds what each species' opportunity profile already implies, while preserving the empirical support distribution at every C and total M. Sparse C classes cannot be shuffled meaningfully; movable fractions are reported. These are two marginal role tests, not simultaneous consumer/resource variance components: a species signal can reflect its partners, shared sites, observation effort or other network structure. It is not causal evidence for abundance or compatibility.

Confirmed input versions and existing occurrence flags are retained. Uncertain excluded records: Gottin PP 474 missing-identity records, Garraf PP2 four zero columns, Montseny one zero row. Much occupancy is interaction-derived. Prior Gottin HP V1 correction, Garraf HP/Olot one-link discrepancies and undocumented sites remain limitations. Gottin datasets share sites and other collection overlaps remain unresolved; counts of significant datasets are descriptive, not independent replication probabilities.