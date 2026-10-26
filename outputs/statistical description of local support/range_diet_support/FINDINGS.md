# Species occupancy, partner breadth and local support

Occupancy means number of sampled sites, not geographic range area. Degree is observed regional partner breadth: diet/host breadth for consumers, enemy/visitor breadth for resources. Degree-zero species are listed but excluded from specialist comparisons. Degree and occupancy both depend on observation effort.

- consumer: positive degree–occupancy association in 10/10 datasets; median Spearman rho 0.79 (range 0.69–0.96).
  Species with <=2 partners are a majority in 8/10; median shares: 55.8% of linked species, 19.4% of regional links, 22.3% of singleton links. Equal-dataset summaries; no pooling.
- resource: positive degree–occupancy association in 10/10 datasets; median Spearman rho 0.73 (range 0.41–0.93).
  Species with <=2 partners are a majority in 3/10; median shares: 38.9% of linked species, 5.3% of regional links, 8.1% of singleton links. Equal-dataset summaries; no pooling.

## Dataset results

|Dataset|Consumer rho|Resource rho|K=1: one-site endpoint|K=1: one-site overlap only|K=1 despite C>1|
|---|---:|---:|---:|---:|---:|
|Quercus|0.70|0.66|28.8%|10.6%|60.6%|
|Nahuel|0.79|0.67|56.7%|5.6%|37.8%|
|Salix_Galpar|0.96|0.60|11.0%|32.1%|56.9%|
|Gottin_HP|0.79|0.70|9.4%|5.7%|84.9%|
|Gottin_PP|0.91|0.88|21.3%|7.7%|71.1%|
|Garraf_HP|0.75|0.85|33.3%|19.0%|47.6%|
|Garraf_PP|0.87|0.93|54.1%|9.7%|36.2%|
|Garraf_PP2|0.69|0.41|63.0%|21.5%|15.5%|
|Olot|0.74|0.75|28.3%|15.2%|56.5%|
|Montseny|0.88|0.81|37.3%|12.9%|49.8%|

## Beyond opportunities

Consumers with <=2 partners have a smaller singleton-rate excess than the within-C null in 10/10 datasets (median difference -12.0 percentage points). Thus opportunity restriction alone predicts a stronger specialist disadvantage than observed; the conditional pattern offsets rather than amplifies it. This is conditional on known interactions, not evidence that specialists intrinsically interact more efficiently.
- consumer, degree_K_rho: 7/10 departures after BH correction across 40 comparisons. Quercus, Nahuel, Gottin_HP, Garraf_HP, Garraf_PP, Olot, Montseny.
- consumer, low_degree_singleton_excess: 8/10 departures after BH correction across 40 comparisons. Quercus, Salix_Galpar, Gottin_HP, Gottin_PP, Garraf_HP, Garraf_PP, Garraf_PP2, Montseny.
- resource, degree_K_rho: 5/10 departures after BH correction across 40 comparisons. Garraf_HP, Garraf_PP, Garraf_PP2, Olot, Montseny.
- resource, low_degree_singleton_excess: 0/10 departures after BH correction across 40 comparisons. .

The within-C permutation assigns observed K among links with identical C. It retains endpoint degrees, pair opportunities, total support, and the entire support distribution. It asks whether low-degree endpoints are associated with low support beyond these opportunities; it cannot explain or predict the distribution itself. Null bands are central 95% allocation-reference intervals, not confidence intervals. A departure does not isolate a causal diet effect. Shared species make links dependent; permutation results are conditional diagnostics under exchangeability within C, not independent-link population tests.

## Interpretation

Restricted occupancy necessarily bounds support: K<=C<=min(endpoint occupancies). A one-site endpoint forces every one of its observed links to be a singleton. But numerous specialist species need not generate most links: each has few partners, so species fractions must be compared with link fractions. The decomposition measures which rare links are forced to be rare by limited opportunities; it is not a causal attribution. K=1 despite C>1 requires something beyond one-site opportunity restriction, potentially unequal encounters, partner compatibility, local context or imperfect detection. It does not identify which. Thresholds K<=2 and K<=3 show whether the result extends beyond singletons. Degree thresholds 1,2,3 are sensitivity descriptions, not universal specialist definitions. Partial rank correlations are descriptive only and do not establish mediation.

## Provenance and limitations

Confirmed restored occurrence records and existing validated interaction records are used. Every K and C was reconstructed and checked; the site universe and prior input hashes were verified. Original input hashes and current source hashes accompany the outputs. Confirmed restorations are retained. Unconfirmed Gottin_PP missing-identity records, Garraf_PP2 zero columns and Montseny zero row remain excluded. Garraf_HP/Olot historical one-link discrepancies remain unresolved. Many species presences originate from interaction sampling, so occupancy–degree associations can include observation coupling. Gottin datasets share sites; Catalan datasets have collection overlap limitations (see copied metadata tables). No independent landscape replication is claimed. Binary records cannot separate ecological absence from nondetection, and spatial occupancy is not habitat specialism. No new gap predictions or ecological mechanisms are fitted.

