# Effect of restoring zero-interaction presence records

The correction changes model calibration and absolute prediction errors, but does not remove the main modelling limitation: heterogeneous pair propensities improve on a shared probability, yet neither model jointly predicts absolute support, regional richness and the retention gap adequately.

The original pilot already included K=0 pairs among its recorded co-occurrences. This analysis restores species-site presences discarded before co-occurrences were constructed. It does not fill arbitrary zero matrix cells with assumed presences.

## What was included

The confirmed version restores 629 species-site presences: 156 Quercus resources supported by positive unparasitised-host counts; 308 Gottin HP resources with positive host-cell counts; and 115 Garraf HP and 50 Olot resources explicitly retained as non-parasitised hosts according to workbook notes.

The inclusive sensitivity version also adds 474 Gottin PP consumer-site records with named insects but missing plant identities, four named all-zero Garraf PP2 plant columns, and one Montseny plant row. These 479 additional presences have uncertain sampling/presence interpretations and are not promoted to confirmed observations. Missing plant names are never treated as species. Nahuel's regionally aligned matrix padding is not evidence of local presence. Salix's existing independent host records are retained; Garraf PP has no additional zero rows or columns in the prior audit.

All ten datasets are refitted under both versions. Interactions and the original site universes are held fixed to isolate occupancy changes; no P2/P3 parasitoid interaction expansion or unrelated interaction correction is mixed into this comparison. New resource species may enter the co-occurring pool even if they never interacted regionally. All their co-occurring pairs receive K=0.

| Dataset | Confirmed presences added | Opportunities before → after | Empirical sum(K)/sum(n), before → after |
|---|---:|---:|---:|
| Quercus | 156 | 2,886 → 5,161 | 0.149 → 0.083 |
| Gottin HP | 308 | 1,754 → 3,728 | 0.204 → 0.096 |
| Garraf HP | 115 | 529 → 1,079 | 0.301 → 0.147 |
| Olot | 50 | 1,048 → 1,511 | 0.241 → 0.167 |

## Predictive consequences

Equal-dataset mean errors under leave-one-dataset-out prediction:

| Model and metric | Original | Confirmed | Inclusive sensitivity |
|---|---:|---:|---:|
| M0 total-support error | 37.4% | 50.2% | 57.8% |
| M1 total-support error | 33.5% | 40.8% | 48.4% |
| M0 regional-richness error | 100.2% | 119.7% | 130.8% |
| M1 regional-richness error | 28.6% | 38.0% | 44.3% |
| M0 mean-gap error | 6.33 pp | 6.76 pp | 6.89 pp |
| M1 mean-gap error | 3.47 pp | 3.22 pp | 3.21 pp |

Count errors are mean absolute relative errors; gap errors are mean absolute errors in percentage points. M1 improves predictive likelihood over M0 in all ten datasets under both corrected versions, also under grouped holdout. Average NLPD falls across input versions, but this is **not directly interpretable as improved fit** because the evaluation pool gains K=0 pairs. Compare M0 with M1 within the same version; absolute count and gap targets remain comparable across versions.

Effects are uneven. With confirmed additions, M1 overpredicts Quercus total support by 73.4% and richness by 104.9%, while its gap error is only -1.49 pp. Gottin HP has support and richness errors of +51.0% and +75.2%, despite a gap error of just -0.30 pp. Conversely, Garraf HP total-support error improves from -46.9% to -6.3%. Corrections therefore alter the apparent success of individual datasets, not just a common intercept. Even datasets without new presences change predictions because training parameters change.

The inclusive scenario is especially consequential for Gottin PP: its co-occurrence opportunities rise from 19,463 to 29,831 and M1 overpredicts support by 112.1% and richness by 137.1%. Its missing-plant records deserve provenance review before adopting that version as primary.

## Is this only a transfer problem?

No. With confirmed inputs, within-dataset M1 fits have mean absolute richness error of 3.8%, but support error of 25.4% and gap error of 3.76 pp. Inclusive values are 3.7%, 26.0% and 3.88 pp. Within-dataset fits continue to underestimate support and the gap in every dataset. Local likelihood optimisation still does not reproduce the repeated-support structure sufficiently; parameter transfer alone cannot explain the mismatch.

Grouped holdouts give M1 support/richness/gap errors of 49.3%/40.8%/3.44 pp for confirmed inputs and 57.9%/50.1%/3.42 pp for inclusive inputs. The qualitative conclusion is unchanged.

## Implication for the paper

The empirical local-versus-regional interaction-loss gap is unchanged **by construction**: K and N are unchanged. This is not evidence that the gap is immune to sampling or detection. What changes is the inferred opportunity-to-interaction relationship and how well it predicts that gap. Raw empirical realisation rates from discarded-presence inputs were too high in the affected datasets. The corrections do not rescue the shared-parameter models; they make absolute prediction less reliable on average, while normalisation conceals some of that mismatch.

Use the confirmed restoration as the defensible comparison, retain the inclusive version as sensitivity, and resolve missing-partner provenance before treating it as a corrected primary dataset. Neither version establishes complete independent occupancy for both guilds.

Start with `comparison_01_inputs.png`, `comparison_02_transfer_accuracy.png`, and `comparison_03_dataset_errors.png`. Each scenario also contains the seven original diagnostic figure types and complete tables.
