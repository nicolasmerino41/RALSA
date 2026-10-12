# Does the common parameter baseline transfer?

**Transfer retains the broad hump-shaped gap but loses much of the accuracy of the own-fit reconstruction.** With each dataset weighted equally, mean absolute gap error rises from **0.56 pp in-sample to 3.29 pp under transfer**, a deterioration of 2.73 pp. MAE increases in all ten datasets. Median dataset MAE rises from 0.52 to 2.68 pp. No alternative averaging method was selected or tested.

Each target receives the arithmetic means of s and r from the other nine datasets, with its own N used only to bound and normalise the support distribution. All expectations are evaluated analytically at integer removal levels. This is a held-out common-parameter baseline conditional on N, not prediction from co-occurrence or regional network structure.

## Dataset variation

All errors are percentage points over integer m/N<=0.8. The signed maximum below is the error at the level with greatest absolute error.

| Dataset | Transfer MAE | Own-fit MAE | Transfer signed mean | Transfer signed maximum | Removal at maximum |
|---|---:|---:|---:|---:|---:|
| Quercus | 2.35 | 0.09 | -2.35 | -4.81 | 77.3% |
| Nahuel | 1.06 | 0.05 | +1.06 | +2.87 | 75.0% |
| Salix–Galpar | 3.00 | 0.63 | -3.00 | -4.53 | 72.4% |
| Gottin HP | 1.49 | 1.16 | -1.49 | -2.16 | 78.1% |
| Gottin PP | 1.84 | 0.73 | +1.84 | +3.19 | 68.8% |
| Garraf HP | 6.46 | 0.38 | +6.46 | +10.53 | 72.0% |
| Garraf PP | 4.20 | 1.22 | +4.20 | +7.05 | 72.5% |
| Garraf PP2 | 6.07 | 0.83 | -6.07 | -11.58 | 76.2% |
| Olot | 1.61 | 0.07 | -1.61 | -2.63 | 78.6% |
| Montseny | 4.76 | 0.42 | +4.76 | +7.85 | 72.2% |

Nahuel has the smallest transfer MAE, although its maximum deviation approaches 2.9 pp. Gottin HP, Olot and Gottin PP have MAE below 2 pp. The weakest transfers are Garraf HP and Garraf PP2, with mean errors above 6 pp and maximum errors above 10 pp. Errors generally grow with removal, reaching their main-range maxima around 69–79%. Some maxima occur at the last included integer level; these are not assertions about the full 0–100% range.

The equal-dataset mean signed error is only **+0.38 pp**, but this is cancellation: five datasets are overestimated and five underestimated. Mean absolute dataset bias is **3.29 pp**. The average dataset maximum error is 5.72 pp, compared with 1.23 pp in-sample; the largest transfer error is 11.58 pp. Gottin HP is an exception on the maximum-error metric alone: its maximum decreases from 2.86 to 2.16 pp, while its MAE still increases and bias changes sign.

## What transfer does and does not preserve

Transferred s ranges narrowly from 0.531 to 0.566 and r from 0.715 to 0.741, whereas own-fit estimates span s=0.394–0.708 and r=0.611–0.842. The shared baseline therefore averages over meaningful differences in observed support. For example, Garraf HP receives fewer singletons and a slower tail decline than its own fit, overestimating the gap; Garraf PP2 receives more singletons and a faster decline, underestimating it. These are statistical consequences of the specified model, not ecological explanations for the parameter differences.

The earlier close reconstruction showed that two fitted parameters can preserve this particular expected curve. It did not imply those parameters were transferable. The present result distinguishes those questions: the common baseline gives a rough description, but dataset-specific support information remains important for accurate gap values. It does not show that every possible transfer approach would fail.

## Validation limitations

This is **not ten independent landscape validations**. Gottin HP and PP share all 32 documented sites; one remains in training when the other is tested. The Garraf HP/PP/PP2, Olot and Montseny datasets share a documented provider collection; exact campaign/site overlap is not fully resolved, particularly among the Garraf datasets. Their folds retain related collection members in training. The metadata audit and source hashes are saved alongside the results. No grouped validation was substituted for the requested procedure.

Primary metrics use equally weighted integer removal levels within each dataset, ending at the last m/N<=0.8. Different site counts give different grids and endpoints. The analysis does not standardise sampling extent, effort or detection, infer a mechanism, or estimate wholly undetected interactions.

All existing data-quality flags remain: unresolved one-link discrepancies for Garraf HP and Olot, removal of the fictitious Gottin HP V1 consumer, and outstanding primary-occupancy revisions. Occurrence changes alone do not change these calculations when K and N are fixed.

Verification confirmed the exact nine-dataset arithmetic means and equal weights in every fold, no target parameter in training, agreement with every previous own-fit curve, and agreement of transferred curves with an independent exact-binomial calculation within 4.45e-16. Summary metrics were independently checked.
