# Zero-presence sensitivity analysis

This directory is separate from the preserved `outputs/support_prediction` pilot. Its `confirmed/` and `inclusive/` folders each contain inputs, fitted parameters, metrics, seven PNG/PDF diagnostic figures and run logs. `FINDINGS.md` explains the comparison. Original data, scripts and manuscript figures are unchanged.

## Reproduce

From the RALSA root, using Python with pandas and openpyxl installed:

```
python scripts/prepare_zero_presence.py . ../InteractionExtinctionDebt
julia --startup-file=no --project=. scripts/support_prediction_zero.jl
julia --startup-file=no --project=. scripts/support_prediction_zero.jl --inclusive
julia --startup-file=no --project=. scripts/compare_zero_presence.jl
```

`--plot-only` can be combined with either Julia scenario. The parallel modelling script freezes the original pilot's numerical method and weighting; only its input/output paths and occurrence-audit labels change. The original pilot README documents fitting, exact random-removal expectations, weighting, grouped blocks and metrics. No additional ecological models are introduced here.

`presence_record_provenance.csv` records source file, row/sheet, species, site, guild, evidence, uncertainty and whether the presence already existed. Repeated source observations are deduplicated at species-site-guild level. `source_checksums.csv` fingerprints source workbooks/CSVs and Quercus metadata. Names are trimmed as in the existing prepared data, and Montseny resources retain the Plant- prefix. Records with blank or out-of-universe site IDs are not added; this comparison deliberately retains the ten original site universes.

`input_changes.csv` records additions and opportunity changes. Co-occurrences are rebuilt as the within-site Cartesian product of recorded consumers and resources, retaining guild identity. Interactions are copied byte-for-byte. Assertions verify that original co-occurrences and interactions are contained in the rebuilt tables. Source documents establish the confirmed additions; zero padding alone does not establish presence. The inclusive version is an assumption sensitivity, not a complete or validated occupancy census.

`comparison_metrics.csv` stacks all three input versions and all fitting schemes; `comparison_summary.csv` contains equal-dataset averages. Original and new observed T, L and mean gaps match exactly (within numerical tolerance for gaps). All 100 new fits converged without boundary solutions; multistart objectives agreed within 2.38e-14. Independent log-gamma likelihood calculations agreed within 1.89e-15. Local perturbations did not improve fitted likelihoods. Model assertions also check probability normalisation and retention identities.

## Comparison figures

- **01 inputs:** opportunity increases and changes in empirical sum(K)/sum(n), by dataset. Grey is original, blue confirmed, orange inclusive. Connecting segments show the change from baseline.
- **02 transfer accuracy:** equal-dataset average errors under LODO across versions; orange M0, blue M1. NLPD is shown for model comparison within versions: the evaluated pair population changes across versions, so a lower value across versions is not by itself an improvement.
- **03 dataset errors:** signed support and richness errors (%) and gap error (percentage points), for every dataset/version/model under LODO. Blue cells are underestimates, orange overestimates; colour scales differ by metric.

The scenario-specific figures retain the original guide's definitions and colours. Predictions are ratios of expected regional counts, not simulations or confidence intervals. These analyses do not correct interaction detection, sampling effort, or all uncertain occupancy provenance.
