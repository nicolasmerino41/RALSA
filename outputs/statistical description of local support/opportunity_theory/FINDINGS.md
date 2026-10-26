# What opportunities and realisation rules explain

**Opportunities place useful, sometimes strong limits on the gap. Persistent pair differences improve prediction within systems, but none of these small generative rules adequately transfers the complete support pattern across all ten datasets.** The exact bounds are established; the empirical model comparison remains exploratory.

## Exact constraints and empirical positions

The gap definition agrees with the existing analysis. Both envelopes are exact: concentrating support minimises the gap; capped equalisation maximises it. One configuration achieves each extreme throughout the curve and for the integrated gap. Exhaustive checks passed 67,024 small cases; independent integer-binomial calculations agreed to within 6×10⁻¹⁶.

Adding observed opportunity capacities to fixed N, L and M reduces mean envelope width by **8.5–79.1%**. The reduction is smallest in Gottin HP (8.5%) and Gottin PP (11.7%), intermediate in most systems (23–36%), and largest in Garraf PP2 (79.1%). Opportunities therefore constrain the range very differently among systems; they do not generally determine the gap.

Empirical integrated positions in the opportunity-constrained envelope lie between **0.435 and 0.621**, averaging **0.511**. They occupy a relatively similar middle portion of their feasible ranges. This is a descriptive regularity, not evidence for a universal mechanism. Membership in the envelope is guaranteed, and both envelopes use empirical M.

## What the controlled rules distinguish

At equal **expected** activity, persistent pair differences substantially increase singleton frequency and reduce the gap in most systems. At h=3, the reduction relative to the common rule ranges from 0.03 pp in Garraf PP2 to 6.44 pp in Gottin HP. Shared site differences change mean gaps much less: 0.02–0.69 pp. In these specified rules, repeated differences among pairs redistribute support more strongly than independent differences among sites.

The near-absence of this response in Garraf PP2 is informative: opportunities and activity already leave little room to redistribute support. This connects a conditional mathematical restriction to a scenario's response; it is not validation of that scenario. Expected M is matched analytically. Across 96 replicates, average realised activity differs from its target by at most 2.45%; each replicate's variation is retained in the tables.

Recovery uses the actual Nahuel, Gottin HP and Salix–Galpar opportunities. Strong pair/site scenarios (h=2) are correctly classified in all 24 replicates per matrix. Weak pair heterogeneity (h=0.5) is recovered in 71–92%, and weak site heterogeneity in 83–100%. However, the true common rule is recognised in only 25% of Nahuel replicates, versus 79% for Gottin HP and 96% for Salix–Galpar. Small samples and unpenalised model comparison readily mistake no heterogeneity for weak heterogeneity. Parameter recovery is also imprecise in Nahuel. These on-grid recovery experiments establish limited distinguishability of the chosen scenarios, not identification of arbitrary ecological mechanisms.

## Predictions with withheld observations

Equal-dataset averages:

| Rule | Entire-dataset gap MAE (pp) | Within-system gap MAE (pp) | Within-system cell Brier score |
|---|---:|---:|---:|
| Common probability | 4.56 | 1.21 | 0.202 |
| Persistent pair differences | 2.88 | 0.59 | 0.180 |
| Shared site differences | 3.42 | 1.19 | 0.205 |
| Previous common-parameter transfer | 3.29 | — | — |

Within-system pair predictions improve gap MAE in **9/10 datasets**; Olot worsens. Cell Brier scores improve in 8/10. This is useful evidence that observations of a pair elsewhere in the same system help predict its withheld occurrences. Approximately 80% of sites are already observed, so the resulting full-landscape gap completion is an easier task than transfer to an unseen system.

For entirely withheld datasets, the pair rule improves on the previous transfer baseline in **6/10**. Its average gain is 0.40 pp (12.3%), with mean relative error 21.5%. The common and site rules average 35.0% and 24.6% relative error. The pair rule's held-out cell Brier score (0.241) is almost identical to the common rule (0.241): better gap accuracy does not imply a comparably large gain in predicting individual occurrences.

| Dataset | Previous transfer MAE | Pair-rule MAE | Change (pp) |
|---|---:|---:|---:|
| Quercus | 2.35 | 0.98 | −1.37 |
| Nahuel | 1.06 | 1.75 | +0.69 |
| Salix–Galpar | 3.00 | 1.05 | −1.95 |
| Gottin HP | 1.49 | 6.36 | +4.87 |
| Gottin PP | 1.84 | 5.08 | +3.23 |
| Garraf HP | 6.46 | 1.39 | −5.07 |
| Garraf PP | 4.20 | 0.63 | −3.57 |
| Garraf PP2 | 6.07 | 9.02 | +2.95 |
| Olot | 1.61 | 0.43 | −1.19 |
| Montseny | 4.76 | 2.13 | −2.63 |

## Specific failures and interpretation

Gottin HP and PP are predicted to have too few singletons, too many repeated links and too much activity. The pair rule predicts roughly 587 versus 357 observed occurrences in Gottin HP, and 3,181 versus 1,987 in Gottin PP. Garraf PP2 fails in the opposite direction: about 2,275 versus 4,288 occurrences, with 66.7% singletons predicted versus 39.4% observed. Its broad repeated-support tail is missed. These are failures of transferring a shared propensity distribution, not simply unexplained differences in gap size.

Even a close gap curve can coexist with wrong activity or support tails. The pair rule's conditional prediction intervals contain observed M in only 5/10 systems and singleton frequency in 4/10. The site rule usually predicts excessive clustering after support is controlled, despite matching that conspicuous feature in Garraf PP2. Thus broad site-wide propensity differences are not a general explanation of these records. Fitting a combined model is deferred: the present results do not yet support a clean general attribution to both components.

The new contribution is a separation of **what opportunities permit**, **how persistent versus site-wide realisation differences behave at controlled activity**, and **which information transfers**. Pair-level information has within-system predictive value, but a shared pair distribution is not a sufficient cross-system theory. The findings neither identify the ecological cause of those differences nor reject every opportunity-based explanation.

## Limits and discriminating observations

All models condition explicitly on known regional link identities and positivity; none predicts wholly undetected interactions. Occurrence provenance, shared Gottin sites and unresolved provider-collection overlap limit independence. The finite grid, approximate synthetic likelihood, one site split per dataset, and prediction intervals conditional on point estimates make this a preliminary test. Three within-system estimates hit grid boundaries; no leave-one-dataset-out estimate does. Monte Carlo error in mean predicted gaps is at most about 0.09 pp, separate from parameter and sampling uncertainty.

The most useful additional observations would be independent species-presence surveys with repeated interaction sampling; repeat observations of the same pairs across environmental contrasts; and measured sampling effort, phenology or abundance at each site. These would help distinguish persistent pair differences, shared site conditions and detection effects. More accurate gap reconstruction alone would not distinguish them.
