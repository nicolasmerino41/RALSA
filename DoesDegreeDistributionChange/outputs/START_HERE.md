# Does the degree distribution change?

Read figures 01–02 for proportional retention, 03 for shape changes, and 04 for the effect of assigning support to particular species. Figure 05 separates species disconnection. Figures 6–9 compare mean-scaled cumulative distributions.

## What was tested

For each dataset, consumers and resources are analysed separately. Exact expected fractions of partners retained use every link's K. Then 2000 paired replicates compare the observed site–link matrix with randomly reassigned complete columns. Both scenarios use the same removed sites, and the same assignment is used at both removal levels within each replicate. Removal levels are the nearest integer to 40% and 80%; exact fractions are reported.

The shuffle preserves the regional graph, the entire K distribution, total support, every site's interaction count, and the number of surviving regional links for each removal sequence. It changes which species receive those support patterns, including their grouping within species. It is a structural null, not an occupancy-feasible alternative landscape. Thus differences cannot be attributed exclusively to the slope of mean K versus degree.

## Results

|Role|Target removal|Positive degree–retention correlation|Observed mean Gini increases / decreases|Observed mean Gini higher / lower than shuffle|
|---|---:|---:|---:|---:|
|consumer|40%|8/10|0 / 10|8 / 2|
|consumer|80%|9/10|0 / 10|8 / 2|
|resource|40%|10/10|2 / 8|7 / 3|
|resource|80%|10/10|1 / 9|8 / 2|

Counts describe directions, not independent significant replications. Gini is the prespecified shape statistic: higher values indicate greater inequality among connected species, not necessarily a power-law exponent change or a heavier upper tail specifically. Mean-scaled curves show where distributions differ.

|Dataset|Role|Actual removal|Degree–retention rho|Observed ΔGini|Shuffled ΔGini|Observed − shuffled|
|---|---|---:|---:|---:|---:|---:|
|Quercus|consumer|40.9%|0.066|-0.030|-0.037|0.007|
|Quercus|consumer|81.8%|0.268|-0.106|-0.131|0.024|
|Quercus|resource|40.9%|0.553|-0.015|-0.015|0.000|
|Quercus|resource|81.8%|0.594|-0.064|-0.071|0.008|
|Nahuel|consumer|37.5%|0.390|-0.017|-0.038|0.020|
|Nahuel|consumer|75.0%|0.498|-0.066|-0.119|0.053|
|Nahuel|resource|37.5%|0.269|-0.008|-0.004|-0.004|
|Nahuel|resource|75.0%|0.162|-0.014|-0.026|0.012|
|Salix_Galpar|consumer|39.9%|0.681|-0.010|-0.024|0.014|
|Salix_Galpar|consumer|80.1%|0.713|-0.051|-0.096|0.044|
|Salix_Galpar|resource|39.9%|0.439|0.014|0.002|0.011|
|Salix_Galpar|resource|80.1%|0.490|0.034|0.009|0.025|
|Gottin_HP|consumer|40.6%|-0.091|-0.033|-0.039|0.006|
|Gottin_HP|consumer|81.2%|0.014|-0.118|-0.141|0.023|
|Gottin_HP|resource|40.6%|0.076|-0.008|-0.022|0.014|
|Gottin_HP|resource|81.2%|0.237|-0.060|-0.094|0.034|
|Gottin_PP|consumer|40.6%|0.426|-0.019|-0.032|0.013|
|Gottin_PP|consumer|81.2%|0.526|-0.080|-0.123|0.043|
|Gottin_PP|resource|40.6%|0.676|-0.016|-0.027|0.012|
|Gottin_PP|resource|81.2%|0.724|-0.066|-0.104|0.038|
|Garraf_HP|consumer|40.0%|0.288|-0.043|-0.053|0.010|
|Garraf_HP|consumer|80.0%|0.408|-0.112|-0.173|0.060|
|Garraf_HP|resource|40.0%|0.800|-0.022|-0.046|0.024|
|Garraf_HP|resource|80.0%|0.799|-0.127|-0.167|0.040|
|Garraf_PP|consumer|40.0%|0.553|-0.025|-0.050|0.025|
|Garraf_PP|consumer|80.0%|0.592|-0.088|-0.160|0.072|
|Garraf_PP|resource|40.0%|0.839|-0.042|-0.028|-0.013|
|Garraf_PP|resource|80.0%|0.855|-0.136|-0.110|-0.026|
|Garraf_PP2|consumer|38.1%|0.146|-0.020|-0.015|-0.005|
|Garraf_PP2|consumer|81.0%|0.339|-0.084|-0.061|-0.023|
|Garraf_PP2|resource|38.1%|0.171|0.005|-0.005|0.010|
|Garraf_PP2|resource|81.0%|0.263|-0.008|-0.022|0.014|
|Olot|consumer|42.9%|-0.079|-0.027|-0.015|-0.011|
|Olot|consumer|78.6%|-0.045|-0.103|-0.074|-0.030|
|Olot|resource|42.9%|0.724|-0.010|-0.033|0.022|
|Olot|resource|78.6%|0.798|-0.030|-0.099|0.069|
|Montseny|consumer|38.9%|0.585|-0.030|-0.036|0.007|
|Montseny|consumer|77.8%|0.669|-0.087|-0.128|0.041|
|Montseny|resource|38.9%|0.766|-0.024|-0.015|-0.009|
|Montseny|resource|77.8%|0.932|-0.096|-0.070|-0.026|

## Reading the comparison

- If initially high-degree species retain a greater fraction of links, that is more than merely retaining more links in absolute numbers.
- If observed and shuffled shapes differ, assignment of support to the graph matters beyond the marginal support distribution.
- Similar shapes can coexist with unequal species-level retention; distributional stability is not evidence of equal erosion.
- Higher inequality is not synonymous with flattening every part of a cumulative curve. No distribution-family fits or exponent tests were selected after viewing results.

## Definitions and uncertainty

Degree zero is excluded from the connected-species shape calculation but retained in disconnection statistics and full raw histograms. Each replicate's degrees are divided by its own connected-species mean before averaging CCDFs. This is not the distribution of expected degrees. Original connected species define the denominator of disconnection. Low degree means initial degree <= the within-dataset/guild median, with ties kept together; high degree is above that median.

Figure 03 intervals are central 95% removal/allocation ranges. Figure 04 grey intervals are central 95% paired outcome ranges; blue intervals are mean +/-1.96 Monte Carlo standard errors. Those show simulation precision only, not ecological confidence or independent-dataset inference. No permutation significance claim is made from one shuffle per removal sequence.

Empty-remnant replicate count across scenarios: 0. Gini/CCDF condition on nonempty remnants; disconnection includes all replicates. CCDF display ends at five times the mean; unscaled integer-degree histograms include the full tail.

## Input version and scope

All K values were reconstructed from current validated interactions and checked against the retained confirmed/pairs.csv; prior input hashes were verified. The site universe is reconstructed from confirmed occupancy and checked against the saved N for every dataset; the removed intermediate site-universe table is not required. Outstanding audit issues are copied, including Salix's 33 additional sampled zero-interaction localities and unincorporated Gottin_HP parasite slots. These results do not silently correct those inputs. The two Göttingen datasets share sites; other collection dependencies remain. Sampling and interaction-derived occupancy limitations still apply. No inference of intrinsic diet specialism or independent prediction of a new system is made.

This is conditional prediction from observed spatial support. The exact formula determines species' expected degrees; the full site patterns determine joint loss and realised distributions. A K distribution plus a correlation coefficient alone does not uniquely determine these outcomes.

