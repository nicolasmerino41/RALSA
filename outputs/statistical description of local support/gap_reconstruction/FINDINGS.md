# Does the existing support fit preserve the removal gap?

**The fitted distribution reconstructs the gap closely in several datasets and captures its broad shape throughout, despite its previously identified distributional discrepancies.** Mean absolute errors range from 0.05 to 1.22 percentage points over the integer removal levels within 0–80%. Residual bias is mainly upward: fitted regional retention, and therefore the gap, is usually too high.

These are reconstructions using parameters already estimated from observed support, not independent predictions. No parameters were refitted.

| Dataset | Mean absolute error (pp) | Maximum absolute error (pp) | Signed mean error (pp) | Removal at maximum error |
|---|---:|---:|---:|---:|
| Quercus | 0.095 | 0.209 | +0.090 | 54.5% |
| Nahuel | 0.049 | 0.105 | +0.049 | 50.0% |
| Salix–Galpar | 0.630 | 1.341 | +0.630 | 73.9% |
| Gottin HP | 1.163 | 2.865 | +1.163 | 78.1% |
| Gottin PP | 0.726 | 1.550 | +0.726 | 71.9% |
| Garraf HP | 0.381 | 0.701 | +0.381 | 64.0% |
| Garraf PP | 1.219 | 2.593 | +1.219 | 75.0% |
| Garraf PP2 | 0.831 | 1.967 | +0.831 | 76.2% |
| Olot | 0.070 | 0.152 | +0.047 | 50.0% |
| Montseny | 0.417 | 0.853 | +0.417 | 72.2% |

## Where the approximation works and where it misses

- **Nahuel, Quercus and Olot:** very close reconstruction throughout the main range; maximum deviations are below 0.21 pp. Quercus and Olot cross into slight underestimation late in that range, ending at -0.040 and -0.136 pp respectively.
- **Garraf HP and Montseny:** maximum errors remain below 0.86 pp, with modest overestimation concentrated at intermediate-to-high removal.
- **Salix–Galpar, Gottin PP and Garraf PP2:** the overall curve is retained, but overestimation reaches 1.34–1.97 pp around 72–76% removal.
- **Gottin HP and Garraf PP:** the largest departures, reaching +2.865 pp at 25/32 sites removed and +2.593 pp at 30/40 sites removed. Their mean absolute errors are approximately 1.2 pp. These fits should not be treated as exact reconstructions, particularly late in the main range.

For Gottin HP and Garraf PP2, the largest main-range error occurs at the final included integer level; the table does not claim those are the largest errors beyond 80%. The full curves CSV includes all removal levels through 100%.

## Why distributional misfit need not imply a large gap error

Regional retention averages survival probabilities over the support distribution. This weighted average is less detailed than the full frequency distribution. At low removal, many repeatedly supported links have survival probabilities close to one, so differences among larger K values contribute little. As removal increases, differences in the repeated-support distribution matter more. The reconstructed singleton fraction is exact, making the expectation after the first site removal exact as well. The fitted mean also matches the observed mean, fixing agreement when only one site remains; both curves reach zero when all sites are removed.

These constraints help explain close reconstruction without requiring the full distribution to be adequate. In particular, Montseny and Salix–Galpar showed detectable distributional misfit previously, yet their main-range mean gap errors are only 0.42 and 0.63 pp. Conversely, Gottin HP and Garraf PP retain noticeable high-removal bias. The direction of that bias reflects the complete redistribution of probability across K, not merely whether the rarest high-support tail is underestimated.

## Scope and caveats

Metrics weight all integer removal levels with m/N<=0.8 equally, including zero. The last included fraction differs by N (75% for Nahuel, 76.2% for Garraf PP2, versus exactly 80% for Garraf HP and Garraf PP). No fractional-removal interpolation is used. Numerical differences are not hypothesis tests or uncertainty estimates.

This establishes how well the existing two-parameter approximation preserves the expected gap for these observed networks. It does not establish transferable prediction or robustness to sampling, detection or different removal geometries. Existing occurrence corrections affect co-occurrence-based quantities but leave these curves unchanged unless interaction records or the site universe change. The current Garraf HP/Olot one-link quality flags and limits of the recorded site universe remain in place.

Every empirical and fitted integer-level regional-retention value was independently checked using exact-integer binomial coefficients; the maximum absolute difference from the stable Julia calculation was 1.67e-15. Summary metrics were independently verified. See `README.md` for definitions and provenance.
