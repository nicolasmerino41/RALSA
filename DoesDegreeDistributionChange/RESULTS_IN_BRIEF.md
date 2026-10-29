# Does support assignment change the degree distribution?

**Yes, but a positive support–degree relationship does not imply universal flattening.**

- Generalists retain a larger expected fraction of their partners in 8/10 consumer datasets at approximately 40% removal and 9/10 at approximately 80%. The association is positive in all ten resource datasets at both levels.
- Nevertheless, degree inequality among connected consumers decreases in all ten datasets at both removal levels. The remaining network becomes less unequal by this measure.
- In eight consumer datasets, inequality remains higher than under shuffled support assignment, at both removal levels. Thus, the observed assignment usually moderates the reduction in inequality. Garraf PP2 and Olot show the opposite direction.

These are counts of directional results, not counts of statistically significant ecological effects. Resource networks are more variable: observed inequality exceeds the shuffled comparison in 7/10 datasets at 40% removal and 8/10 at 80%.

## Read these figures first

1. [Generalists' proportional partner retention](outputs/01_consumer_expected_retention.png).
2. [How degree inequality changes](outputs/03_degree_shape_change.png): negative values mean lower inequality than initially.
3. [Effect of support assignment](outputs/04_support_assignment_effect.png): positive values mean higher inequality than the shuffled comparison.

The full cumulative distributions are in figures 6–9; species disconnection is separate in figure 05. Gini measures inequality, not a fitted tail exponent, so the results do not establish a universal tail shape.

## What the comparison means

The test compares actual site removal with the same removal after shuffling complete support patterns among regional links. This preserves the regional graph, support distribution and interaction richness of every site. It changes which species receive which support patterns, including how their interactions occur together. It does not preserve species-occurrence feasibility.

This is a conditional structural comparison, not independent prediction or proof that the support–degree slope alone causes the result. Full methods, numerical results, input checks and data limitations are in [START_HERE](outputs/START_HERE.md).
