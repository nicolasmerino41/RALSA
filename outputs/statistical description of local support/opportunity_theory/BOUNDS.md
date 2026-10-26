# Exact bounds on the expected retention gap

Let f_m(k)=choose(m,k)/choose(N,k) for k≤m and zero otherwise. Then G(m)=m/N−mean(f_m(K)). For independent removal f_r(k)=r^k instead.

## Why the fixed-count result holds

For k≤m, f_m(k+1)/f_m(k)=(m−k)/(N−k). For 1≤k≤N−2, its second difference is

f_m(k+2)−2f_m(k+1)+f_m(k)
= f_m(k)(N−m)(N−m−1)/[(N−k)(N−k−1)] ≥ 0,

with the zero continuation when k>m. For m=N the sequence is constant; for m=N−1 it is linear. Thus the extinction probability is a discrete convex function of support for every fixed-count removal level, including its transition to zero. The independent-removal function r^k is also convex. Endpoints give zero gaps and degenerate envelopes.

This result concerns convex ordering/majorisation, **not variance alone**. Two arbitrary support distributions with equal means need not be ordered by all these convex functions merely because one has greater variance.

## Envelope A

Fix N, L, M and integer supports 1≤K≤N. To maximise G, minimise the sum of convex extinction probabilities: allocate support as equally as possible (floor(M/L) and ceiling(M/L)). To minimise G, maximise that sum: start every link at one and fill links to N sequentially, with at most one intermediate link.

## Envelope B

Also fix each link's individual capacity C_l. For the upper gap bound, use capped equalisation: start all links at one, increment the smallest nonsaturated supports, and stop when their sum is M. Links whose capacity is below the common water level stay capped; all uncapped supports differ by at most one. A feasible transfer from a support at least two larger to a smaller nonsaturated support cannot increase a convex sum. Repeating such transfers reaches this balanced allocation.

For the lower gap bound, sort capacities decreasingly, start all links at one, and fill those capacities in that order. Its descending partial sums are as large as feasibility permits: for the largest t supports their sum cannot exceed either the sum of the t largest capacities or M−(L−t). The greedy vector achieves this bound and therefore majorises every feasible vector. It maximises every convex sum of support values.

The lower and upper configurations depend on N, M and capacities, **not on m**. Consequently one vector attains each pointwise bound across the entire curve. They also optimise any nonnegative weighted sum of gap values, including the reported discrete main-range integrated gap. Pointwise optimisations do not switch incompatible vectors as removal increases.

The opportunity matrix introduces no further restriction on these support vectors: each link independently chooses K_l sites among its C_l eligible sites. Links may share sites. Fixing site totals as well would be a different allocation problem and is not done here.

Envelope A needs empirical total activity M; envelope B additionally needs observed opportunities and link identities. They are nested conditional feasible ranges, not independently predicted intervals. An empirical curve's membership is guaranteed once the constraints are satisfied.

## Positions and checks

Z(m)=(G_emp−G_min)/(G_max−G_min). Z is undefined (NaN) when width≤10⁻¹⁰. The integrated position is the ratio formed from integrated empirical and extremal gaps, **not the average of pointwise Z**. The integration measure is the uniform distribution on integer removal counts m with m/N≤0.8, including zero.

The Julia implementation enumerates every labelled support configuration for every capacity vector with N=2,…,5 and L=1,…,4: 67,024 feasible capacity/support combinations. It checks every fixed-count removal level and independent probabilities 0.1, 0.37 and 0.8. Every case passes. Independent Python calculations use exact integer binomial coefficients to check the empirical and extremal exported curves; the largest difference is below 6×10⁻¹⁶.

These are exact mathematical results for the stated removal process and constraints. They do not validate any assumed realisation mechanism, geographical removal model, or ecological interpretation.
