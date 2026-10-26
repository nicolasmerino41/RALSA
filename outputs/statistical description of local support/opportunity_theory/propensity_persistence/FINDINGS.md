# Findings

All inputs, known positive links and total support are fixed. No parameter fitting or prediction validation is performed.

Independent propensity draws per opportunity integrate exactly to the common-probability model. Persistent pair propensities are compared at equal pre-conditioning mean (0.5) and variance (ρ/4), using two-state and beta distributions. Positivity and M are then jointly conditioned upon. ρ is the unconditioned within-pair correlation, not the correlation after conditioning.

Central 95% bands describe 300 independent conditional allocations, not confidence intervals. CCDFs are exact conditional expectations; ρ = 0.6 is a prespecified illustration, not a selected best fit.

- Two-state, singleton: increasing ρ from 0 to 0.9 increases the mean in 8/10 datasets; endpoint is closer to empirical in 8/10; empirical lies inside endpoint allocation interval in 4/10.
- Beta, singleton: increasing ρ from 0 to 0.9 increases the mean in 5/10 datasets; endpoint is closer to empirical in 5/10; empirical lies inside endpoint allocation interval in 1/10.
- Two-state, gini: increasing ρ from 0 to 0.9 increases the mean in 8/10 datasets; endpoint is closer to empirical in 8/10; empirical lies inside endpoint allocation interval in 1/10.
- Beta, gini: increasing ρ from 0 to 0.9 increases the mean in 5/10 datasets; endpoint is closer to empirical in 5/10; empirical lies inside endpoint allocation interval in 3/10.
- Two-state, within_C_variance: increasing ρ from 0 to 0.9 increases the mean in 10/10 datasets; endpoint is closer to empirical in 5/10; empirical lies inside endpoint allocation interval in 3/10.
- Beta, within_C_variance: increasing ρ from 0 to 0.9 increases the mean in 10/10 datasets; endpoint is closer to empirical in 6/10; empirical lies inside endpoint allocation interval in 3/10.

The experiment tests whether persistence of differences can generate concentrated support beyond unequal opportunities alone. It does not identify abundance, compatibility or environment as its biological cause, and does not establish a universal distribution family. Read intermediate levels and both families before claiming robustness; an endpoint overshoot is informative, not a fitted-model rejection.

The previous Pair model mixed separately zero-truncated binomials. Here latent probabilities are mixed BEFORE imposing positive support and total M, matching the stated generative experiment. Thus these are new theoretical scenarios, not refits of that model.

Known-link identities and M are supplied; unknown interactions and support amount are not predicted. Only opportunity counts matter for these support statistics because sites are exchangeable within each pair; a full binary matrix need not be generated.

Confirmed occurrence restorations are retained. Uncertain excluded records remain: Gottin PP 474 missing-identity records, Garraf PP2 four zero columns, Montseny one zero row. Interaction-derived occurrence, undocumented sites, previous Gottin HP V1 correction and unresolved one-link discrepancies in Garraf HP/Olot remain limitations. Gottin datasets share sites and other collection relationships remain unresolved. These datasets are not independent landscape replicates. Site count, spatial coverage and sampling effort can affect the distributions.

## Interpretation of the experiment

The shared response is increased conditional variation: at rho=0.9 both families give greater variation at fixed C than the baseline in all ten datasets. Singleton prevalence and overall concentration are less robust: they can decrease or turn down as heterogeneity increases, particularly for the beta family. Garraf HP and Garraf PP show decreasing concentration under both families. Therefore persistent pair differences can generate additional within-opportunity variation, but are not by themselves a general guarantee of the empirical marginal shape.

The conditioning matters. Low-propensity pairs must still realise at least one interaction, and the total support budget must remain fixed; the retained probabilities and their association with C need not resemble the initial latent distribution. This is a mathematical selection effect of the experiment, not evidence for an ecological filtering process. The latent mean of 0.5 is a fixed theoretical control, not a fitted realisation rate; robustness to other latent means has not been tested.

All exhaustive small-case comparisons, unconditioned moment identities, saved Common baseline comparisons, input hashes and conditional support constraints passed. The four PNGs were visually checked for legibility and clipping.