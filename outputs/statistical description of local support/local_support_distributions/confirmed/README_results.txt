LOCAL-SUPPORT DISTRIBUTION ANALYSIS

Input: outputs/support_prediction_zero/confirmed/, validated against site-level tables; stored SHA-256 file hashes. Source SHA checks are enforced when input_versions.csv is present.
Primary inference: all co-occurring pairs, K>=0. Fit binomial and beta-binomial K|n using exact observed opportunity counts.
Beta-binomial p is the across-pair mean realisation probability; phi is concentration (alpha=p*phi, beta=(1-p)*phi).
LOO: parameters trained on nine whole datasets and applied without recalibration to the tenth; target n is supplied, target K is never used to fit.
Likelihood uses independent pair-level contributions. Shared sites/species can create dependence: predictive scores are useful, but naive standard errors and AIC certainty would be optimistic.
Empirical PMFs include zeros or condition on K>0. Descriptive families (finite geometric, power law, shifted negative binomial) condition on K>0 and truncate at 1:N.
Power-law fits are descriptive candidates only, NOT proof of power-law behaviour or a formal Clauset-style tail test.
Descriptive parameter_1: geometric lambda, power-law alpha, or negative-binomial log(shape). Descriptive parameter_2: negative-binomial logit(success probability).
LOO parametric envelopes hold n, model parameters and candidate pairs fixed. They exclude model-parameter uncertainty and site/species dependence.
Site subsampling uses 40 replicates at min dataset site count=8 (unless disabled). It does not standardise spatial grain or study extent.
Zero support means unobserved interaction among an inferred co-occurring pair, not necessarily a confirmed biological absence.
Binomial and beta-binomial models abstract from network degree structure, species identities, spatial correlations, and detection differences.
Simulations use seed=20261014 and 100 replicates per dataset per LOO model.
All CSV and PNG outputs are generated inside the named scenario output folder.
