parameter-summary.csv contains posterior means and equal-tailed 95% Bayesian credible intervals.
empirical_bias and interval_covers_reference are blank unless --reference-beta and/or --reference-gamma are supplied.
Starting values are algorithm inputs, not true values, and therefore are not used as bias/coverage references.
For a simulation study, aggregate empirical bias and coverage over independent generated datasets, not only MCMC chains.
Inspect split_rhat (ideally near 1), effective_sample_size, trace plots, and autocorrelation plots before interpreting estimates.
