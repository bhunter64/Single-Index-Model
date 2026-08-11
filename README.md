# Single-Index-Model
Package containing all of the tools for new in-development Single Index Model framework

## HER-2, CA19-9, and Axl MCMC analysis

The Compute Canada runner fits four MCMC chains to `bio30` (HER-2), `bio3`
(CA19-9), and `bio37` (Axl), starting gamma at `-0.6, 0.5, 0.5` in that
order. The three beta coefficients default to Cox-derived starting values at
that gamma; set `BETA_START="b1,b2,b3"` to override them.

Submit from the project root:

```sh
sbatch scripts/submit-biomarker-mcmc.sbatch
```

If the data is stored at an absolute `/data` path, specify it at submission:

```sh
DATA_FILE=/data/cleaned_data.rda sbatch scripts/submit-biomarker-mcmc.sbatch
```

The output directory contains the parameter estimates, equal-tailed 95%
credible intervals, split-Rhat and effective sample sizes, starting values,
compressed posterior traces, serialized chains, and downloadable PNG/PDF trace
and autocorrelation plots.

Empirical bias and coverage require known reference values; initial values are
not treated as truths. For simulated or validation data, pass references to the
R runner directly:

```sh
Rscript scripts/run-biomarker-mcmc.R \
  --reference-gamma "-0.6,0.5,0.5" \
  --reference-beta "0.1,0.2,0.3"
```

The same references can be supplied to the batch job as `REFERENCE_GAMMA`
and `REFERENCE_BETA` environment variables.
