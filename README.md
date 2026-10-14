# Southern California Bight (SCB) Sediment Core Correlation State-Space Model (SSM)

This repository contains an R and Stan workflow for comparing fish scale deposition rates (SDRs) from four sediment cores and assessing their correlation through time. These cores were taken from Santa Monica and Santa Barbara Basins.

## Project layout
All files live in a single project folder, `SCB_Sed_Core_SSM/`:
- `SCB_Sed_Core_SSM.Rproj`: RStudio project file; marks the project root for `here()`.
- `SCB_Sed_Core_SSM_fit.R`: data preparation and model fitting.
- `SCB_Sed_Core_SSM_plot.R`: posterior summaries and figures.
- `scale_ssm_model.stan`: the Stan model (and its compiled executable, created on first run).
- Four Excel data files, one per core (see the top rows of each file for metadata):
  - `Santa_Monica_Basin.xlsx`
  - `SH2015_Santa_Barbara_Basin.xlsx`
  - `Baumgartner1992_Santa_Barbara_Basin_214.xlsx`
  - `Baumgartner1992_Santa_Barbara_Basin_224.xlsx`
- `stan_run/`: outputs produced by `SCB_Sed_Core_SSM_fit.R` (created automatically).
- `stan_run/figs/`: the three output figures produced by `SCB_Sed_Core_SSM_plot.R` (created automatically).

All file paths are built with the [`here`](https://here.r-lib.org/) package, which finds the project root from the `.Rproj` file. You do not need to set a working directory or edit any paths.

## Dependencies
- R
- [CmdStan](https://mc-stan.org/cmdstanr/)
- Packages: here, readxl, dplyr, tidyr, tibble, posterior, ggplot2, patchwork, scales, cmdstanr

```r
install.packages(c("here", "readxl", "dplyr", "tidyr", "tibble", "posterior", "ggplot2", "patchwork", "scales"))

install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan(cores = 4)
```
## Recommended workflow
Run the scripts in this order:
1. Verify that the Excel files and `scale_ssm_model.stan` are in the project folder. Inspect the header rows as needed.
2. Open `SCB_Sed_Core_SSM.Rproj` in RStudio. This sets the project root for `here()`. If not using RStudio, run `here::set_here()` once in the project folder.
3. Run `SCB_Sed_Core_SSM_fit.R`. This loads the data, builds the Stan data list, compiles and samples the model, and saves fit objects and diagnostics to `stan_run/`.
4. Run `SCB_Sed_Core_SSM_plot.R`. This loads the saved fit and produces three figures in `stan_run/figs/`.

## File details
### `scale_ssm_model.stan`
The Stan source code defining the model and its associated likelihoods, detection limits, and priors.

### `SCB_Sed_Core_SSM_fit.R`
The data pre-processing and execution script that contains code to:
- load the four cores' SDR records from Excel
- build a common time grid; handle left-censored zeros
- construct the Stan data list
- configure the sampler parameters
- compile and sample `scale_ssm_model.stan` via cmdstanr (falls back to rstan if cmdstanr is unavailable)
- print sampler diagnostics and a summary of the main parameters

### `SCB_Sed_Core_SSM_plot.R`
Post-processing script that reads the posterior draws, calculates posterior summaries, and generates three figures in `stan_run/figs/`:
- **`fig_state_trajectories.png`**: Latent log-SDR trajectories per core with 50% and 90% credible-interval ribbons. Observed data points are overlaid, as are the censored below-detection-limit observations. Periods before and after each core's observed range are grayed out.
- **`fig_typical_and_pairs.png`**: Posterior histogram of the typical cross-core correlation (`typical_r`), alongside a 4 × 4 pairwise correlation matrix showing posterior medians and 90% credible intervals for each core pair.
- **`fig_scalar_forest.png`**: Per-core forest plot for the parameters μ, ρ, σ_proc, and τ_obs.

## Note
A repository with all of the photographed scales from SMB, as well as files containing their surface area measurements and SDR calculation, can be found [here](10.6084/m9.figshare.32988020).
