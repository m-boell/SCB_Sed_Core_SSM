# SCB_Sed_Core_SSM_plot.R

# Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(posterior)
  library(scales)
  library(here)
})

# Path settings
OUT_DIR <- here("stan_run")
FIG_DIR <- file.path(OUT_DIR, "figs")
dir.create(FIG_DIR, showWarnings = FALSE)

# Load fit
core_levels  <- readRDS(file.path(OUT_DIR, "core_levels.rds")) 
grid_centers <- readRDS(file.path(OUT_DIR, "grid_centers.rds")) 
dat_used     <- readRDS(file.path(OUT_DIR, "dat_used.rds"))
pair_grid    <- readRDS(file.path(OUT_DIR, "pair_grid.rds")) 
stan_data    <- readRDS(file.path(OUT_DIR, "stan_data.rds"))

C <- length(core_levels) 
T_grid <- length(grid_centers) 
t_c_lo <- stan_data$t_lo 

t_c_hi <- stan_data$t_hi 
L_per_core <- stan_data$L_det 

cmdstan_path <- file.path(OUT_DIR, "fit_cmdstan.rds")
rstan_path   <- file.path(OUT_DIR, "fit_rstan.rds")

if (file.exists(cmdstan_path)) {
  fit <- readRDS(cmdstan_path)
  draws_array <- fit$draws(format = "draws_array")
} else if (file.exists(rstan_path)) {
  fit <- readRDS(rstan_path)
  draws_array <- posterior::as_draws_array(rstan::extract(fit, permuted = FALSE))
} else {
  stop("No fit object found in stan_run/. Run SCB_Sed_Core_SSM_fit.R first.")
}

post_array <- function(name) {
  posterior::subset_draws(draws_array, variable = name) |>
    posterior::merge_chains() |>
    posterior::as_draws_matrix()
}

# Plot theme
theme_report <- theme_classic(base_size = 12) + theme(
  panel.grid.minor = element_blank(),
  strip.background = element_rect(fill = "grey92", color = NA),
  strip.text = element_text(face = "bold", size = 11.3),
  plot.title = element_text(face = "bold", size = 13.3)
)

# SDR trajectories plot
x_post <- post_array("x") 
n_draws <- nrow(x_post) 

x_arr   <- array(NA_real_, dim = c(n_draws, T_grid, C))
for (c in seq_len(C)) {
  for (t in seq_len(T_grid)) {
    col_name <- sprintf("x[%d,%d]", t, c)
    if (col_name %in% colnames(x_post)) {
      x_arr[, t, c] <- x_post[, col_name]
    }
  }
}

# Posterior summaries per (t, c)
x_summary <- expand.grid(t = seq_len(T_grid), c = seq_len(C)) |>
  mutate(year = grid_centers[t],
         core = core_levels[c],
         med  = mapply(function(t, c) median(x_arr[, t, c]), t, c),
         q05  = mapply(function(t, c) quantile(x_arr[, t, c], 0.05), t, c),
         q25  = mapply(function(t, c) quantile(x_arr[, t, c], 0.25), t, c),
         q75  = mapply(function(t, c) quantile(x_arr[, t, c], 0.75), t, c),
         q95  = mapply(function(t, c) quantile(x_arr[, t, c], 0.95), t, c)) |>
  mutate(core = factor(core, levels = core_levels))

# Cores' temporal range
range_df <- tibble(core = core_levels,
                   year_lo = grid_centers[t_c_lo],
                   year_hi = grid_centers[t_c_hi]) |>
  mutate(core = factor(core, levels = core_levels))
x_summary <- x_summary |>
  left_join(range_df, by = "core") |>
  filter(year >= year_lo & year <= year_hi)

# Observed points
obs_df <- dat_used |>
  mutate(core = factor(core_levels[core_idx], levels = core_levels),
         year = bin_median_year_AD,
         kind = ifelse(obs_indic == 1, "observed", "censored"),
         y_plot = ifelse(obs_indic == 1, y_log, L_per_core[core_idx]))

# Gray-outed areas beyond a core's temporal range
mask_df <- tibble(core = core_levels,
                  left_xmin  = grid_centers[1] - 10,
                  left_xmax  = grid_centers[t_c_lo],
                  right_xmin = grid_centers[t_c_hi],
                  right_xmax = grid_centers[T_grid] + 10) |>
  mutate(core = factor(core, levels = core_levels))

p_states <- ggplot() +
  geom_rect(data = mask_df, aes(xmin = left_xmin, xmax = left_xmax,ymin = -Inf, ymax = Inf), fill = "grey88", inherit.aes = FALSE) +
  geom_rect(data = mask_df, aes(xmin = right_xmin, xmax = right_xmax,ymin = -Inf, ymax = Inf), fill = "grey88", inherit.aes = FALSE) +
  geom_ribbon(data = x_summary, aes(x = year, ymin = q05, ymax = q95),fill = "#2c7fb8", alpha = 0.22) +
  geom_ribbon(data = x_summary, aes(x = year, ymin = q25, ymax = q75), fill = "#2c7fb8", alpha = 0.35) +
  geom_line(data   = x_summary, aes(x = year, y = med), color = "#08519c", linewidth = 0.55) + 
  geom_hline(data = tibble(core = factor(core_levels, levels = core_levels), L = L_per_core), aes(yintercept = L), color = "red", linetype = "dotted", linewidth = 0.3, alpha = 0.6) +
  geom_point(data = obs_df |> filter(kind == "observed"),aes(x = year, y = y_plot), color = "black", size = 0.7) +
  geom_point(data = obs_df |> filter(kind == "censored"), aes(x = year, y = y_plot), color = "red", shape = 25, fill = "red", size = 1.1) + 
  facet_wrap(~ core, ncol = 1, scales = "free_y", strip.position = "right") + 
  scale_x_continuous(breaks = scales::breaks_width(200)) +
  labs(x = "Year AD", y = "ln(SDR)", title = "Latent ln(SDR) posterior (median, 50% and 90% CrI)") +
  theme_report

ggsave(file.path(FIG_DIR, "fig_state_trajectories.png"), p_states,
       width = 9.5, height = 9, dpi = 600)

# Typical correlation plot
z_mean <- as.numeric(post_array("z_mean")) 
typ_r  <- tanh(z_mean) 

p_typr <- ggplot(tibble(typ_r), aes(typ_r)) +
  geom_histogram(bins = 45, fill = "#2c7fb8", alpha = 0.7, color = NA) +
  geom_vline(xintercept = median(typ_r), color = "#c0392b", linewidth = 0.7) +
  geom_vline(xintercept = quantile(typ_r, c(0.05, 0.95)),
             color = "#c0392b", linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
  coord_cartesian(xlim = c(-0.6, 0.9)) +
  labs(x = expression(tanh(italic(z[0]))~"  = typical cross-core correlation"),
       y = "posterior density",
       title = sprintf("(A) Typical correlation:\nmedian %.3f, 90%% CrI [%.3f, %.3f]",
                       median(typ_r), quantile(typ_r, 0.05), quantile(typ_r, 0.95))) +
  theme_report

# Correlation heat map plot
R_off <- post_array("R_offdiag")   
R_off_summary <- pair_grid |>
  mutate(k = row_number(),
         med = apply(R_off, 2, median),
         q05 = apply(R_off, 2, quantile, 0.05),
         q95 = apply(R_off, 2, quantile, 0.95))

R_full <- diag(C); rownames(R_full) <- core_levels; colnames(R_full) <- core_levels
R_full_lo <- R_full_hi <- R_full
for (kk in seq_len(nrow(R_off_summary))) {
  i <- R_off_summary$i[kk]; j <- R_off_summary$j[kk]
  R_full[i, j]    <- R_full[j, i]    <- R_off_summary$med[kk]
  R_full_lo[i, j] <- R_full_lo[j, i] <- R_off_summary$q05[kk]
  R_full_hi[i, j] <- R_full_hi[j, i] <- R_off_summary$q95[kk]
}

heat_df <- as.data.frame(as.table(R_full)) |>
  rename(core_i = Var1, core_j = Var2, R_med = Freq) |>
  mutate(R_lo = as.numeric(as.table(R_full_lo)),
         R_hi = as.numeric(as.table(R_full_hi)),
         is_diag = as.character(core_i) == as.character(core_j),
         label_val = ifelse(is_diag, "1", sprintf("%.2f", R_med)),
         label_ci  = ifelse(is_diag, "",  sprintf("[%.2f, %.2f]", R_lo, R_hi)),
         txt_color = ifelse(abs(R_med) > 0.55, "white", "black"),
         i_idx = match(as.character(core_i), core_levels),
         j_idx = match(as.character(core_j), core_levels)) |>
  filter(i_idx <= j_idx) |>
  mutate(core_i = factor(core_i, levels = core_levels),
         core_j = factor(core_j, levels = rev(core_levels)))

p_heat <- ggplot(heat_df, aes(x = core_i, y = core_j, fill = R_med)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label_val, color = txt_color),
            size = 4.4, fontface = "bold", vjust = -0.3) +
  geom_text(aes(label = label_ci, color = txt_color),
            size = 2.7, vjust = 1.55) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-1, 1),
                       breaks = seq(-1, 1, 0.5), name = "correlation") +
  scale_color_identity() +
  labs(x = NULL, y = NULL,
       title = "(B) Pair-level innovation correlations\n(posterior median, 90% CrI)") +
  theme_report +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank(),
        legend.position = "right")

combined_corr <- (p_typr | p_heat) + plot_layout(widths = c(0.95, 1.5))
ggsave(file.path(FIG_DIR, "fig_typical_and_pairs.png"), combined_corr,
       width = 13, height = 5.0, dpi = 600)

#  Per-core scalar forest plot
make_post_summary <- function(name) {
  m <- post_array(name)  
  tibble(core = core_levels,
         med  = apply(m, 2, median),
         lo   = apply(m, 2, quantile, 0.05),
         hi   = apply(m, 2, quantile, 0.95)) |>
    mutate(core = factor(core, levels = rev(core_levels)))
}

scalar_summaries <- bind_rows(
  make_post_summary("mu")         |> mutate(param = "mu"),
  make_post_summary("rho")        |> mutate(param = "rho"),
  make_post_summary("sigma_proc") |> mutate(param = "sigma_proc"),
  make_post_summary("tau_obs")    |> mutate(param = "tau_obs")
) |> mutate(param = factor(param, levels = c("mu","rho","sigma_proc","tau_obs")))

p_scalar <- ggplot(scalar_summaries, aes(x = med, y = core)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.18,
                color = "#2c7fb8", orientation = "y") +
  geom_point(color = "#08519c", size = 2) +
  facet_wrap(~ param, scales = "free_x", nrow = 1) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  labs(x = "posterior median (90% CrI)", y = NULL,
       title = "Per-core scalar parameter posteriors") +
  theme_report

ggsave(file.path(FIG_DIR, "fig_scalar_forest.png"), p_scalar,
       width = 11.5, height = 3.6, dpi = 600)

# Diagnostics
summary_tbl <- posterior::summarise_draws(
  posterior::subset_draws(draws_array,
    variable = c("mu","rho","sigma_proc","tau_obs","z_mean","z_sd",
                 "typical_r","R_offdiag")),
  mean, median, sd, ~quantile(.x, c(0.05, 0.95)), default_convergence_measures()
)
# write.csv(summary_tbl, file.path(OUT_DIR, "posterior_summary.csv"), row.names = FALSE)
