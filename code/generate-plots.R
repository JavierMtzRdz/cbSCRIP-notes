# ==============================================================================
# cbSCRIP Simulation Plot Generation Helper Script (code/generate-plots.R)
# ==============================================================================

if (!require("pacman")) install.packages("pacman")
if (!require("remotes")) install.packages("remotes")

pacman::p_load(tidyverse, here, glue)

# Attempt to install mytidyfunctions from GitHub if not found
if (!requireNamespace("mytidyfunctions", quietly = TRUE))
    remotes::install_github("JavierMtzRdz/mytidyfunctions")

# Load custom fonts and set theme only if the package is available
if (requireNamespace("mytidyfunctions", quietly = TRUE)) {
    mytidyfunctions::set_mytheme(text = element_text(family = "Times New Roman"))
}

# Setting names
settings_tbl <- tibble::tribble(
    ~setting, ~desc,
    1L, "1. Single effects on end-point of interest",
    2L, "2. Single effects on both end-points",
    3L, "3. Opposing effects",
    4L, "4. Mixture of effects",
    5L, "5. Non-proportional hazards"
)

# Model colors
cpl_palette <- c(
    "cbSCRIP"             = "#277DA1",
    "Aalen-Johansen"      = "#43AA8B",
    "Penalized Fine-Gray" = "#F9C74F",
    "enet-iCox"           = "#f94144",
    "Random Forest"       = "#90BE6D",
    "SHBoost"             = "#F9844A"
)

# Load pre-computed summaries (produced by make eval / code/02_evaluation.R)
brier_file <- here::here("results", "simulation", "brier_summary.csv")
cif_file   <- here::here("results", "simulation", "cif_summary.csv")
auc_file   <- here::here("results", "simulation", "auc_summary.csv")

brier_data <- if (file.exists(brier_file)) {
    read.csv(brier_file) |>
        mutate(title = settings_tbl$desc[match(setting, settings_tbl$setting)])
} else {
    stop("brier_summary.csv not found. Run: make eval")
}

cif_data <- if (file.exists(cif_file)) {
    read.csv(cif_file) |>
        mutate(title = settings_tbl$desc[match(setting, settings_tbl$setting)])
} else {
    stop("cif_summary.csv not found. Run: make eval")
}

auc_data <- if (file.exists(auc_file)) {
    read.csv(auc_file) |>
        mutate(title = settings_tbl$desc[match(setting, settings_tbl$setting)])
} else {
    stop("auc_summary.csv not found. Run: make eval")
}

# Helper: filter and plot Brier
plot_brier <- function(p_val, k_val) {
    df <- brier_data |> filter(p == p_val, k == k_val)
    if (is.null(df) || nrow(df) == 0) {
        cli::cli_alert_info(paste0(glue("No Brier data for p={p_val}, k={k_val}"))); return(invisible(NULL))
    }
    df |>
        ggplot(aes(times, Brier_mean, colour = model)) +
        geom_ribbon(aes(ymin = Brier_mean - 1.96 * (Brier_sd / sqrt(n_reps)),
                        ymax = Brier_mean + 1.96 * (Brier_sd / sqrt(n_reps)),
                        fill = model), alpha = 0.15, colour = NA) +
        geom_line(aes(linewidth = ifelse(model == "cbSCRIP", 0.5, 0.1))) +
        facet_wrap(. ~ title, nrow = 2, scales = "free_x") +
        scale_colour_manual(values = cpl_palette) +
        scale_fill_manual(values = cpl_palette) +
        scale_linewidth_continuous(range = c(0.5, 0.8), guide = "none") +
        labs(x = "Follow-up time (years)", colour = "Models", fill = "Models",
             y = "Brier Score for Cause 1 Predictions",
             title = glue("p = {p_val}, k = {k_val}"))
}

# Helper: filter and plot CIF
plot_cif <- function(p_val, k_val) {
    df <- cif_data |> filter(p == p_val, k == k_val)
    if (is.null(df) || nrow(df) == 0) {
        cli::cli_alert_info(paste0(glue("No CIF data for p={p_val}, k={k_val}"))); return(invisible(NULL))
    }
    df |>
        ggplot(aes(Time, Risk_mean, colour = Method)) +
        geom_ribbon(aes(ymin = Risk_mean - 1.96 * (Risk_sd / sqrt(n_reps)),
                        ymax = Risk_mean + 1.96 * (Risk_sd / sqrt(n_reps)),
                        fill = Method), alpha = 0.15, colour = NA) +
        geom_line(aes(linewidth = ifelse(Method == "cbSCRIP", 0.5, 0.1))) +
        facet_wrap(. ~ title, nrow = 2, scales = "free_x") +
        scale_colour_manual(values = cpl_palette) +
        scale_fill_manual(values = cpl_palette) +
        scale_linewidth_continuous(range = c(0.5, 0.8), guide = "none") +
        labs(x = "Follow-up time (years)", colour = "Models", fill = "Models",
             y = "Absolute Risk (CIF)",
             title = glue("p = {p_val}, k = {k_val}"))
}

# Helper: filter and plot AUC
plot_auc <- function(p_val, k_val) {
    df <- auc_data |> filter(p == p_val, k == k_val)
    if (is.null(df) || nrow(df) == 0) {
        cli::cli_alert_info(paste0(glue("No AUC data for p={p_val}, k={k_val}"))); return(invisible(NULL))
    }
    df |>
        ggplot(aes(times, AUC_mean, colour = model)) +
        geom_ribbon(aes(ymin = AUC_mean - 1.96 * (AUC_sd / sqrt(n_reps)),
                        ymax = AUC_mean + 1.96 * (AUC_sd / sqrt(n_reps)),
                        fill = model), alpha = 0.15, colour = NA) +
        geom_line(aes(linewidth = ifelse(model == "cbSCRIP", 0.5, 0.1))) +
        facet_wrap(. ~ title, nrow = 2, scales = "free_x") +
        scale_colour_manual(values = cpl_palette) +
        scale_fill_manual(values = cpl_palette) +
        scale_linewidth_continuous(range = c(0.5, 0.8), guide = "none") +
        labs(x = "Follow-up time (years)", colour = "Models", fill = "Models",
             y = "Time-dependent AUC",
             title = glue("p = {p_val}, k = {k_val}"))
}

# Output directory
res_dir <- here::here("figs")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# Helper to safely save plots only when they are not NULL
save_plot <- function(filepath, plot_obj) {
    if (!is.null(plot_obj)) {
        ggsave(filepath, plot_obj, width = 20, height = 12, units = "cm", dpi = 300)
        cli::cli_alert_info(paste0("Saved: ", basename(filepath)))
    } else {
        cli::cli_alert_info(paste0("Skipped (no data): ", basename(filepath)))
    }
}

# Generate and save Brier plots
cli::cli_alert_info(paste0("Generating Brier plots..."))
save_plot(file.path(res_dir, "brier_grid_N300_p120_k20.png"), plot_brier(120, 20))
save_plot(file.path(res_dir, "brier_grid_N300_p500_k20.png"), plot_brier(500, 20))
save_plot(file.path(res_dir, "brier_grid_N300_p500_k84.png"), plot_brier(500, 84))
save_plot(file.path(res_dir, "brier_grid_N300_p1000_k20.png"), plot_brier(1000, 20))
save_plot(file.path(res_dir, "brier_grid_N300_p1000_k168.png"), plot_brier(1000, 168))

# Generate and save CIF plots
cli::cli_alert_info(paste0("Generating CIF plots..."))
save_plot(file.path(res_dir, "cif_grid_N300_p120_k20.png"), plot_cif(120, 20))
save_plot(file.path(res_dir, "cif_grid_N300_p500_k20.png"), plot_cif(500, 20))
save_plot(file.path(res_dir, "cif_grid_N300_p500_k84.png"), plot_cif(500, 84))
save_plot(file.path(res_dir, "cif_grid_N300_p1000_k20.png"), plot_cif(1000, 20))
save_plot(file.path(res_dir, "cif_grid_N300_p1000_k168.png"), plot_cif(1000, 168))

# Generate and save AUC plots
cli::cli_alert_info(paste0("Generating AUC plots..."))
save_plot(file.path(res_dir, "auc_grid_N300_p120_k20.png"), plot_auc(120, 20))
save_plot(file.path(res_dir, "auc_grid_N300_p500_k20.png"), plot_auc(500, 20))
save_plot(file.path(res_dir, "auc_grid_N300_p500_k84.png"), plot_auc(500, 84))
save_plot(file.path(res_dir, "auc_grid_N300_p1000_k20.png"), plot_auc(1000, 20))
save_plot(file.path(res_dir, "auc_grid_N300_p1000_k168.png"), plot_auc(1000, 168))

cli::cli_alert_info(paste0("All simulation plots successfully generated and saved to figs/."))
