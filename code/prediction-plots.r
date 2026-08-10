# ==============================================================================
# cbSCRIP Predictive Performance Plots (code/prediction-plots.r)
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(glue)
library(here)

# Set global theme if mytidyfunctions package is available
if (requireNamespace("mytidyfunctions", quietly = TRUE)) {
    mytidyfunctions::set_mytheme(text = element_text(family = "Times New Roman"))
}

# Setting names
settings_tbl <- tibble::tribble(
    ~setting, ~desc,
    1L, "Setting 1: Single effects on event 1",
    2L, "Setting 2: Single effects on both events",
    3L, "Setting 3: Opposing effects",
    4L, "Setting 4: Mixture of effects",
    5L, "Setting 5: Non-proportional hazards"
)

# Model colors (matching selection-plots.r and manuscript)
cpl_palette <- c(
    "cbSCRIP" = "#277DA1",
    "Aalen-Johansen" = "#43AA8B",
    "Penalized Fine-Gray" = "#F9C74F",
    "enet-iCox" = "#f94144",
    "Random Forest" = "#90BE6D",
    "SHBoost" = "#F9844A"
)

# Output directory
res_dir <- here::here("figs")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# Load pre-computed summaries (produced by make eval / code/02_evaluation.R)
brier_file <- here::here("results", "simulation", "brier_summary.csv")
cif_file <- here::here("results", "simulation", "cif_summary.csv")
auc_file <- here::here("results", "simulation", "auc_summary.csv")

brier_data <- if (file.exists(brier_file)) {
    read.csv(brier_file) |>
        mutate(
            model = case_match(model,
                "enet-CR" ~ "enet-iCox",
                "fastcmprsk" ~ "Penalized Fine-Gray",
                "penalized Fine-Gray method (fastcmprsk)" ~ "Penalized Fine-Gray",
                .default = as.character(model)
            ),
            title = settings_tbl$desc[match(setting, settings_tbl$setting)]
        )
} else NULL

cif_data <- if (file.exists(cif_file)) {
    read.csv(cif_file) |>
        mutate(
            Method = case_match(Method,
                "enet-CR" ~ "enet-iCox",
                "fastcmprsk" ~ "Penalized Fine-Gray",
                "penalized Fine-Gray method (fastcmprsk)" ~ "Penalized Fine-Gray",
                .default = as.character(Method)
            ),
            title = settings_tbl$desc[match(setting, settings_tbl$setting)]
        )
} else NULL

auc_data <- if (file.exists(auc_file)) {
    read.csv(auc_file) |>
        mutate(
            model = case_match(model,
                "enet-CR" ~ "enet-iCox",
                "fastcmprsk" ~ "Penalized Fine-Gray",
                "penalized Fine-Gray method (fastcmprsk)" ~ "Penalized Fine-Gray",
                .default = as.character(model)
            ),
            title = settings_tbl$desc[match(setting, settings_tbl$setting)]
        )
} else NULL

# Helper: Create placeholder when data is missing
create_empty_plot <- function(metric_name, p_val, k_val) {
    ggplot() +
        annotate(
            "text", x = 0.5, y = 0.5,
            label = glue("{metric_name}\n(p = {p_val}, k = {k_val})\n[Simulation Pending / In Progress]"),
            size = 4.5, color = "grey50", hjust = 0.5, vjust = 0.5
        ) +
        labs(
            title = glue("{metric_name}: p = {p_val}, k = {k_val}"),
            x = "Follow-up time (years)",
            y = metric_name
        ) +
        theme_minimal() +
        theme(
            panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(hjust = 0.5, face = "bold")
        )
}

# Helper: Plot Brier Score
plot_brier <- function(p_val, k_val) {
    if (is.null(brier_data)) return(create_empty_plot("Brier Score (Cause 1)", p_val, k_val))
    df <- brier_data |> filter(p == p_val, k == k_val)
    if (nrow(df) == 0) return(create_empty_plot("Brier Score (Cause 1)", p_val, k_val))
    
    df <- df |>
        mutate(
            sd_val = if_else(is.na(Brier_sd) | n_reps <= 1, 0, Brier_sd),
            ymin = pmax(0, Brier_mean - 1.96 * (sd_val / sqrt(n_reps))),
            ymax = Brier_mean + 1.96 * (sd_val / sqrt(n_reps))
        )
    
    ggplot(df, aes(x = times, y = Brier_mean, color = model)) +
        geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = model), alpha = 0.15, colour = NA) +
        geom_line(linewidth = 0.7) +
        facet_wrap(~ title, nrow = 2, scales = "free_x") +
        scale_color_manual(values = cpl_palette) +
        scale_fill_manual(values = cpl_palette) +
        labs(
            x = "Follow-up time (years)",
            y = "Brier Score (Cause 1)",
            color = "Model",
            fill = "Model"
        )
}

# Helper: Plot CIF Curves
plot_cif <- function(p_val, k_val) {
    if (is.null(cif_data)) return(create_empty_plot("Absolute Risk (CIF)", p_val, k_val))
    df <- cif_data |> filter(p == p_val, k == k_val)
    if (nrow(df) == 0) return(create_empty_plot("Absolute Risk (CIF)", p_val, k_val))
    
    df <- df |>
        mutate(
            sd_val = if_else(is.na(Risk_sd) | n_reps <= 1, 0, Risk_sd),
            ymin = pmax(0, Risk_mean - 1.96 * (sd_val / sqrt(n_reps))),
            ymax = pmin(1, Risk_mean + 1.96 * (sd_val / sqrt(n_reps)))
        )
    
    ggplot(df, aes(x = Time, y = Risk_mean, color = Method)) +
        geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = Method), alpha = 0.15, colour = NA) +
        geom_line(linewidth = 0.7) +
        facet_wrap(~ title, nrow = 2, scales = "free_x") +
        scale_color_manual(values = cpl_palette) +
        scale_fill_manual(values = cpl_palette) +
        labs(
            x = "Follow-up time (years)",
            y = "Absolute Risk (CIF)",
            color = "Model",
            fill = "Model"
        )
}

# Helper: Plot Time-dependent AUC
plot_auc <- function(p_val, k_val) {
    if (is.null(auc_data)) return(create_empty_plot("Time-dependent AUC", p_val, k_val))
    df <- auc_data |> filter(p == p_val, k == k_val)
    if (nrow(df) == 0) return(create_empty_plot("Time-dependent AUC", p_val, k_val))
    
    df <- df |>
        mutate(
            sd_val = if_else(is.na(AUC_sd) | n_reps <= 1, 0, AUC_sd),
            ymin = pmax(0, AUC_mean - 1.96 * (sd_val / sqrt(n_reps))),
            ymax = pmin(1, AUC_mean + 1.96 * (sd_val / sqrt(n_reps)))
        )
    
    ggplot(df, aes(x = times, y = AUC_mean, color = model)) +
        geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = model), alpha = 0.15, colour = NA) +
        geom_line(linewidth = 0.7) +
        facet_wrap(~ title, nrow = 2, scales = "free_x") +
        scale_color_manual(values = cpl_palette) +
        scale_fill_manual(values = cpl_palette) +
        labs(
            x = "Follow-up time (years)",
            y = "Time-dependent AUC",
            color = "Model",
            fill = "Model"
        )
}

# Helper to save plot
save_plot <- function(filepath, plot_obj) {
    if (!is.null(plot_obj)) {
        ggsave(filepath, plot_obj, width = 200, height = 120, units = "mm", dpi = 300, bg = "transparent")
        cli::cli_alert_success(paste0("Generated: ", basename(filepath)))
    }
}

# Standard parameter combinations to generate
target_configs <- tribble(
    ~p, ~k,
    120, 20,
    500, 20,
    500, 84,
    1000, 20,
    1000, 168
)

# Also include any other combinations present in data
all_configs <- distinct(bind_rows(
    target_configs,
    if (!is.null(brier_data)) select(brier_data, p, k) else NULL,
    if (!is.null(cif_data)) select(cif_data, p, k) else NULL,
    if (!is.null(auc_data)) select(auc_data, p, k) else NULL
)) |> filter(!is.na(p) & !is.na(k)) |> distinct()

cli::cli_alert_info("Generating simulation plots...")

for (i in seq_len(nrow(all_configs))) {
    p_val <- all_configs$p[i]
    k_val <- all_configs$k[i]
    
    # 1. Brier Plot
    p_brier <- plot_brier(p_val, k_val)
    if (!is.null(p_brier)) {
        save_plot(file.path(res_dir, glue("brier-p{p_val}-k{k_val}.png")), p_brier)
    }
    
    # 2. CIF Plot
    p_cif <- plot_cif(p_val, k_val)
    if (!is.null(p_cif)) {
        save_plot(file.path(res_dir, glue("cif-p{p_val}-k{k_val}.png")), p_cif)
    }
    
    # 3. AUC Plot
    p_auc <- plot_auc(p_val, k_val)
    if (!is.null(p_auc)) {
        save_plot(file.path(res_dir, glue("auc-p{p_val}-k{k_val}.png")), p_auc)
    }
}

cli::cli_alert_info("Simulation plot generation complete.")
