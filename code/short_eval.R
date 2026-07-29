# ==============================================================================
# Short Evaluation Script (code/short_eval.R)
# 
# Performs a quick evaluation of the custom simulation runs by printing
# summary statistics directly to the console and generating a performance plot.
# ==============================================================================

library(ggplot2)
library(dplyr)
library(purrr)
library(glue)
library(here)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
settings <- if (length(args) >= 1 && nzchar(args[1])) as.numeric(strsplit(args[1], ",")[[1]]) else c(1)
p_vals <- if (length(args) >= 2 && nzchar(args[2])) as.numeric(strsplit(args[2], ",")[[1]]) else c(120)
reps <- if (length(args) >= 3 && nzchar(args[3])) as.numeric(args[3]) else 5
k_vals <- if (length(args) >= 4 && nzchar(args[4])) as.numeric(strsplit(args[4], ",")[[1]]) else NULL

varsel_dir <- here("results", "simulation", "varsel")
preds_dir <- here("results", "simulation", "preds")
benchmark_dir <- here("figs")
dir.create(benchmark_dir, recursive = TRUE, showWarnings = FALSE)

cli::cli_alert_info(paste0("\n=========================================="))
cli::cli_alert_info(paste0("   Running Short Simulation Evaluation    "))
cli::cli_alert_info(paste0("==========================================\n"))

# 1. Variable Selection Evaluation
varsel_files <- list.files(varsel_dir, pattern = "^varsel_s.*_p.*_rep.*\\.rds$", full.names = TRUE)
matching_varsel <- varsel_files[sapply(basename(varsel_files), function(fn) {
    s <- as.numeric(sub(".*_s([0-9]+)_.*", "\\1", fn))
    p <- as.numeric(sub(".*_p([0-9]+)_.*", "\\1", fn))
    k <- as.numeric(sub(".*_k([0-9]+)_.*", "\\1", fn))
    r <- as.numeric(sub(".*_rep([0-9]+)\\.rds", "\\1", fn))
    k_match <- if (is.null(k_vals)) TRUE else (k %in% k_vals)
    (s %in% settings) && (p %in% p_vals) && k_match && (r <= reps)
})]

if (length(matching_varsel) > 0) {
    varsel_results <- map_df(matching_varsel, readRDS)
    
    varsel_summary <- varsel_results %>%
        group_by(model, p, k, setting) %>%
        summarize(
            Mean_Sensitivity = mean(Sensitivity, na.rm = TRUE),
            Mean_Specificity = mean(Specificity, na.rm = TRUE),
            Mean_CSR = mean(exact_match, na.rm = TRUE),
            .groups = "drop"
        )
    
    cli::cli_alert_info(paste0("--- Variable Selection Performance ---"))
    print(varsel_summary)
    cli::cli_alert_info(paste0(""))
    
    # Generate Sensitivity vs Model Size plot
    plot_df <- varsel_results %>%
        group_by(model, model_size) %>%
        summarize(Mean_Sensitivity = mean(Sensitivity, na.rm = TRUE), .groups = "drop")
    
    cbPalette <- c("cbSCRIP" = "#56B4E9", "enet-iCox" = "#D55E00", "fastcmprsk" = "#CC79A7")
    
    p_varsel <- ggplot(plot_df, aes(x = model_size, y = Mean_Sensitivity, color = model)) +
        geom_line(linewidth = 1) +
        geom_point() +
        theme_bw() +
        scale_color_manual(values = cbPalette) +
        labs(
            title = "Custom Run: Variable Selection Sensitivity",
            subtitle = glue("Settings: {paste(settings, collapse=',')} | P: {paste(p_vals, collapse=',')} | Reps: {reps}"),
            x = "Selected Model Size (number of non-zero coefficients)",
            y = "Sensitivity (TPR)",
            color = "Method"
        )
    
    plot_path <- file.path(benchmark_dir, "short_varsel_plot.png")
    ggsave(plot_path, p_varsel, width = 6, height = 4, dpi = 150)
    cli::cli_alert_info(paste0(glue("Saved custom run performance plot to: {plot_path}\n")))
} else {
    cli::cli_alert_info(paste0("No matching variable selection files found.\n"))
}

# 2. Timing Aggregation
pred_files <- list.files(preds_dir, pattern = "^models_s.*_p.*_rep.*\\.rds$", full.names = TRUE)
matching_preds <- pred_files[sapply(basename(pred_files), function(fn) {
    s <- as.numeric(sub(".*_s([0-9]+)_.*", "\\1", fn))
    p <- as.numeric(sub(".*_p([0-9]+)_.*", "\\1", fn))
    k <- as.numeric(sub(".*_k([0-9]+)_.*", "\\1", fn))
    r <- as.numeric(sub(".*_rep([0-9]+)\\.rds", "\\1", fn))
    k_match <- if (is.null(k_vals)) TRUE else (k %in% k_vals)
    (s %in% settings) && (p %in% p_vals) && k_match && (r <= reps)
})]

if (length(matching_preds) > 0) {
    timing_list <- list()
    for (f in matching_preds) {
        try({
            bundle <- readRDS(f)
            if (!is.null(bundle$timing)) {
                df <- tibble(
                    method = names(bundle$timing),
                    elapsed = as.numeric(bundle$timing),
                    setting = bundle$meta$setting,
                    p = bundle$meta$p
                )
                timing_list[[length(timing_list) + 1]] <- df
            }
        }, silent = TRUE)
    }
    
    if (length(timing_list) > 0) {
        timing_data <- bind_rows(timing_list)
        timing_summary <- timing_data %>%
            group_by(method, p) %>%
            summarize(Mean_Fit_Time_Sec = mean(elapsed, na.rm = TRUE), .groups = "drop")
        
        cli::cli_alert_info(paste0("--- Computational Timing Performance ---"))
        print(timing_summary)
        cli::cli_alert_info(paste0(""))
    }
} else {
    cli::cli_alert_info(paste0("No matching timing data files found.\n"))
}

cli::cli_alert_info(paste0("=========================================="))
