# ==============================================================================
# cbSCRIP Empirical Timing Benchmark (03_benchmarks.R)
# 
# Reads and aggregates the recorded timing information directly from the
# saved simulation results under paper/results/preds/
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(glue)
library(here)

save_dir_preds <- here("results", "simulation", "preds")
benchmark_dir <- here("results", "simulation", "benchmark")
dir.create(benchmark_dir, recursive = TRUE, showWarnings = FALSE)

# Find all saved model RDS files
pred_files <- list.files(save_dir_preds, pattern = "^models_s.*_p.*_rep.*\\.rds$", full.names = TRUE)

if (length(pred_files) == 0) {
    stop("No simulation output files found in paper/results/preds/. Please run simulations first.")
}

message(glue("Found {length(pred_files)} simulation result files. Aggregating timing info..."))

# Extract timing data
timing_list <- list()
for (file_path in pred_files) {
    try({
        bundle <- readRDS(file_path)
        if (!is.null(bundle$timing)) {
            df <- tibble(
                method = names(bundle$timing),
                elapsed = as.numeric(bundle$timing),
                setting = bundle$meta$setting,
                p = bundle$meta$p,
                rep = bundle$meta$seed
            )
            timing_list[[length(timing_list) + 1]] <- df
        }
    })
}

if (length(timing_list) == 0) {
    stop("No timing information found inside the RDS files.")
}

timing_data <- bind_rows(timing_list)

# Aggregate timing results across settings/replicates
timing_summary <- timing_data %>%
    group_by(method, p) %>%
    summarize(
        mean_time = mean(elapsed, na.rm = TRUE),
        sd_time = sd(elapsed, na.rm = TRUE),
        n = n(),
        se_time = sd_time / sqrt(n),
        lower_time = mean_time - 1.96 * se_time,
        upper_time = mean_time + 1.96 * se_time,
        .groups = "drop"
    )

# Print summary table
print(timing_summary)

# Define matching color palette
cbPalette <- c("cbSCRIP" = "#56B4E9", "enet-iCox" = "#D55E00", "fastcmprsk" = "#CC79A7", "rfsrc" = "#E69F00")

# Plot Results (Mean with Standard Error Bars)
plot_benchmark <- ggplot(timing_summary, aes(x = as.factor(p), y = mean_time, fill = method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = pmax(0, mean_time - se_time), ymax = mean_time + se_time), 
                  position = position_dodge(width = 0.9), width = 0.25) +
    theme_bw() +
    scale_fill_manual(values = cbPalette) +
    labs(
        title = "Empirical Model Fitting Execution Time",
        subtitle = "Averaged across simulation settings and replications (N = 300)",
        x = "Number of Predictors (p)",
        y = "Mean Computation Time (Seconds)",
        fill = "Method"
    ) +
    theme(
        plot.title = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12)
    )

# Save Plot & Data
ggsave(here("figs", "timing_benchmark.png"), plot_benchmark, width = 8, height = 5, dpi = 300)
saveRDS(timing_summary, file.path(benchmark_dir, "timing_results.rds"))
message("Timing benchmark plot and results successfully saved.")
