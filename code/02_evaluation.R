# ==============================================================================
# cbSCRIP Evaluation and Plotting Pipeline (02_evaluation.R)
#
# Aggregates pre-computed metrics from results/simulation/
# Output: results/simulation/{varsel,brier,cif}_summary.csv
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(here)

varsel_dir  <- here("results", "simulation", "varsel")
preds_dir   <- here("results", "simulation", "preds")
results_dir <- here("results", "simulation")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. VARIABLE SELECTION SUMMARY
# ==============================================================================
cli::cli_alert_info(paste0("=== Variable Selection ==="))

varsel_files <- list.files(varsel_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(varsel_files) > 0) {
    varsel_list <- lapply(varsel_files, function(fp) {
        df <- tryCatch(readRDS(fp), error = function(e) NULL)
        if (is.null(df) || nrow(df) == 0) return(NULL)
        if (!"k" %in% colnames(df)) {
            fn <- basename(fp)
            if (grepl("_k([0-9]+)_", fn)) {
                df$k <- as.numeric(gsub(".*_k([0-9]+)_.*", "\\1", fn))
            } else {
                df$k <- case_match(df$p[1], 120 ~ 20, 500 ~ 84, 1000 ~ 168, .default = 20)
            }
        }
        df
    })
    varsel_all <- bind_rows(Filter(Negate(is.null), varsel_list))
    varsel_summary <- varsel_all |>
        group_by(model, setting, p, k, lambda) |>
        summarise(
            mean_sens       = mean(Sensitivity, na.rm = TRUE),
            sd_sens         = sd(Sensitivity,   na.rm = TRUE),
            n               = n(),
            se_sens         = sd_sens / sqrt(n),
            lower_sens      = mean_sens - 1.96 * se_sens,
            upper_sens      = mean_sens + 1.96 * se_sens,
            mean_spec       = mean(Specificity,  na.rm = TRUE),
            mean_model_size = mean(model_size,   na.rm = TRUE),
            csr             = mean(exact_match,  na.rm = TRUE),
            .groups = "drop"
        )
    write.csv(varsel_summary, file.path(results_dir, "varsel_summary.csv"),
              row.names = FALSE)
    cli::cli_alert_info(paste0(sprintf("Saved varsel_summary.csv  [%d rows]", nrow(varsel_summary))))

    # --- Correct Selection Rate (R2: exact recovery of the true support) ---
    # csr_path = proportion of replicates where SOME lambda on the path recovers
    # the true support exactly. This is an oracle/best-case rate: it assumes
    # lambda is chosen perfectly, so it upper-bounds what CV would achieve.
    #
    # A "best fixed lambda" version is NOT computable here: lambda_max is derived
    # from the data, so each replicate has its own grid (589 distinct lambdas
    # across 20 replicates for s1/p120/k20, only 11 shared by >1 replicate).
    # Grouping by lambda across replicates therefore yields groups of size ~1,
    # and any statistic over them describes a single replicate, not a rate.
    # The same caveat applies to the `csr` column of varsel_summary above.
    csr_summary <- varsel_all |>
        group_by(model, setting, p, k, rep) |>
        summarise(hit = any(exact_match, na.rm = TRUE), .groups = "drop") |>
        group_by(model, setting, p, k) |>
        summarise(csr_path = mean(hit), n_reps = n(), .groups = "drop")

    write.csv(csr_summary, file.path(results_dir, "csr_summary.csv"),
              row.names = FALSE)
    cli::cli_alert_info(paste0(sprintf("Saved csr_summary.csv  [%d rows]", nrow(csr_summary))))
} else {
    cli::cli_alert_info(paste0("No varsel files found in ", varsel_dir))
}


# ==============================================================================
# 2. BRIER SCORE SUMMARY  (pre-computed by 01_simulations.R)
# ==============================================================================
cli::cli_alert_info(paste0("\n=== Brier Score ==="))

preds_files <- list.files(preds_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(preds_files) > 0) {
    brier_list <- lapply(preds_files, function(fp) {
        res <- readRDS(fp)
        if (is.null(res$brier_table)) {
            cli::cli_alert_info(paste0("  Skipping ", basename(fp), " (no brier_table)"))
            return(NULL)
        }
        bt <- res$brier_table
        if (!"k" %in% colnames(bt)) {
            k_val <- if (!is.null(res$meta$k)) res$meta$k else case_match(res$meta$p, 120 ~ 20, 500 ~ 84, 1000 ~ 168, .default = 20)
            if (is.null(k_val)) {
                p_val <- if (!is.null(bt$p)) bt$p[1] else 120
                k_val <- case_match(p_val, 120 ~ 20, 500 ~ 84, 1000 ~ 168, .default = 20)
            }
            bt$k <- k_val
        }
        bt
    })
    brier_list <- Filter(Negate(is.null), brier_list)

    if (length(brier_list) > 0) {
        brier_all <- bind_rows(brier_list)
        # Interpolate onto a common grid of 100 points per setting, p, k
        brier_interpolated <- brier_all |>
            group_by(setting, p, k) |>
            group_split() |>
            lapply(function(sub_df) {
                t_min <- min(sub_df$times, na.rm = TRUE)
                t_max <- max(sub_df$times, na.rm = TRUE)
                if (t_min == t_max) t_max <- t_min + 1
                common_times <- seq(t_min, t_max, length.out = 100)
                
                sub_df |>
                    group_by(model, rep) |>
                    group_split() |>
                    lapply(function(rep_df) {
                        if (nrow(rep_df) < 2) return(NULL)
                        interp_val <- tryCatch({
                            approx(x = rep_df$times, y = rep_df$Brier, xout = common_times, rule = 2)$y
                        }, error = function(e) rep(NA, 100))
                        tibble(
                            model   = rep_df$model[1],
                            times   = common_times,
                            Brier   = interp_val,
                            rep     = rep_df$rep[1],
                            setting = rep_df$setting[1],
                            p       = rep_df$p[1],
                            k       = rep_df$k[1]
                        )
                    }) |>
                    bind_rows()
            }) |>
            bind_rows()

        brier_summary <- brier_interpolated |>
            group_by(model, times, setting, p, k) |>
            summarise(
                Brier_mean = mean(Brier, na.rm = TRUE),
                Brier_sd   = sd(Brier,   na.rm = TRUE),
                n_reps     = sum(!is.na(Brier)),
                .groups = "drop"
            )
        write.csv(brier_summary, file.path(results_dir, "brier_summary.csv"),
                  row.names = FALSE)
        cli::cli_alert_info(paste0(sprintf("Saved brier_summary.csv  [%d rows, settings=%s, p=%s, k=%s]",
                        nrow(brier_summary),
                        paste(sort(unique(brier_summary$setting)), collapse = ","),
                        paste(sort(unique(brier_summary$p)),       collapse = ","),
                        paste(sort(unique(brier_summary$k)),       collapse = ","))))
    } else {
        cli::cli_alert_info(paste0("No Brier tables found (run simulations first or check 01_simulations.R)"))
    }
} else {
    cli::cli_alert_info(paste0("No prediction files found in ", preds_dir))
}


# ==============================================================================
# 3. CIF CURVE SUMMARY  (pre-computed by 01_simulations.R)
# ==============================================================================
cli::cli_alert_info(paste0("\n=== CIF Curves ==="))

if (length(preds_files) > 0) {
    cif_list <- lapply(preds_files, function(fp) {
        res <- readRDS(fp)
        if (is.null(res$cif_table)) {
            cli::cli_alert_info(paste0("  Skipping ", basename(fp), " (no cif_table)"))
            return(NULL)
        }
        ct <- res$cif_table
        if (!"k" %in% colnames(ct)) {
            k_val <- if (!is.null(res$meta$k)) res$meta$k else case_match(res$meta$p, 120 ~ 20, 500 ~ 84, 1000 ~ 168, .default = 20)
            if (is.null(k_val)) {
                p_val <- if (!is.null(ct$p)) ct$p[1] else 120
                k_val <- case_match(p_val, 120 ~ 20, 500 ~ 84, 1000 ~ 168, .default = 20)
            }
            ct$k <- k_val
        }
        ct
    })
    cif_list <- Filter(Negate(is.null), cif_list)

    if (length(cif_list) > 0) {
        cif_all <- bind_rows(cif_list)
        # Interpolate onto a common grid of 100 points per setting, p, k
        cif_interpolated <- cif_all |>
            group_by(setting, p, k) |>
            group_split() |>
            lapply(function(sub_df) {
                t_min <- min(sub_df$Time, na.rm = TRUE)
                t_max <- max(sub_df$Time, na.rm = TRUE)
                if (t_min == t_max) t_max <- t_min + 1
                common_times <- seq(t_min, t_max, length.out = 100)
                
                sub_df |>
                    group_by(Method, rep) |>
                    group_split() |>
                    lapply(function(rep_df) {
                        if (nrow(rep_df) < 2) return(NULL)
                        interp_val <- tryCatch({
                            approx(x = rep_df$Time, y = rep_df$Risk, xout = common_times, rule = 2)$y
                        }, error = function(e) rep(NA, 100))
                        tibble(
                            Method  = rep_df$Method[1],
                            Time    = common_times,
                            Risk    = interp_val,
                            rep     = rep_df$rep[1],
                            setting = rep_df$setting[1],
                            p       = rep_df$p[1],
                            k       = rep_df$k[1]
                        )
                    }) |>
                    bind_rows()
            }) |>
            bind_rows()

        cif_summary <- cif_interpolated |>
            group_by(Method, Time, setting, p, k) |>
            summarise(
                Risk_mean = mean(Risk, na.rm = TRUE),
                Risk_sd   = sd(Risk,   na.rm = TRUE),
                n_reps    = sum(!is.na(Risk)),
                .groups = "drop"
            )
        write.csv(cif_summary, file.path(results_dir, "cif_summary.csv"),
                  row.names = FALSE)
        cli::cli_alert_info(paste0(sprintf("Saved cif_summary.csv  [%d rows, settings=%s, p=%s, k=%s]",
                        nrow(cif_summary),
                        paste(sort(unique(cif_summary$setting)), collapse = ","),
                        paste(sort(unique(cif_summary$p)),       collapse = ","),
                        paste(sort(unique(cif_summary$k)),       collapse = ","))))
    } else {
        cli::cli_alert_info(paste0("No CIF tables found (run simulations first or check 01_simulations.R)"))
    }
}

# ==============================================================================
# 4. TIME-DEPENDENT AUC SUMMARY (pre-computed by 01_simulations.R)
# ==============================================================================
cli::cli_alert_info(paste0("\n=== Time-dependent AUC ==="))

if (length(preds_files) > 0) {
    auc_list <- lapply(preds_files, function(fp) {
        res <- readRDS(fp)
        if (is.null(res$auc_table)) {
            cli::cli_alert_info(paste0("  Skipping ", basename(fp), " (no auc_table)"))
            return(NULL)
        }
        at <- res$auc_table
        if (!"k" %in% colnames(at)) {
            k_val <- if (!is.null(res$meta$k)) res$meta$k else case_match(res$meta$p, 120 ~ 20, 500 ~ 84, 1000 ~ 168, .default = 20)
            if (is.null(k_val)) {
                p_val <- if (!is.null(at$p)) at$p[1] else 120
                k_val <- case_match(p_val, 120 ~ 20, 500 ~ 84, 1000 ~ 168, .default = 20)
            }
            at$k <- k_val
        }
        at
    })
    auc_list <- Filter(Negate(is.null), auc_list)

    if (length(auc_list) > 0) {
        auc_all <- bind_rows(auc_list)
        # Interpolate onto a common grid of 100 points per setting, p, k
        auc_interpolated <- auc_all |>
            group_by(setting, p, k) |>
            group_split() |>
            lapply(function(sub_df) {
                t_min <- min(sub_df$times, na.rm = TRUE)
                t_max <- max(sub_df$times, na.rm = TRUE)
                if (t_min == t_max) t_max <- t_min + 1
                common_times <- seq(t_min, t_max, length.out = 100)
                
                sub_df |>
                    group_by(model, rep) |>
                    group_split() |>
                    lapply(function(rep_df) {
                        if (nrow(rep_df) < 2) return(NULL)
                        interp_val <- tryCatch({
                            approx(x = rep_df$times, y = rep_df$AUC, xout = common_times, rule = 2)$y
                        }, error = function(e) rep(NA, 100))
                        tibble(
                            model   = rep_df$model[1],
                            times   = common_times,
                            AUC     = interp_val,
                            rep     = rep_df$rep[1],
                            setting = rep_df$setting[1],
                            p       = rep_df$p[1],
                            k       = rep_df$k[1]
                        )
                    }) |>
                    bind_rows()
            }) |>
            bind_rows()

        auc_summary <- auc_interpolated |>
            group_by(model, times, setting, p, k) |>
            summarise(
                AUC_mean = mean(AUC, na.rm = TRUE),
                AUC_sd   = sd(AUC,   na.rm = TRUE),
                n_reps   = sum(!is.na(AUC)),
                .groups = "drop"
            )
        write.csv(auc_summary, file.path(results_dir, "auc_summary.csv"),
                  row.names = FALSE)
        cli::cli_alert_info(paste0(sprintf("Saved auc_summary.csv  [%d rows, settings=%s, p=%s, k=%s]",
                        nrow(auc_summary),
                        paste(sort(unique(auc_summary$setting)), collapse = ","),
                        paste(sort(unique(auc_summary$p)),       collapse = ","),
                        paste(sort(unique(auc_summary$k)),       collapse = ","))))
    } else {
        cli::cli_alert_info(paste0("No AUC tables found (run simulations first or check 01_simulations.R)"))
    }
}

# ==============================================================================
# 5. FIT TIMING SUMMARY (pre-computed by 01_simulations.R)
# ==============================================================================
# Timings are only ever compared within a stage. "path" is a bare regularization
# path used for the variable-selection comparison; "tuned" includes the
# cross-validation needed to get a prediction-ready model. Mixing the two would
# compare a no-CV fit against a 5-fold CV fit.
cli::cli_alert_info(paste0("\n=== Fit timings ==="))

if (length(preds_files) > 0) {
    timing_list <- lapply(preds_files, function(fp) {
        res <- readRDS(fp)
        tm <- res$timing
        # Bundles written before timings were grouped by stage carry a flat
        # list; they are not comparable and are skipped rather than guessed at.
        if (is.null(tm) || !all(c("path", "tuned") %in% names(tm))) return(NULL)

        stage_rows <- function(stage) {
            vals <- tm[[stage]]
            tibble(
                model   = names(vals),
                stage   = stage,
                seconds = as.numeric(unlist(lapply(vals, function(x) if (is.null(x)) NA_real_ else x))),
                rep     = res$meta$seed,
                setting = res$meta$setting,
                p       = res$meta$p,
                k       = res$meta$k
            )
        }
        bind_rows(stage_rows("path"), stage_rows("tuned"))
    })
    timing_list <- Filter(Negate(is.null), timing_list)

    if (length(timing_list) > 0) {
        timing_summary <- bind_rows(timing_list) |>
            group_by(model, stage, setting, p, k) |>
            summarise(
                seconds_mean = mean(seconds, na.rm = TRUE),
                seconds_sd   = sd(seconds,   na.rm = TRUE),
                n_reps       = sum(!is.na(seconds)),
                .groups = "drop"
            )
        write.csv(timing_summary, file.path(results_dir, "timing_summary.csv"),
                  row.names = FALSE)
        cli::cli_alert_info(paste0(sprintf("Saved timing_summary.csv  [%d rows, stages=%s]",
                        nrow(timing_summary),
                        paste(sort(unique(timing_summary$stage)), collapse = ","))))
    } else {
        cli::cli_alert_info(paste0("No stage-grouped timings found. Bundles predating the ",
                                   "timing fix in 01_simulations.R carry a flat, non-comparable ",
                                   "timing list; rerun the simulations to populate this."))
    }
}

cli::cli_alert_info(paste0("\nEvaluation complete."))
