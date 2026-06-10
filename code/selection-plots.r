library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(forcats)
library(glue)
library(here)
library(zoo)

# Define directories
varsel_dir <- here("results", "simulation", "varsel")
figs_dir <- here("figs")
dir.create(figs_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Find and read all varsel files
varsel_files <- list.files(varsel_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(varsel_files) == 0) {
    stop("No varsel files found in ", varsel_dir)
}

cli::cli_alert_info(paste0("Reading varsel files..."))
select_mod <- map_df(varsel_files, function(f) {
    res <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) return(NULL)
    
    if (!"k" %in% names(res)) {
        fn <- basename(f)
        if (grepl("_k([0-9]+)_", fn)) {
            res$k <- as.numeric(sub(".*_k([0-9]+)_.*", "\\1", fn))
        } else {
            p_val <- if (!is.null(res$p)) res$p[1] else 120
            res$k <- case_match(p_val, 120 ~ 20, 500 ~ 84, 1000 ~ 168, .default = 20)
        }
    }
    
    res |> mutate(
        model = as.character(model),
        rep = as.numeric(rep),
        setting = as.numeric(setting),
        p = as.numeric(p),
        k = as.numeric(k)
    )
})

# 2. Add dim label and clean model names
select_mod <- select_mod |>
    mutate(
        dim = glue("p = {p}, k = {k}"),
        model = case_match(model,
            "enet-CR" ~ "enet-iCox",
            "fastcmprsk" ~ "Penalized Fine-Gray",
            "penalized Fine-Gray method (fastcmprsk)" ~ "Penalized Fine-Gray",
            .default = model
        )
    )

# Palette matching cpl_palette
cpl_palette <- c(
    "cbSCRIP"             = "#277DA1",
    "Aalen-Johansen"      = "#43AA8B",
    "Penalized Fine-Gray" = "#F9C74F",
    "enet-iCox"           = "#f94144",
    "Random Forest"       = "#90BE6D",
    "SHBoost"             = "#F9844A"
)

# Set global theme if mytidyfunctions package is available
if (requireNamespace("mytidyfunctions", quietly = TRUE)) {
    mytidyfunctions::set_mytheme(text = element_text(family = "Times New Roman"))
}

# ==============================================================================
# Plot 1: Sensitivity (TPR) vs. Model Size (with extrapolation & 95% CI)
# ==============================================================================
cli::cli_alert_info("Interpolating Sensitivity vs. Model Size...")

sum_tpr <- select_mod |>
    group_by(model, setting, dim, rep, p, k) |>
    group_split() |>
    map_dfr(~ {
        max_p <- .x$p[1]
        
        # Deduplicate model_size per replicate by taking maximum sensitivity, minimum specificity, maximum MCC
        clean_rep <- .x |>
            group_by(model_size) |>
            summarize(
                Sensitivity = max(Sensitivity, na.rm = TRUE),
                Specificity = min(Specificity, na.rm = TRUE),
                MCC         = max(MCC, na.rm = TRUE),
                .groups = "drop"
            )
        
        # Anchor at model_size = 0 (0 sensitivity, 1 specificity, 0 MCC)
        if (!0 %in% clean_rep$model_size) {
            clean_rep <- bind_rows(tibble(model_size = 0, Sensitivity = 0, Specificity = 1, MCC = 0), clean_rep)
        }
        
        grid <- tibble(model_size = 0:max_p)
        res <- clean_rep |>
            full_join(grid, by = "model_size") |>
            arrange(model_size)
        
        # Extrapolate using rule = 2
        res <- res |> mutate(
            Sensitivity = zoo::na.approx(Sensitivity, model_size, na.rm = FALSE, rule = 2),
            Specificity = zoo::na.approx(Specificity, model_size, na.rm = FALSE, rule = 2),
            MCC         = zoo::na.approx(MCC, model_size, na.rm = FALSE, rule = 2),
            model = .x$model[1],
            setting = .x$setting[1],
            dim = .x$dim[1],
            rep = .x$rep[1],
            p = .x$p[1],
            k = .x$k[1]
        )
        res
    }) |>
    filter(model_size > 0)

plot_data_tpr <- sum_tpr |>
    group_by(model_size, p, k, setting, model, dim) |>
    summarise(
        Sens_mean = mean(Sensitivity, na.rm = TRUE),
        Sens_sd   = sd(Sensitivity,   na.rm = TRUE),
        n_reps    = sum(!is.na(Sensitivity)),
        .groups = "drop"
    ) |>
    mutate(
        Sens_sd = if_else(is.na(Sens_sd), 0, Sens_sd),
        ymin = pmax(0, Sens_mean - 1.96 * (Sens_sd / sqrt(n_reps))),
        ymax = pmin(1, Sens_mean + 1.96 * (Sens_sd / sqrt(n_reps)))
    ) |>
    arrange(p, model_size)

p1 <- ggplot(plot_data_tpr, aes(x = model_size, y = Sens_mean, color = model)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = model), alpha = 0.15, colour = NA) +
    geom_path(linewidth = 0.7) +
    facet_grid(paste0("Setting ", setting) ~ fct_inorder(dim), scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_color_manual(values = cpl_palette) +
    scale_fill_manual(values = cpl_palette) +
    labs(x = "Model Size", y = "Sensitivity (TPR)", color = "Model", fill = "Model")

ggsave(
    filename = "selection-tpr.png",
    plot = p1,
    path = figs_dir,
    width = 200, height = 120, units = "mm", dpi = 300, bg = "transparent"
)
cli::cli_alert_success("Generated selection-tpr.png")


# ==============================================================================
# Plot 2: Matthews Correlation Coefficient (MCC) vs. Model Size (with 95% CI)
# ==============================================================================
cli::cli_alert_info("Interpolating MCC vs. Model Size...")

plot_data_mcc <- sum_tpr |>
    group_by(model_size, p, k, setting, model, dim) |>
    summarise(
        MCC_mean = mean(MCC, na.rm = TRUE),
        MCC_sd   = sd(MCC,   na.rm = TRUE),
        n_reps   = sum(!is.na(MCC)),
        .groups = "drop"
    ) |>
    mutate(
        MCC_sd = if_else(is.na(MCC_sd), 0, MCC_sd),
        ymin = pmax(-1, MCC_mean - 1.96 * (MCC_sd / sqrt(n_reps))),
        ymax = pmin(1, MCC_mean + 1.96 * (MCC_sd / sqrt(n_reps)))
    ) |>
    arrange(p, model_size)

p2 <- ggplot(plot_data_mcc, aes(x = model_size, y = MCC_mean, color = model)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = model), alpha = 0.15, colour = NA) +
    geom_path(linewidth = 0.7) +
    facet_grid(paste0("Setting ", setting) ~ fct_inorder(dim), scales = "free_x") +
    coord_cartesian(ylim = c(-0.1, 1)) +
    scale_color_manual(values = cpl_palette) +
    scale_fill_manual(values = cpl_palette) +
    labs(x = "Model Size", y = "Matthews Correlation Coefficient (MCC)", color = "Model", fill = "Model")

ggsave(
    filename = "selection-mcc.png",
    plot = p2,
    path = figs_dir,
    width = 200, height = 120, units = "mm", dpi = 300, bg = "transparent"
)
cli::cli_alert_success("Generated selection-mcc.png")


# ==============================================================================
# Plot 3: Variable Selection ROC Curve (Sensitivity vs. 1 - Specificity with 95% CI)
# ==============================================================================
cli::cli_alert_info("Interpolating Variable Selection ROC curves...")

common_fpr <- seq(0, 1, length.out = 101)

sum_roc <- select_mod |>
    group_by(model, setting, dim, rep, p, k) |>
    group_split() |>
    map_dfr(~ {
        # Calculate FPR and TPR
        clean_rep <- .x |>
            mutate(FPR = pmax(0, pmin(1, 1 - Specificity)),
                   TPR = pmax(0, pmin(1, Sensitivity))) |>
            select(FPR, TPR) |>
            filter(!is.na(FPR) & !is.na(TPR))
        
        # Add anchors (0,0) and (1,1)
        clean_rep <- bind_rows(
            tibble(FPR = 0, TPR = 0),
            clean_rep,
            tibble(FPR = 1, TPR = 1)
        ) |>
            arrange(FPR, TPR) |>
            group_by(FPR) |>
            summarize(TPR = max(TPR), .groups = "drop")
        
        # Interpolate onto common FPR grid
        interp_tpr <- approx(x = clean_rep$FPR, y = clean_rep$TPR, xout = common_fpr, rule = 2)$y
        
        tibble(
            FPR = common_fpr,
            TPR = interp_tpr,
            model = .x$model[1],
            setting = .x$setting[1],
            dim = .x$dim[1],
            rep = .x$rep[1],
            p = .x$p[1],
            k = .x$k[1]
        )
    })

plot_data_roc <- sum_roc |>
    group_by(FPR, p, k, setting, model, dim) |>
    summarise(
        TPR_mean = mean(TPR, na.rm = TRUE),
        TPR_sd   = sd(TPR,   na.rm = TRUE),
        n_reps   = sum(!is.na(TPR)),
        .groups = "drop"
    ) |>
    mutate(
        TPR_sd = if_else(is.na(TPR_sd), 0, TPR_sd),
        ymin = pmax(0, TPR_mean - 1.96 * (TPR_sd / sqrt(n_reps))),
        ymax = pmin(1, TPR_mean + 1.96 * (TPR_sd / sqrt(n_reps)))
    )

p3 <- ggplot(plot_data_roc, aes(x = FPR, y = TPR_mean, color = model)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = model), alpha = 0.15, colour = NA) +
    geom_path(linewidth = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    facet_grid(paste0("Setting ", setting) ~ fct_inorder(dim)) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_x_continuous(breaks = c(0, .25, .5, .75, 1)) +
    scale_y_continuous(breaks = c(0, .25, .5, .75, 1)) +
    scale_color_manual(values = cpl_palette) +
    scale_fill_manual(values = cpl_palette) +
    labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", color = "Model", fill = "Model")

ggsave(
    filename = "selection-curve.png",
    plot = p3,
    path = figs_dir,
    width = 180, height = 180, units = "mm", dpi = 300, bg = "transparent"
)
cli::cli_alert_success("Generated selection-curve.png")

cli::cli_alert_info("Successfully generated selection-tpr.png, selection-mcc.png, and selection-curve.png in figs/.")
