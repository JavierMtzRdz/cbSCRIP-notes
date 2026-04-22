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
    res <- readRDS(f)
    # Ensure columns match expected types
    res |> mutate(
        model = as.character(model),
        rep = as.numeric(rep),
        setting = as.numeric(setting),
        p = as.numeric(p)
    )
})

# 2. Add k and dim
select_mod <- select_mod |>
    mutate(
        k = case_match(p,
            120 ~ 20,
            500 ~ 84,
            1000 ~ 168,
            .default = 20
        ),
        dim = glue("p = {p}, k = {k}")
    )

# 3. Clean model names
select_mod <- select_mod |>
    mutate(model = case_match(model,
        "enet-CR" ~ "enet-iCox",
        "fastcmprsk" ~ "Penalized Fine-Gray",
        "penalized Fine-Gray method (fastcmprsk)" ~ "Penalized Fine-Gray",
        .default = model
    ))

# 4. Interpolate metrics for missing model sizes
cli::cli_alert_info(paste0("Interpolating metrics..."))
sum_selec <- select_mod |>
    group_by(model, setting, dim, rep, p, k) |>
    group_split() |>
    map_dfr(~ {
        # Determine unique model size range up to p
        max_p <- .x$p[1]
        grid <- tibble(model_size = 1:max_p)
        res <- .x |>
            full_join(grid, by = "model_size") |>
            arrange(model_size) |>
            fill(model, setting, dim, rep, p, k, .direction = "downup")
        
        # Only interpolate if we have at least 2 non-NA points
        tryCatch({
            if (sum(!is.na(res$Sensitivity)) >= 2 && sum(!is.na(res$Specificity)) >= 2) {
                res <- res |> mutate(
                    Sensitivity = zoo::na.approx(Sensitivity, model_size, na.rm = FALSE, rule = 2),
                    Specificity = zoo::na.approx(Specificity, model_size, na.rm = FALSE, rule = 2)
                )
            }
        }, error = function(e) { })
        res
    })

# 5. Summarize metrics across simulations
plot_data <- sum_selec |>
    group_by(model_size, p, k, setting, model, dim) |>
    summarise(
        Sens_mean = mean(Sensitivity, na.rm = TRUE),
        Sens_sd   = sd(Sensitivity,   na.rm = TRUE),
        Spec_mean = mean(Specificity, na.rm = TRUE),
        n_reps    = sum(!is.na(Sensitivity)),
        .groups = "drop"
    ) |>
    arrange(p, model_size)

# Palette matching cpl_palette in the supplement
cpl_palette <- c(
    "cbSCRIP"             = "#277DA1",
    "Aalen-Johansen"      = "#43AA8B",
    "Penalized Fine-Gray" = "#F9C74F",
    "enet-iCox"           = "#f94144",
    "Random Forest"       = "#90BE6D",
    "SHBoost"             = "#F9844A"
)

# Set global Times New Roman theme if mytidyfunctions package is available
if (requireNamespace("mytidyfunctions", quietly = TRUE)) {
    mytidyfunctions::set_mytheme(text = element_text(family = "Times New Roman"))
}

# Plot 1: Sensitivity vs. Model Size
p1 <- ggplot(plot_data, aes(x = model_size, y = Sens_mean, color = model)) +
    geom_ribbon(aes(ymin = Sens_mean - 1.96 * (Sens_sd / sqrt(n_reps)),
                    ymax = Sens_mean + 1.96 * (Sens_sd / sqrt(n_reps)),
                    fill = model), alpha = 0.15, colour = NA) +
    geom_path() +
    facet_grid(paste0("Setting ", setting) ~ fct_inorder(dim), scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    scale_color_manual(values = cpl_palette) +
    scale_fill_manual(values = cpl_palette) +
    labs(x = "Model Size", color = "Model", fill = "Model")

ggsave(
    filename = "selection-tpr.png",
    plot = p1,
    path = figs_dir,
    width = 200, height = 120, units = "mm", dpi = 300, bg = "transparent"
)

# Plot 3: ROC-like Curve (Sensitivity vs. 1 - Specificity)
p3 <- ggplot(plot_data, aes(x = 1 - Spec_mean, y = Sens_mean, color = model)) +
    geom_path() +
    facet_grid(paste0("Setting ", setting) ~ fct_inorder(dim)) +
    coord_equal() +
    scale_x_continuous(breaks = c(.25, .5, .75, 1)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    scale_color_manual(values = cpl_palette) +
    labs(x = "1 - Specificity", color = "Model")

ggsave(
    filename = "selection-curve.png",
    plot = p3,
    path = figs_dir,
    width = 180, height = 180, units = "mm", dpi = 300, bg = "transparent"
)

cli::cli_alert_info(paste0("Successfully generated selection-tpr.png and selection-curve.png in figs/."))
