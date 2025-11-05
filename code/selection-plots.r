#################################################################
##                       Selection plots                       ##
#################################################################
##
## Description:    
##                 
##
## Author:         Javier Mtz.-Rdz.  
##
## Creation date:  2025-09-23
##
## Email:          javier.mr@stat.ubc.ca
##
## ---------------------------
## Notes:          
## ---------------------------

# Setup ----
## Packages to use ----
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
    tidyverse, here, glue,
    readxl, writexl, janitor, kableExtra,
    survival, cmprsk, pec, casebase, survminer,
    glmnet, patchwork, scales, zoo,
    future.apply, foreach, parallel, tictoc
)

## Set theme ------
# mytidyfunctions::set_mytheme(text = element_text(family = "Times New Roman"))

# Source local helper functions.
source(here("notes_jmr", "code", "fitting_functionsV2.R"))

# Load data

# Define file paths.
results_dir <- here("paper", "results")
figs_dir <- here("paper", "figs")
select_mod_path <- file.path(results_dir, "select_mod.rds")

# Create directories if they don't exist.
dir.create(figs_dir, showWarnings = FALSE, recursive = TRUE)

# Read processed data if it exists, otherwise generate it.
if (file.exists(select_mod_path)) {
    select_mod <- readRDS(select_mod_path)
} else {
    # Find all raw model result files.
    model_files <- list.files(
        path = file.path(results_dir, "lambda_paths"),
        pattern = "models_",
        recursive = TRUE,
        full.names = TRUE
    )
    
    # Tidy file metadata.
    data_proj <- tibble(file_path = model_files) %>%
        mutate(
            name = str_extract(file_path, "(?<=models_)(.*?)(?=.rds)"),
            name = str_replace(name, "selec", "setting")
        ) %>%
        separate(
            name,
            into = c("setting", "sim", "p", "k"),
            sep = "_",
            remove = FALSE,
            convert = TRUE # Automatically converts to numeric
        ) %>% 
        mutate(
            across(c(setting, sim, p, k), ~ as.numeric(str_remove(., "^.*-"))),
            dim = glue("p = {p}, k = {k}")
        ) %>%
        arrange(p, k)
    
    # Read and combine individual model coefficients.
    
    select_mod <- map_df(1:nrow(data_proj),
                         function(i){data_proj[i,-1] %>%
                                 select(-sim) %>%
                                 bind_cols(suppressMessages(readRDS(data_proj$file_path[i])) %>%
                                               mutate(model_size = TP + FP,
                                                      model = if_else(model == "enet-casebase-Acc", "cbSCRIP", model))
                                 )})
    
    
    # Save the processed data frame.
    saveRDS(select_mod, select_mod_path)
}

select_mod <- select_mod %>%
    filter(!str_detect(model, "SCAD")) |> 
    mutate(model = case_match(model, "enet-CR" ~ "enet-iCox", 
                               .default = model))


# Interpolate metrics for missing model sizes.
sum_selec <- select_mod %>%
    group_by(model, setting, dim, sim, p, k) %>%
    group_split() %>%
    map_dfr(~ {
        .x %>%
            full_join(tibble(model_size = 1:.x$p[1]), by = "model_size") %>%
            arrange(model_size) %>%
            mutate(across(c(Sensitivity, Specificity, MCC),
                          ~ zoo::na.approx(.x, model_size, na.rm = FALSE)
            )) %>%
            fill(model, setting, dim, sim, p, k, .direction = "downup")
    })

# Summarize metrics across simulations.
plot_data <- sum_selec %>%
    group_by(model_size, p, k, setting, model, dim) %>%
    summarise(
        across(c(Sensitivity, Specificity, MCC), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
    ) %>%
    arrange(p, model_size)


# Plot 1: Sensitivity vs. Model Size
(p1 <- ggplot(plot_data, aes(x = model_size, y = Sensitivity, color = model)) +
    geom_path() +
    facet_grid(paste0("Setting ", setting) ~ fct_inorder(dim), scales = "free_x") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Model Size", color = "Model"))

ggsave(
    filename = "selection-tpr.png",
    plot = p1,
    path = figs_dir,
    width = 200, height = 120, units = "mm", dpi = 300, bg = "transparent"
)

# Plot 2: MCC vs. Model Size
(p2 <- ggplot(plot_data, aes(x = model_size, y = MCC, color = model)) +
    geom_line() +
    facet_grid(paste0("Setting ", setting) ~ fct_inorder(dim), scales = "free_x") +
    coord_cartesian(ylim = c(-0.1, 1)) +
    labs(x = "Model Size", color = "Model"))

ggsave(
    filename = "selection-mcc.png",
    plot = p2,
    path = figs_dir,
    width = 200, height = 120, units = "mm", dpi = 300, bg = "transparent"
)

# Plot 3: ROC-like Curve (Sensitivity vs. 1 - Specificity)
(p3 <- ggplot(plot_data, aes(x = 1 - Specificity, y = Sensitivity, color = model)) +
    geom_path() +
    facet_grid(paste0("Setting ", setting) ~ fct_inorder(dim)) +
    coord_equal() +
    scale_x_continuous(breaks = c(.25, .5, .75, 1)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "1 - Specificity", color = "Model"))

ggsave(
    filename = "selection-curve.png",
    plot = p3,
    path = figs_dir,
    width = 180, height = 180, units = "mm", dpi = 300, bg = "transparent"
)

