# ==============================================================================
# cbSCRIP Simulation Pipeline (01_simulations.R)
# 
# Replications: Managed by bash scripts (e.g., run-exp_p-120.sh)
# Test size: 10,000
# Models: cbSCRIP, enet-iCox, Penalized Fine-Gray (fastcmprsk), rfsrc, Aalen-Johansen
# ==============================================================================

# --- Configuration ---
n_train <- 300
n_test <- 10000

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Requires 4 arguments: p, num_true, setting, iter")
}
p <- as.numeric(args[1])
num_true <- as.integer(args[2])
setting <- as.numeric(args[3])
iter <- as.numeric(args[4])
overwrite <- as.logical(Sys.getenv("OVERWRITE", unset = "TRUE"))
if (is.na(overwrite)) overwrite <- TRUE
rm(args)

save_dir_varsel <- here::here("results", "simulation", "varsel")
save_dir_preds <- here::here("results", "simulation", "preds")
dir.create(save_dir_varsel, recursive = TRUE, showWarnings = FALSE)
dir.create(save_dir_preds, recursive = TRUE, showWarnings = FALSE)

save_file_varsel <- file.path(save_dir_varsel, glue::glue("varsel_s{setting}_p{p}_k{num_true}_rep{iter}.rds"))
save_file_preds <- file.path(save_dir_preds, glue::glue("models_s{setting}_p{p}_k{num_true}_rep{iter}.rds"))

if (!overwrite && file.exists(save_file_varsel) && file.exists(save_file_preds)) {
    cli::cli_alert_success(glue::glue("Skipping setting {setting}, p={p}, rep={iter} as files exist."))
    quit(save = "no")
}


suppressPackageStartupMessages({
  library(casebase)
  library(glmnet)
  library(parallel)
  library(tictoc)
  library(dplyr)
  library(tidyr)
  library(survival)
  library(cmprsk)
  library(here)
  library(purrr)
  library(ggplot2)
  library(glue)
  if (dir.exists(here::here("cbSCRIP"))) {
    pkgload::load_all(here::here("cbSCRIP"))
  } else {
    library(cbSCRIP)
  }
  library(fastcmprsk)
  library(randomForestSRC)
  library(riskRegression)
  library(pec)
  library(CoxBoost)
})



cli::cli_alert_info(paste0(glue("Starting setting {setting}, p={p}, rep={iter}")))

# Generate data
data <- gen_data(n_train = n_train, n_test = n_test, p = p, 
                 num_true = num_true, setting = setting, iter = iter)
train <- data$train
test <- data$test
beta1 <- data$beta1
beta2 <- data$beta2

# Outcomes for standard survival packages
y_train_cox <- Surv(time = train$ftime, event = train$fstatus == 1)
x_train <- model.matrix(~ . - ftime - fstatus, data = train)[, -1]

y_test_cox <- Surv(time = test$ftime, event = test$fstatus == 1)
x_test <- model.matrix(~ . - ftime - fstatus, data = test)[, -1]

lmr <- ifelse(nrow(x_train) < ncol(x_train), 0, 0.001)

fit_all_models <- function(train, x_train, y_train_cox, x_test, y_test_cox) {
    # 1. enet-iCox (two.cv.CSlassos for prediction/model comparison)
    cox_mod_cv <- NULL
    t_cox <- NA
    tryCatch({
        t_cox <- system.time({
            cox_mod_cv <- cbSCRIP::two.cv.CSlassos(data = train, nfold = 5,
             lambda.select = "lambda.1se", nlambda = 30, lambda.min.ratio = lmr)
        })["elapsed"]
    }, error = function(e) {
        cli::cli_alert_info(paste0("   two.cv.CSlassos failed: ", e$message))
    })

    # Direct glmnet path model for variable selection
    y_train_cox1 <- survival::Surv(time = train$ftime, event = train$fstatus == 1)
    lm_ratio <- ifelse(length(y_train_cox) < ncol(x_train), 0, 0.001)
    cox_mod_path <- NULL
    tryCatch({
        cox_mod_path <- glmnet::glmnet(
            x = x_train,
            y = y_train_cox1,
            family = "cox",
            alpha = 0.5,
            lambda.min.ratio = lm_ratio,
            nlambda = 30
        )
    }, error = function(e) {
        cli::cli_alert_info(paste0("   glmnet path fitting failed: ", e$message))
    })

    # 2. cbSCRIP (Elastic Net Case-Base)
    cb_mod <- NULL
    t_cb <- NA
    tryCatch({
        t_cb <- system.time({
            cb_mod <- cbSCRIP(Surv(ftime, fstatus) ~ ., train, coeffs = "original", nlambda = 30, lambda.min.ratio = lmr)
        })["elapsed"]
    }, error = function(e) {
        cli::cli_alert_info(paste0("   cbSCRIP fitting failed: ", e$message))
    })

    # 3. Penalized Fine-Gray (fastcmprsk)
    fg_mod <- NULL
    t_fg <- NA
    tryCatch({
        t_fg <- system.time({
            lambdas <- if (!is.null(cox_mod_path)) cox_mod_path$lambda else seq(0.1, 0.0001, length.out = 50)
            fg_mod <- fastCrrp(Crisk(train$ftime, train$fstatus) ~ x_train,
                               penalty = "ENET", alpha = 0.5, lambda = lambdas)
        })["elapsed"]
    }, error = function(e) {
        cli::cli_alert_info(paste0("   fastCrrp fitting failed: ", e$message))
    })

    # 4. Random Forest SRC
    rf_mod <- NULL
    t_rf <- NA
    tryCatch({
        t_rf <- system.time({
            rf_mod <- rfsrc(Surv(ftime, fstatus) ~ ., data = train, 
                            ntree = 500, splitrule = "logrank", cause = 1)
        })["elapsed"]
    }, error = function(e) {
        cli::cli_alert_info(paste0("   rfsrc fitting failed: ", e$message))
    })

    # 5. Aalen-Johansen (for reference prediction)
    aj_mod <- NULL
    tryCatch({
        aj_mod <- prodlim(Hist(ftime, fstatus) ~ 1, data = train)
    }, error = function(e) {
        cli::cli_alert_info(paste0("   prodlim fitting failed: ", e$message))
    })
    
    # 6. SHBoost (iCoxBoost)
    sh_mod <- NULL
    t_sh <- NA
    tryCatch({
        t_sh <- system.time({
            optim_res <- tryCatch({
                iCoxBoost(Hist(ftime, fstatus) ~ ., data = train, cause = 1, cv = TRUE, maxstepno = 60)
            }, error = function(e) NULL)
            
            opt_step <- if (!is.null(optim_res)) optim_res$cv.res$optimal.step else 30
            
            sh_mod <- iCoxBoost(Hist(ftime, fstatus) ~ ., data = train, cause = 1, stepno = opt_step)
        })["elapsed"]
    }, error = function(e) {
        cli::cli_alert_info(paste0("   iCoxBoost fitting failed: ", e$message))
    })
    
    # Clean formula environments to prevent circular references
    if (!is.null(cb_mod) && !is.null(cb_mod$formula)) attr(cb_mod$formula, ".Environment") <- globalenv()
    if (!is.null(aj_mod) && !is.null(aj_mod$formula)) attr(aj_mod$formula, ".Environment") <- globalenv()
    if (!is.null(rf_mod) && !is.null(rf_mod$formula)) attr(rf_mod$formula, ".Environment") <- globalenv()
    if (!is.null(sh_mod) && !is.null(sh_mod$formula)) attr(sh_mod$formula, ".Environment") <- globalenv()
    
    return(list(
        cox_mod_path = cox_mod_path,
        cox_mod_cv = cox_mod_cv,
        cb_mod = cb_mod,
        fg_mod = fg_mod,
        rf_mod = rf_mod,
        aj_mod = aj_mod,
        sh_mod = sh_mod,
        t_cox = t_cox,
        t_cb = t_cb,
        t_fg = t_fg,
        t_rf = t_rf,
        t_sh = t_sh
    ))
}

# Run the clean fit
fits <- fit_all_models(train, x_train, y_train_cox, x_test, y_test_cox)

cox_mod_path <- fits$cox_mod_path
cox_mod_cv <- fits$cox_mod_cv
cb_mod <- fits$cb_mod
fg_mod <- fits$fg_mod
rf_mod <- fits$rf_mod
aj_mod <- fits$aj_mod
sh_mod <- fits$sh_mod
t_cox <- fits$t_cox
t_cb <- fits$t_cb
t_fg <- fits$t_fg
t_rf <- fits$t_rf
t_sh <- fits$t_sh

cli::cli_alert_info(paste0(sprintf("   - Fit times: enet-iCox=%s | cbSCRIP=%s | Penalized Fine-Gray=%s | Random Forest=%s | SHBoost=%s", 
                if (is.na(t_cox)) "FAILED" else sprintf("%.2fs", t_cox),
                if (is.na(t_cb)) "FAILED" else sprintf("%.2fs", t_cb),
                if (is.na(t_fg)) "FAILED" else sprintf("%.2fs", t_fg),
                if (is.na(t_rf)) "FAILED" else sprintf("%.2fs", t_rf),
                if (is.na(t_sh)) "FAILED" else sprintf("%.2fs", t_sh))))



# --- Variable Selection Metrics (for penalized models) ---

# enet-iCox path
varsel_cox <- NULL
if (!is.null(cox_mod_path)) {
    varsel_cox <- tryCatch({
        map_df(cox_mod_path$lambda, function(x) {
            coefs <- as.numeric(coef(cox_mod_path, s = x)[, 1])
            res <- as_tibble(varsel_perc(coefs, beta1))
            res %>% mutate(lambda = x, model = "enet-iCox", 
                           model_size = TP + FP, 
                           exact_match = (TP == num_true && FP == 0))
        })
    }, error = function(e) NULL)
}

# cbSCRIP path
varsel_cb <- NULL
if (!is.null(cb_mod)) {
    varsel_cb <- tryCatch({
        map_df(seq_along(cb_mod$lambdagrid), function(x) {
            coefs <- cb_mod$coefficients[[x]][seq_along(beta1), 1]
            res <- as_tibble(varsel_perc(coefs, beta1))
            res %>% mutate(lambda = cb_mod$lambdagrid[x], model = "cbSCRIP", 
                           model_size = TP + FP, 
                           exact_match = (TP == num_true && FP == 0))
        })
    }, error = function(e) NULL)
}

# Penalized Fine-Gray path
varsel_fg <- NULL
if (!is.null(fg_mod)) {
    varsel_fg <- tryCatch({
        purrr::map_df(seq_along(fg_mod$lambda.path), function(x) {
            coefs <- as.numeric(fg_mod$coef[, x])
            res <- as_tibble(varsel_perc(coefs, beta1))
            res %>% mutate(lambda = fg_mod$lambda.path[x], model = "Penalized Fine-Gray", 
                           model_size = TP + FP, 
                           exact_match = (TP == num_true && FP == 0))
        })
    }, error = function(e) {
        cli::cli_alert_info(paste0("   varsel_fg failed: ", e$message))
        NULL
    })
}

varsel_results <- bind_rows(varsel_cox, varsel_cb, varsel_fg)
if (nrow(varsel_results) > 0) {
    varsel_results <- varsel_results %>% 
        mutate(rep = iter, setting = setting, p = p, k = num_true)
}

# --- Variable Selection ROC Plot (untracked output) ---
figs_dir <- here("figs", "varsel_plots")
dir.create(figs_dir, showWarnings = FALSE, recursive = TRUE)

plot_varsel <- varsel_results |>
    ggplot(aes(
        x = 1 - Specificity,
        y = Sensitivity,
        color = model,
        linetype = model
    )) +
    geom_path() +
    coord_equal() +
    labs(
        title = glue("Variable Selection ROC"),
        subtitle = glue("Setting {setting}, p={p}, k={num_true}, rep={iter}"),
        x = "1 - Specificity (FPR)",
        y = "Sensitivity (TPR)",
        color = "Model", linetype = "Model"
    ) +
    theme_bw()

ggsave(
    file.path(figs_dir,
              glue("varsel_roc_s{setting}_p{p}_k{num_true}_rep{iter}.png")),
    plot = plot_varsel,
    bg = "transparent",
    width = 200,
    height = 120,
    units = "mm",
    dpi = 300
)


# Performance updates for console output
max_tpr_cb <- if (!is.null(varsel_cb)) max(varsel_cb$Sensitivity, na.rm = TRUE) else NA
max_tpr_cox <- if (!is.null(varsel_cox)) max(varsel_cox$Sensitivity, na.rm = TRUE) else NA
max_tpr_fg <- if (!is.null(varsel_fg)) max(varsel_fg$Sensitivity, na.rm = TRUE) else NA

cli::cli_alert_info(paste0(sprintf("   - Max Sensitivity (TPR): cbSCRIP=%s | enet-iCox=%s | Penalized Fine-Gray=%s", 
                if (is.na(max_tpr_cb)) "FAILED" else sprintf("%.1f%%", max_tpr_cb * 100), 
                if (is.na(max_tpr_cox)) "FAILED" else sprintf("%.1f%%", max_tpr_cox * 100), 
                if (is.na(max_tpr_fg)) "FAILED" else sprintf("%.1f%%", max_tpr_fg * 100))))
cli::cli_alert_info(paste0(sprintf("   - Exact target selection: cbSCRIP=%s | enet-iCox=%s | Penalized Fine-Gray=%s", 
                if (!is.null(varsel_cb) && any(varsel_cb$exact_match)) "YES" else "NO",
                if (!is.null(varsel_cox) && any(varsel_cox$exact_match)) "YES" else "NO",
                if (!is.null(varsel_fg) && any(varsel_fg$exact_match)) "YES" else "NO")))

# Save Variable Selection Output
saveRDS(varsel_results, file.path(save_dir_varsel, glue("varsel_s{setting}_p{p}_k{num_true}_rep{iter}.rds")))

# # --- Evaluate Predictive Metrics (Brier Score, AUC) ---
cli::cli_alert_info(paste0("   Computing Brier scores and AUC..."))

# Time grid: 20 quantile points within observed event range
brier_tp <- sort(unique(test$ftime))
brier_tp <- brier_tp[is.finite(brier_tp) & brier_tp > 0]
brier_tmax <- max(test$ftime[test$fstatus != 0], na.rm = TRUE)
brier_tp <- brier_tp[brier_tp <= brier_tmax]
brier_tp <- sort(unique(quantile(brier_tp, probs = seq(0, 1, length.out = 20), type = 1)))

# --- CV-based unpenalized refitting (original paper approach) ---
cli::cli_alert_info(paste0("   Running cross-validation to select optimal lambda..."))



fit_cv <- tryCatch(
    cv_cbSCRIP(Surv(ftime, fstatus) ~ ., cb_data = cb_mod$cb_data, nlambda = 30),
    error = function(e) {
        cli::cli_alert_info(paste0("   cv_cbSCRIP failed: ", e$message))
        NULL
    }
)

cb_mod_refitted <- NULL
if (!is.null(fit_cv)) {
    # Save the CV plot (since plot.cbSCRIP.cv returns a ggplot object, we use ggsave)
    cv_plots_dir <- here("figs", "cv_plots")
    dir.create(cv_plots_dir, showWarnings = FALSE, recursive = TRUE)
    
    cv_plot <- tryCatch(plot(fit_cv), error = function(e) NULL)
    if (!is.null(cv_plot)) {
        ggsave(file.path(cv_plots_dir, glue("cv_s{setting}_p{p}_k{num_true}_rep{iter}.png")),
               plot = cv_plot, width = 200, height = 150, units = "mm", dpi = 300)
    }
    
    # Fit adjusted (unpenalized) model at lambda.1se
    cli::cli_alert_info(paste0(glue("   Fitting adjusted model at lambda.1se: {fit_cv$lambda.1se}")))
    cb_mod_refitted <- tryCatch(
        cbSCRIP(cb_data = fit_cv$cb_data,
                lambda = fit_cv$lambda.1se),
        error = function(e) {
            cli::cli_alert_info(paste0("   Adjusted model refit failed at lambda.1se: ", e$message))
            NULL
        }
    )
}

# enet-iCox (cox_mod_cv) is already cross-validated via two.cv.CSlassos

# Put together models list for riskRegression::Score()
model_list_pred <- Filter(Negate(is.null), list(
    "Aalen-Johansen" = aj_mod,
    "Random Forest"  = rf_mod,
    "enet-iCox"      = cox_mod_cv,
    "SHBoost"        = sh_mod
))
if (!is.null(cb_mod_refitted))
    model_list_pred[["cbSCRIP"]] <- cb_mod_refitted$refitted_models

# Run Score for Brier and AUC
sc <- tryCatch(
    Score(model_list_pred,
          formula = Hist(ftime, fstatus) ~ 1,
          data    = test, times = brier_tp,
          summary = "ibs", se.fit = FALSE,
          metrics = c("Brier", "AUC"), cause = 1),
    error = function(e) {
        cli::cli_alert_info(paste0("   Score() failed: ", e$message))
        NULL
    }
)

brier_table <- NULL
auc_table <- NULL

if (!is.null(sc)) {
    brier_table <- tryCatch({
        sc$Brier$score |>
            filter(model != "Null model") |>
            mutate(model = as.character(model),
                   rep = iter, setting = setting, p = p, k = num_true)
    }, error = function(e) NULL)
    
    auc_table <- tryCatch({
        sc$AUC$score |>
            mutate(model = as.character(model),
                   rep = iter, setting = setting, p = p, k = num_true)
    }, error = function(e) NULL)
}

if (!is.null(brier_table))
    cli::cli_alert_info(paste0(sprintf("   Brier: %d rows, models: %s",
                    nrow(brier_table), paste(unique(brier_table$model), collapse = ", "))))
if (!is.null(auc_table))
    cli::cli_alert_info(paste0(sprintf("   AUC: %d rows, models: %s",
                    nrow(auc_table), paste(unique(auc_table$model), collapse = ", "))))


# --- CIF Curves ---
cli::cli_alert_info(paste0("   Computing CIF curves..."))

cif_table <- tryCatch({
    rows <- list()

    # 1. cbSCRIP: absoluteRisk via refitted CompRisk object
    if (!is.null(cb_mod_refitted)) {
        r_cb <- predictRisk(
            object = cb_mod_refitted$refitted_models,
            newdata = test,
            times = brier_tp,
            cause = 1
        )
        risk_cb <- as.numeric(colMeans(r_cb, na.rm = TRUE))
        n_tp <- length(brier_tp)
        if (length(risk_cb) == n_tp) {
            rows[["cbSCRIP"]] <- tibble(Time = brier_tp, Risk = risk_cb,
                                        Method = "cbSCRIP", rep = iter,
                                        setting = setting, p = p, k = num_true)
        }
    }

    # 2. Aalen-Johansen: population-level curve
    r_aj <- predictRisk(aj_mod, newdata = test[1, , drop = FALSE],
                        times = brier_tp, cause = 1)
    risk_aj <- as.numeric(r_aj)
    n_tp <- length(brier_tp)
    if (length(risk_aj) == n_tp) {
        rows[["Aalen-Johansen"]] <- tibble(Time = brier_tp, Risk = risk_aj,
                                           Method = "Aalen-Johansen", rep = iter,
                                           setting = setting, p = p, k = num_true)
    }

    # 3. enet-iCox (twoCSlassos)
    r_cox <- tryCatch({
        predictRisk(cox_mod_cv, newdata = test, times = brier_tp, cause = 1)
    }, error = function(e) NULL)
    if (!is.null(r_cox)) {
        risk_cox <- as.numeric(colMeans(r_cox, na.rm = TRUE))
        if (length(risk_cox) == n_tp) {
            rows[["enet-iCox"]] <- tibble(Time = brier_tp, Risk = risk_cox,
                                          Method = "enet-iCox", rep = iter,
                                          setting = setting, p = p, k = num_true)
        }
    }

    # 4. SHBoost (iCoxBoost)
    r_sh <- tryCatch({
        predictRisk(sh_mod, newdata = test, times = brier_tp, cause = 1)
    }, error = function(e) NULL)
    if (!is.null(r_sh)) {
        risk_sh <- as.numeric(colMeans(r_sh, na.rm = TRUE))
        if (length(risk_sh) == n_tp) {
            rows[["SHBoost"]] <- tibble(Time = brier_tp, Risk = risk_sh,
                                         Method = "SHBoost", rep = iter,
                                         setting = setting, p = p, k = num_true)
        }
    }

    bind_rows(rows)
}, error = function(e) {
    cli::cli_alert_info(paste0("   CIF failed: ", e$message))
    NULL
})

if (!is.null(cif_table))
    cli::cli_alert_info(paste0(sprintf("   CIF: %d rows, methods: %s",
                    nrow(cif_table), paste(unique(cif_table$Method), collapse = ", "))))


# --- Save results bundle ---
results_bundle <- list(
    meta = list(N = n_train, p = p, k = num_true,
                setting = setting, seed = iter),
    data = list(test = test),
    models = list(
        cox_mod = cox_mod_cv,
        cox_mod_path = cox_mod_path,
        cb_mod = cb_mod,
        fg_mod = fg_mod,
        rf_mod = rf_mod,
        aj_mod = aj_mod,
        sh_mod = sh_mod
    ),
    timing = list(
        enet_iCox = t_cox,
        cbSCRIP = t_cb,
        Penalized_Fine_Gray = t_fg,
        Random_Forest = t_rf,
        SHBoost = t_sh
    ),
    brier_table = brier_table,
    auc_table = auc_table,
    cif_table = cif_table
)

saveRDS(results_bundle, file.path(save_dir_preds, glue("models_s{setting}_p{p}_k{num_true}_rep{iter}.rds")))
cli::cli_alert_info(paste0(glue("Completed setting {setting}, p={p}, k={num_true}, rep={iter}")))

# --- Brier Score Plot (untracked output) ---
if (!is.null(brier_table) && nrow(brier_table) > 0) {
    brier_plots_dir <- here("figs", "brier_plots")
    dir.create(brier_plots_dir, showWarnings = FALSE, recursive = TRUE)
    
    plot_brier <- brier_table |>
        ggplot(aes(x = times, y = Brier, color = model, linetype = model)) +
        geom_line() +
        labs(
            title = glue("Brier Score per Run"),
            subtitle = glue("Setting {setting}, p={p}, k={num_true}, rep={iter}"),
            x = "Follow-up time (years)",
            y = "Brier Score",
            color = "Model", linetype = "Model"
        ) +
        theme_bw()
        
    ggsave(
        file.path(brier_plots_dir,
                  glue("brier_s{setting}_p{p}_k{num_true}_rep{iter}.png")),
        plot = plot_brier,
        bg = "transparent",
        width = 200,
        height = 120,
        units = "mm",
        dpi = 300
    )
}

# --- CIF Plot (untracked output) ---
if (!is.null(cif_table) && nrow(cif_table) > 0) {
    cif_plots_dir <- here("figs", "cif_plots")
    dir.create(cif_plots_dir, showWarnings = FALSE, recursive = TRUE)
    
    plot_cif <- cif_table |>
        ggplot(aes(x = Time, y = Risk, color = Method, linetype = Method)) +
        geom_line() +
        labs(
            title = glue("Cumulative Incidence Function (Cause 1)"),
            subtitle = glue("Setting {setting}, p={p}, k={num_true}, rep={iter}"),
            x = "Follow-up time (years)",
            y = "Absolute Risk",
            color = "Method", linetype = "Method"
        ) +
        theme_bw()
        
    ggsave(
        file.path(cif_plots_dir,
                  glue("cif_s{setting}_p{p}_k{num_true}_rep{iter}.png")),
        plot = plot_cif,
        bg = "transparent",
        width = 200,
        height = 120,
        units = "mm",
        dpi = 300
    )
}

# --- AUC Plot (untracked output) ---
if (!is.null(auc_table) && nrow(auc_table) > 0) {
    auc_plots_dir <- here("figs", "auc_plots")
    dir.create(auc_plots_dir, showWarnings = FALSE, recursive = TRUE)
    
    plot_auc <- auc_table |>
        ggplot(aes(x = times, y = AUC, color = model, linetype = model)) +
        geom_line() +
        labs(
            title = glue("Time-dependent AUC per Run"),
            subtitle = glue("Setting {setting}, p={p}, k={num_true}, rep={iter}"),
            x = "Follow-up time (years)",
            y = "AUC",
            color = "Model", linetype = "Model"
        ) +
        theme_bw()
        
    ggsave(
        file.path(auc_plots_dir,
                  glue("auc_s{setting}_p{p}_k{num_true}_rep{iter}.png")),
        plot = plot_auc,
        bg = "transparent",
        width = 200,
        height = 120,
        units = "mm",
        dpi = 300
    )
}


