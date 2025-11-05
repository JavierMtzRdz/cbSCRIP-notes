##################################################################
##             Project: Selection Curves Simulation             ##
##################################################################
##
##
## Author:         Javier Mtz.-Rdz.  
##
## Creation date:  2025-07-10
##
## Email:          javier.mr@stat.ubc.ca
##
## ---------------------------
## Notes:          
## ---------------------------

# Setup ----
## Packages to use ----
# Load packages
pacman::p_load(
    casebase,
    glmnet,
    parallel,
    tictoc,
    dplyr,
    tidyr,
    foreach,
    survival,
    cmprsk,
    here,
    caret,
    pec,
    purrr,
    glue,
    cbSCRIP
)

# Simulation-specific parameters
n_lambda <- 50
n <- 400 # Number of samples

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
p <- as.numeric(args[1])
num_true <- as.integer(args[2])
setting <- as.numeric(args[3])
i <- as.numeric(args[4])
rm(args)

# Generate data
data <- gen_data(
    n = n,
    p = p,
    num_true = num_true,
    setting = setting,
    iter = i
)

# Unpack data
beta1 <- data$beta1
beta2 <- data$beta2
train <- data$train
test <- data$test


## 3.  Cause-Specific Cox model ----

# Prepare data for cause-1 (censor competing event)
y_train <- Surv(time = train$ftime, event = train$fstatus == 1)
x_train <- model.matrix(~ . - ftime - fstatus, data = train)[, -1]

y_test <- Surv(time = test$ftime, event = test$fstatus == 1)
x_test <- model.matrix(~ . - ftime - fstatus, data = test)[, -1]

# Fit cause-specific cox model
tic("cox_mod")
cox_mod <- glmnet(
    x = x_train,
    y = y_train,
    family = "cox",
    alpha = 0.5,
    lambda.min.ratio = ifelse(length(y_test) < ncol(x_test), 0, 0.001),
    nlambda = 50
)
toc()

# Process Cox model results
table_result_cox <- map_df(cox_mod$lambda, \(x) {
    as_tibble(varsel_perc(coef(cox_mod, s = x), beta1)) %>%
        mutate(
            mse_bias = mse_bias(coef(cox_mod, s = x), beta1),
            lambda = x,
            params = list(coef(cox_mod, s = x))
        )
}) %>%
    mutate(
        model = "enet-CR",
        sim = i,
        model_size = TP + FP
    )

# Quick plot (ROC)
table_result_cox %>%
    ggplot(aes(x = 1 - Specificity, y = Sensitivity)) +
    geom_path() +
    coord_equal()


## 4. cbSCRIP - enet) ----

# Fit the cbSCRIP model
tic("cv.lambdaAcc")
set.seed(123)
cv.lambdaAcc <- cbSCRIP(
    Surv(ftime, fstatus) ~ .,
    train,
    coeffs = "original"
)
toc()

# Process cbSCRIP results
table_result_cb <- map_df(seq_along(cv.lambdaAcc$lambdagrid), \(x) {
    as_tibble(varsel_perc(cv.lambdaAcc$coefficients[[x]][1:length(beta1), 1],
                          beta1)) %>%
        mutate(
            mse_bias = mse_bias(cv.lambdaAcc$coefficients[[x]][1:length(beta1), 1], beta1),
            lambda = cv.lambdaAcc$lambdagrid[x],
            params = list(cv.lambdaAcc$coefficients[[x]][1:length(beta1), 1])
        )
}) %>%
    mutate(
        model = "enet-casebase-Acc",
        sim = i,
        model_size = TP + FP
    )

# Plot coefficient paths
plot(cv.lambdaAcc, plot_intercept = F) +
    geom_vline(xintercept = cv.lambdaAcc$lambdagrid,
               alpha = 0.5)

# Quick plots (Sensitivity vs. Model Size / Specificity)
table_result_cb %>%
    ggplot(aes(x = model_size, y = Sensitivity)) +
    geom_path()

table_result_cb %>%
    ggplot(aes(x = 1 - Specificity, y = Sensitivity)) +
    geom_path() +
    coord_equal()



# Combine the results from the active models
table_select <- as_tibble(rbind(
    table_result_cox,
    table_result_cb
))

# Create final ROC plot
(plot <- table_select %>%
        ggplot(aes(
            x = 1 - Specificity,
            y = Sensitivity,
            color = model,
            linetype = model
        )) +
        geom_path() +
        coord_equal())

# Save plot
ggsave(
    here(
        "paper",
        "results",
        glue("models_selec-{setting}_iter-{i}_p-{p}_k-{num_true}.png")
    ),
    plot = plot,
    bg = "transparent",
    width = 200,
    height = 120,
    units = "mm",
    dpi = 300
)

# Save results object
saveRDS(
    table_select,
    here(
        "paper",
        "results",
        "lambda_paths",
        glue("models_selec-{setting}_iter-{i}_p-{p}_k-{num_true}.rds")
    )
)