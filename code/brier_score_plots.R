source("paper/code/brier_cif_fns.R")

#threshold for beta
thr.b <- 1e-6

#generate title

settings_tbl <- tibble::tribble(
    ~setting,    ~desc,
    "setting1",  "1. Single effects on end-point of interest",
    "setting2", "2. Single effects on both end-points",
    "setting3", "3. Opposing effects",
    "setting4", "4. Mixture of effects",
    "setting5", "5. Non-proportional hazards"
)

#----------------------------
# Params

img_dir <- here::here("paper","brier_figs")
res_dir <- here::here("paper","figs")

#---------------------
N <- 400; p <- 120; k <- 20

for (set in 1:5){
    brier_fn(N, p, k, set)
}

# grid of figures
files <- file.path(img_dir, sprintf("brier_setting%d_N%d_p%d_k%d.png", 1:5, N, p, k))
out_name   <- file.path(res_dir, sprintf("brier_grid_N%d_p%d_k%d.png", N, p, k))

plot_grid(files,out_name)
dev.off()  

#---------------------
N <- 400; p <- 500; k <- 20

for (set in 1:5){
    brier_fn(N, p, k, set)
}

# grid of figures
files <- file.path(img_dir, sprintf("brier_setting%d_N%d_p%d_k%d.png", 1:5, N, p, k))
out_name   <- file.path(img_dir, sprintf("brier_grid_N%d_p%d_k%d.png", N, p, k))

files <- file.path(img_dir, sprintf("brier_setting%d_N%d_p%d_k%d.png", 1:5, N, p, k))
out_name   <- file.path(res_dir, sprintf("brier_grid_N%d_p%d_k%d.png", N, p, k))

plot_grid(files,out_name)
dev.off()   

#----------------------------
N <- 400; p <- 500; k <- 84

for (set in 1:5){
    brier_fn(N, p, k, set)
}

# grid of figures
files <- file.path(img_dir, sprintf("brier_setting%d_N%d_p%d_k%d.png", 1:5, N, p, k))
out_name   <- file.path(res_dir, sprintf("brier_grid_N%d_p%d_k%d.png", N, p, k))

plot_grid(files,out_name)
dev.off()  

#----------------------------
N <- 400; p <- 1000; k <- 20

for (set in 1:5){
    brier_fn(N, p, k, set)
}

# grid of figures
files <- file.path(img_dir, sprintf("brier_setting%d_N%d_p%d_k%d.png", 1:5, N, p, k))
out_name   <- file.path(res_dir, sprintf("brier_grid_N%d_p%d_k%d.png", N, p, k))

plot_grid(files,out_name)
dev.off()  

#----------------------------------
N <- 400; p <- 1000; k <- 168

for (set in 1:5){
    brier_fn(N, p, k, set)
}

# grid of figures
files <- file.path(img_dir, sprintf("brier_setting%d_N%d_p%d_k%d.png", 1:5, N, p, k))
out_name   <- file.path(res_dir, sprintf("brier_grid_N%d_p%d_k%d.png", N, p, k))

plot_grid(files,out_name)
dev.off()   
