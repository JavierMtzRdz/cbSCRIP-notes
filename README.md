# cbSCRIP: Case-Base Sampling Competing Risks Penalized Regression

`cbSCRIP` is an R package for high-dimensional variable selection and predictive modeling under competing risks. It implements a parametric continuous-time cause-specific hazard model using multinomial logistic regression under discrete time-domain case-base sampling.

---

## Repository Structure

The project is structured to guarantee a clean, memory-safe reproducible workflow for our extensive simulation studies:

```
cbSCRIP-notes/
├── cbSCRIP/                 # R package source (SAGA solver in C++, multinomial case-base fitting, CV)
├── code/                    # Reproducible scripts
│   ├── run_all_sims.sh      # Master bash wrapper to sequence all experiments
│   ├── run-exp_p-*.sh       # Bash dispatchers for specific predictor dimensions (120, 500, 1000)
│   ├── 01_simulations.R     # Core single-replicate simulation script
│   ├── 02_evaluation.R      # Result aggregation, metric calculation, and plotting (Brier, AUC, C-index, CSR)
│   ├── 03_benchmarks.R      # Empirical timing benchmark script for optimizer comparison
│   └── 04_real_data.R       # Bladder cancer clinical data preprocessing and biomarker discovery
├── paper/                   # Manuscript materials (cbSCRIP.tex, WileyNJDv5.cls, reviewer PDFs, compiled cbSCRIP.pdf)
├── figs/                    # Generated high-resolution figures used by manuscript
├── refs/                    # Bibliography database (competing-risk.bib)
├── suppl/                   # Quarto files for supplementary reports
└── Makefile                 # Master workflow orchestrator
```

---

## Reproducing the Experiments

Because the competing risk baseline models and `cbSCRIP` penalization grid paths are computationally intensive, the outer simulation loop is managed via bash scripts. This memory-safe architecture dispatches independent R sessions for every single replicate. It guarantees memory stability (preventing R memory leaks across the 100 replicates) and will not overwhelm your local compute resources.

### Quick Start with Make
The repository includes a `Makefile` that orchestrates the entire workflow. You can use standard `make` commands to run specific segments of the pipeline:

```bash
make sims       # Run all memory-safe simulations (calls run_all_sims.sh)
make eval       # Aggregate metrics and generate plots
make benchmarks # Run optimizer timing benchmarks
make data       # Run real-world dataset analysis
make paper      # Compile the LaTeX manuscript PDF
make clean      # Clean up LaTeX intermediate files
make all        # Run the entire pipeline
```

## Installation

To install the `cbSCRIP` R package locally:

```R
# From within the repository
remotes::install_github("JavierMtzRdz/cbSCRIP")
```
