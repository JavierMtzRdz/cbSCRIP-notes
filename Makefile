# ==============================================================================
# Makefile for cbSCRIP Reproducibility Pipeline
# ==============================================================================
# This Makefile orchestrates the entire workflow from simulation to paper compilation.
#
# Usage:
#   make all - Run everything (sims, eval, benchmark, data, paper)
#   make sims - Run all parallel memory-safe simulations (takes time)
#   make eval - Aggregate simulation results and plot metrics
#   make benchmarks  - Run empirical timing benchmarks
#   make data - Run real dataset analysis
#   make paper - Compile the LaTeX manuscript
#   make clean - Remove LaTeX intermediate build files
# ==============================================================================

.PHONY: all sims eval benchmarks data sims-custom eval-short suppl paper clean plots real-plots

all: sims eval benchmarks data suppl paper

# 1. Run all simulation experiments via the bash dispatcher
# Default is to overwrite existing simulation results
OVERWRITE ?= TRUE
export OVERWRITE

sims:
	@echo "Dispatching memory-safe simulations (OVERWRITE=$(OVERWRITE))..."
	bash code/run_all_sims.sh

# 2. Aggregate metrics and plot (Time-AUC, Brier, CSR)
eval:
	@echo "Evaluating simulation predictions and generating plots..."
	Rscript code/02_evaluation.R

# 3. Run computational timing benchmarks
benchmarks:
	@echo "Running timing benchmarks..."
	Rscript code/03_benchmarks.R

# 4. Run real-world data analysis (Bladder Cancer dataset)
data:
	@echo "Running real data analysis..."
	quarto render suppl/real-data-analysis.qmd --to pdf
	cp suppl/real-data-analysis.pdf paper/

# Run a custom number of simulation experiments
# Usage: make sims-custom REPS=2 SETTINGS="1" P_VALS="120" K_VALS="20"
REPS ?= 2
SETTINGS ?= 1
P_VALS ?= 120
K_VALS ?= $(K)

sims-custom:
	@echo "Running custom simulations: settings='$(SETTINGS)', p='$(P_VALS)', k='$(if $(K_VALS),$(K_VALS),auto)', reps=$(REPS)..."
	@for p in $(P_VALS); do \
		for s in $(SETTINGS); do \
			for i in $$(seq 1 $(REPS)); do \
				if [ -n "$(K_VALS)" ]; then \
					for k in $(K_VALS); do \
						echo "Running Setting $$s, P $$p, K $$k, Replicate $$i..."; \
						Rscript -e 'source("code/01_simulations.R")' $$p $$k $$s $$i; \
					done; \
				else \
					num_true=20; \
					if [ $$p -eq 500 ]; then num_true=84; fi; \
					if [ $$p -eq 1000 ]; then num_true=168; fi; \
					echo "Running Setting $$s, P $$p, K $$num_true, Replicate $$i..."; \
					Rscript -e 'source("code/01_simulations.R")' $$p $$num_true $$s $$i; \
				fi; \
			done; \
		done; \
	done

# Evaluate the custom simulation run
eval-short:
	@echo "Evaluating custom simulation runs..."
	Rscript code/short_eval.R "$(shell echo $(SETTINGS) | tr ' ' ',')" "$(shell echo $(P_VALS) | tr ' ' ',')" "$(REPS)" "$(shell echo $(K_VALS) | tr ' ' ',')"

# 5. Targets for generating plots and compiling supplementary materials
.PHONY: plots real-plots suppl

plots:
	@echo "Generating simulation plots and rendering supplement..."
	Rscript code/selection-plots.r
	Rscript code/prediction-plots.r

suppl:
	@echo "Generating real dataset plots and rendering supplement..."
	quarto render suppl/preproccessing-bladder-dataset.qmd --to pdf
	quarto render suppl/real-data-analysis.qmd --to pdf
	cp suppl/preproccessing-bladder-dataset.pdf suppl/real-data-analysis.pdf submission/
	quarto render suppl/Brier-Score_CIF_plots.qmd --to pdf
	cp suppl/Brier-Score_CIF_plots.pdf submission/

# 6. Compile the LaTeX manuscript in the submission/ directory
paper:
	@echo "Compiling LaTeX manuscript..."
	cd submission && pdflatex -interaction=nonstopmode cbSCRIP.tex
	cd submission && bibtex cbSCRIP
	cd submission && pdflatex -interaction=nonstopmode cbSCRIP.tex
	cd submission && pdflatex -interaction=nonstopmode cbSCRIP.tex
	@echo "Manuscript compiled successfully: submission/cbSCRIP.pdf"

# Clean up LaTeX artifacts (keeping the final PDF)
clean:
	@echo "Cleaning up LaTeX build artifacts..."
	rm -f submission/*.aux submission/*.log submission/*.out submission/*.bbl submission/*.blg submission/*.synctex.gz
