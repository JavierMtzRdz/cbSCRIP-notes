#!/bin/bash

# ==============================================================================
# Master Simulation Runner
# ==============================================================================

echo "Starting Full Simulation Suite..."

# Run p=120 experiments (Settings 1-5, Iterations 1-100)
echo "----------------------------------------"
echo "Dispatching p=120 experiments"
echo "----------------------------------------"
bash code/run-exp_p-120.sh

# Run p=500 experiments (Settings 1-5, Iterations 1-100)
echo "----------------------------------------"
echo "Dispatching p=500 experiments"
echo "----------------------------------------"
bash code/run-exp_p-500.sh

# Run p=1000 experiments (Settings 1-5, Iterations 1-100)
echo "----------------------------------------"
echo "Dispatching p=1000 experiments"
echo "----------------------------------------"
bash code/run-exp_p-1000.sh

echo "----------------------------------------"
echo "All Simulations Completed Successfully!"
echo "Next step: Run code/02_evaluation.R to aggregate results and plot metrics."
echo "----------------------------------------"
