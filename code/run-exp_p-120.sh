#!/bin/bash

# S script to run an R script multiple times with a different index each time.

echo "Starting analysis loop..."

# --- Configuration ---
# Start and end of index loop.
START_INDEX=1
END_INDEX=30

# Define path
R_SCRIPT_NAME="paper/code/fit-curve-lambda.R"

# --- Additional Arguments ---
# Define static arguments that will be passed to the R script 
P_ARG=120
NUM_TRUE=20

# Create the output directory if it doesn't exist
# mkdir -p "$OUTPUT_DIR"

# --- Loop Execution ---
# Outer loop to iterate through settings from 1 to 5.
for s in $(seq 1 5)
do
  echo "#################################"
  echo "### Processing for SETTING: $s ###"
  echo "#################################"

  # Inner 'for' loop to iterate from START_INDEX to END_INDEX for each setting.
  for i in $(seq $START_INDEX $END_INDEX)
  do
    # Print which iteration is currently running.
    echo "---------------------------------"
    echo "Running iteration with index: $i (Setting: $s)"
    echo "---------------------------------"

    Rscript "$R_SCRIPT_NAME" "$P_ARG" "$NUM_TRUE" "$s" "$i"


    if [ $? -eq 0 ]; then
      echo "Iteration $i completed successfully."
    fi
  done
done

echo "---------------------------------"
echo "All iterations completed."
echo "---------------------------------"
