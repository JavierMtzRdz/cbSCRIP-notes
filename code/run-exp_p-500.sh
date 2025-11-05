#!/bin/bash

# Script to run an R script multiple times with a different index each time.

echo "Starting analysis loop..."

# --- Configuration ---
# Start and end of index loop.
START_INDEX=1
END_INDEX=30

# Path to R script.
R_SCRIPT_NAME="paper/code/fit-curve-lambda.R"

# --- Additional Arguments ---
# Define static P_ARG.
P_ARG=500

# Calculate the dynamic NUM_TRUE
TEMP_VAL=$(echo "$P_ARG * 0.04166667" | bc)
ROUNDED_VAL=$(printf "%.0f" "$TEMP_VAL")
CALCULATED_NUM_TRUE=$((ROUNDED_VAL * 4))

# Define the array of NUM_TRUE 
NUM_TRUE_VALUES=(20 $CALCULATED_NUM_TRUE)

# --- Loop Execution ---
for nt in "${NUM_TRUE_VALUES[@]}"
do
  echo "================================="
  echo "=== Processing for NUM_TRUE: $nt ==="
  echo "================================="

  # Middle loop to iterate through settings from 1 to 5.
  for s in $(seq 1 5)
  do
    echo "#################################"
    echo "### Processing for SETTING: $s ###"
    echo "#################################"

    # Innermost 'for' loop to iterate from START_INDEX to END_INDEX for each setting.
    for i in $(seq $START_INDEX $END_INDEX)
    do
      # Print which iteration is currently running.
      echo "---------------------------------"
      echo "Running iteration with index: $i (NUM_TRUE: $nt, Setting: $s)"
      echo "---------------------------------"
		
      # Execute the R script, passing the current NUM_TRUE ($nt) from the outer loop.
      # The order of arguments matters and must match what the R script expects.
      Rscript "$R_SCRIPT_NAME" "$P_ARG" "$nt" "$s" "$i"

      # Check if the R script ran successfully.
      if [ $? -eq 0 ]; then
        echo "Iteration $i completed successfully."
      fi
    done
  done
done

echo "---------------------------------"
echo "All iterations completed."
echo "---------------------------------"

