#!/bin/bash

# Check if a notebook file is provided as an argument
if [ -z "$1" ]; then
  echo "Usage: $0 <notebook-file.ipynb>"
  exit 1
fi

# Assign the first argument to a variable
NOTEBOOK_FILE=$1

# Check if the file exists
if [ ! -f "$NOTEBOOK_FILE" ]; then
  echo "Error: File '$NOTEBOOK_FILE' not found!"
  exit 1
fi

# Define output notebook name
OUTPUT_NOTEBOOK="output_$(basename $NOTEBOOK_FILE)"

# Execute the notebook using papermill
papermill "$NOTEBOOK_FILE" "$OUTPUT_NOTEBOOK"

# Check if papermill execution was successful
if [ $? -eq 0 ]; then
  echo "Notebook executed successfully. Output saved to '$OUTPUT_NOTEBOOK'."
else
  echo "Error: Failed to execute the notebook."
  exit 1
fi

