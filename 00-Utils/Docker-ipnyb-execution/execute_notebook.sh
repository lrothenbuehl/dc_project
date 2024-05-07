#!/bin/bash

# Check if a file path is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <notebook_file>"
    exit 1
fi

NOTEBOOK_PATH=$1

# Print the notebook path for debugging
echo "Notebook Path: $NOTEBOOK_PATH"

# Print the current working directory
echo "Current Directory: $(pwd)"

# List files in the current directory for debugging
echo "Files in Current Directory:"
ls -la

# Execute the notebook
jupyter nbconvert --to notebook --execute --inplace --output-dir=/output "$NOTEBOOK_PATH"

# Check if jupyter nbconvert command was successful
if [ $? -ne 0 ]; then
    echo "Failed to convert the notebook."
    exit 1
fi

# Ensure the output notebook is moved to a known location

