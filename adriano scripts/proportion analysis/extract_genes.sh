#!/bin/bash


if ! command -v python3 &> /dev/null; then
    echo "Error: Python3 is not installed or not in PATH."
    exit 1
fi
for file in *_expression_proportions_*_*regulated.txt; do
    if [[ -f "$file" ]]; then
        echo "Processing: $file"
        python3 extract_gene_and_musify.py "$file" ExonS_event_Cluster.txt 1
        if [[ $? -ne 0 ]]; then
            echo "Error processing $file"
        fi
    else
        echo "No matching files found"
    fi
done
