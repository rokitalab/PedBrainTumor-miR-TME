#!/bin/bash

# Bicna Song
set -e
set -o pipefail

# Set directory paths
root_dir=$(git rev-parse --show-toplevel)
output_dir="results"

# Set file paths
target_fasta="$output_dir/gencode.v39.3utr.fa"
mirna_fasta="$output_dir/combined_all_novel_mirna.fa"
output_file="$output_dir/miranda_output.txt"

# Run miRanda
echo "Running miRanda target prediction..."
start_time=$(date +%s)

miranda "$mirna_fasta" "$target_fasta" -out "$output_file"

end_time=$(date +%s)
runtime=$((end_time - start_time))

echo "miRanda finished. Output saved to: $output_file"
echo "Runtime: $runtime seconds"
