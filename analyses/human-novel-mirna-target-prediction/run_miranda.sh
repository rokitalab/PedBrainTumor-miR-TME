#!/bin/bash

# Bicna Song
set -e
set -o pipefail

# Set directory paths
root_dir=$(git rev-parse --show-toplevel)
data_dir="$root_dir/data"

# Set file paths
target_fasta_gz="$data_dir/gencode.v44.transcripts.fa.gz"
mirna_fasta="input/de_mirna_sequences.fa"

# Output directory and file
output_dir="results"
output_file="${output_dir}/miranda_output.txt"

# Create output directory if it doesn’t exist
mkdir -p "$output_dir"

# Create a temporary file for uncompressed target FASTA
tmp_fa=$(mktemp /tmp/gencode.XXXXXX.fa)
echo "Decompressing target FASTA to temp file: $tmp_fa"
zcat "$target_fasta_gz" > "$tmp_fa"

# Run miRanda
echo "Running miRanda target prediction..."
start_time=$(date +%s)

miranda "$mirna_fasta" "$tmp_fa" -out "$output_file"

end_time=$(date +%s)
runtime=$((end_time - start_time))

echo "miRanda finished. Output saved to: $output_file"
echo "Runtime: $runtime seconds"

# Clean up temp file
rm "$tmp_fa"
echo "Cleaned up temp file: $tmp_fa"