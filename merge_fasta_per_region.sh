#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------------
# merge_fasta_per_region.sh
# Combine individual FASTA files per haploblock region into merged FASTAs.
#
# Usage:
#   ./merge_fasta_per_region.sh <input_directory> <output_directory> [--clean]
# ----------------------------------------------------------------------

# --- Argument checks ---
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory> [--clean]" >&2
    exit 1
fi

input_dir="$1"
output_dir="$2"
clean_flag="${3:-}"

# --- Safety checks ---
if [ ! -d "$input_dir" ]; then
    echo "Error: input directory '$input_dir' does not exist." >&2
    exit 1
fi

# --- Cleanup old output (optional) ---
if [ -d "$output_dir" ]; then
    echo "Output directory '$output_dir' already exists."
    if [ "$clean_flag" = "--clean" ]; then
        echo "Removing existing directory..."
        rm -rf "$output_dir"
    else
        echo "Refusing to overwrite. Use --clean to remove it." >&2
        exit 1
    fi
fi

mkdir -p "$output_dir"

# --- Logging ---
echo "Starting FASTA merge..."
echo "Input:  $input_dir"
echo "Output: $output_dir"
echo

# --- Merge loop (safe, no subshell issues) ---
found_any=false

# Use find + while loop with process substitution to avoid pipe subshell
while IFS= read -r f; do
    [ -e "$f" ] || continue
    found_any=true

    # Remove extension (.fa, .fasta)
    name_noext=$(basename "$f" | sed -E 's/\.(fa|fasta)$//')

    # Extract region from filename
    region=$(echo "$name_noext" | grep -oE 'chr[0-9XYM]+_region_[0-9]+-[0-9]+' || true)
    [ -z "$region" ] && region="unknown_region"

    output="${output_dir}/${region}.fa"

    # Extract header
    header=$(echo "$name_noext" | grep -oE '[^_]+_chr[0-9XYM]+_region_[0-9]+-[0-9]+_hap[0-9]+' || true)
    [ -z "$header" ] && header="$name_noext"

    {
        echo ">$header"
        grep -v "^>" "$f"
    } >> "$output"

done < <(find "$input_dir" -type f \( -name "*.fa" -o -name "*.fasta" \) | sort)

# --- Check if any FASTAs were found ---
if ! $found_any; then
    echo "No FASTA files found in '$input_dir'." >&2
    exit 1
fi

echo
echo "Merge complete."
echo "One FASTA per region written to: ${output_dir}/"

