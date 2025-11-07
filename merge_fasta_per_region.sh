#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------------
# merge_fasta_per_region.sh
# Combine individual FASTA files per haploblock region into merged FASTAs.
#
# Usage:
#   ./merge_fasta_per_region.sh <input_directory> <output_directory> [--clean]
#
# Example:
#   ./merge_fasta_per_region.sh out_dir/TNFa/tmp/consensus_fasta out_dir/TNFa/haploblock_phased_seq_merged --clean
# ----------------------------------------------------------------------

# --- Argument checks ---------------------------------------------------
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory> [--clean]" >&2
    exit 1
fi

input_dir="$1"
output_dir="$2"
clean_flag="${3:-}"

# --- Safety checks -----------------------------------------------------
if [ ! -d "$input_dir" ]; then
    echo "Error: input directory '$input_dir' does not exist." >&2
    exit 1
fi

# --- Cleanup old output (optional) ------------------------------------
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

# --- Logging setup -----------------------------------------------------
echo "Starting FASTA merge..."
echo "Input:  $input_dir"
echo "Output: $output_dir"
echo

# --- Merge loop --------------------------------------------------------
fasta_files=$(find "$input_dir" -type f \( -name "*.fa" -o -name "*.fasta" \) | sort)

if [ -z "$fasta_files" ]; then
    echo "No FASTA files found in '$input_dir'." >&2
    exit 1
fi

for f in $fasta_files; do
    name_noext=$(basename "$f" | sed 's/\.[^.]*$//')

    # Extract region from filename (e.g., chr6_region_711055-761032)
    region=$(echo "$name_noext" | grep -oE 'chr[0-9XYM]+_region_[0-9]+-[0-9]+' || true)
    [ -z "$region" ] && region="unknown_region"

    output="${output_dir}/${region}.fa"

    # Extract header from filename (e.g., HG001_chr6_region_..._hap1)
    header=$(echo "$name_noext" | grep -oE '[^_]+_chr[0-9XYM]+_region_[0-9]+-[0-9]+_hap[0-9]+' || true)
    [ -z "$header" ] && header="$name_noext"

    {
        echo ">$header"
        grep -v "^>" "$f"
    } >> "$output"
done

# --- Optional parallel version (commented-out for perspective) ---------
# find "$input_dir" -type f \( -name "*.fa" -o -name "*.fasta" \) | \
#   parallel --will-cite -j 4 '
#     f={};
#     name_noext=$(basename "$f" | sed "s/\.[^.]*$//");
#     region=$(echo "$name_noext" | grep -oE "chr[0-9XYM]+_region_[0-9]+-[0-9]+" || echo "unknown_region");
#     output="${output_dir}/${region}.fa";
#     header=$(echo "$name_noext" | grep -oE "[^_]+_chr[0-9XYM]+_region_[0-9]+-[0-9]+_hap[0-9]+" || echo "$name_noext");
#     { echo ">$header"; grep -v "^>" "$f"; } >> "$output";
#   '

# --- Archive output (commented-out perspective) ------------------------
# timestamp=$(date +"%Y%m%d_%H%M%S")
# tar -czf "${output_dir}_${timestamp}.tar.gz" -C "$(dirname "$output_dir")" "$(basename "$output_dir")"
# echo "Archive created: ${output_dir}_${timestamp}.tar.gz"

echo
echo "Merge complete."
echo "One FASTA per region written to: ${output_dir}/"

# --- Optional cleanup of tmp directory (commented for safety) ----------
# echo "Removing temporary directory: ${input_dir}"
# rm -rf "$input_dir"

