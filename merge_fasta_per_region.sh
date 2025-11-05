#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_directory> <output_directory>"
    exit 1
fi

input_dir="$1"
output_dir="$2"

# Check directories
if [ ! -d "$input_dir" ]; then
    echo "Error: input directory '$input_dir' does not exist."
    exit 1
fi

mkdir -p "$output_dir"

for f in "$input_dir"/*.fasta "$input_dir"/*.fa; do
    [ -e "$f" ] || continue

    # Remove extension
    name_noext=$(basename "$f" | sed 's/\.[^.]*$//')

    # Extract region from filename (e.g., chr6_region_711055-761032)
    region=$(echo "$name_noext" | grep -oE 'chr[0-9XYM]+_region_[0-9]+-[0-9]+')

    # Fallback if no region found
    [ -z "$region" ] && region="unknown_region"

    # Output file for this region
    output="${output_dir}/${region}.fa"

    # Extract header from filename
    header=$(echo "$name_noext" | grep -oE '[^_]+_chr[0-9XYM]+_region_[0-9]+-[0-9]+_hap[0-9]+')
    [ -z "$header" ] && header="$name_noext"

    echo ">$header" >> "$output"
    grep -v "^>" "$f" >> "$output"
done

tar -zcvf "${output_dir}.tar.gz" "${output_dir}"

echo "✅ One FASTA per region written to: ${output_dir}"
