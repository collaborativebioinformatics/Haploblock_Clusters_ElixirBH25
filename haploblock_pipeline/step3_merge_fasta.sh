#!/usr/bin/env bash
set -euo pipefail

input_dir="$1"
output_dir="$2"
clean_flag="${3:-}"
threads="${4:-}"      # optional number of threads


[ ! -d "$input_dir" ] && { echo "Input dir not found"; exit 1; }

[ -d "$output_dir" ] && [ "$clean_flag" = "--clean" ] && rm -rf "$output_dir"
mkdir -p "$output_dir"

# Determine number of parallel jobs
if [[ -z "$threads" || "$threads" -le 0 ]]; then
    jobs=$(( $(nproc) - 1 ))
    jobs=$(( jobs > 0 ? jobs : 1 ))
else
    jobs="$threads"
fi

echo "Starting FASTA merge with $jobs jobs..."

merge_region() {
    local region="$1"
    local output_dir="$2"
    shift 2
    local files=("$@")
    output="${output_dir}/${region}.fa"
    > "$output"

    for f in "${files[@]}"; do
        name_noext=$(basename "$f" | sed -E 's/(\.vcf)?\.(fa|fasta)$//')
        header=$(echo "$name_noext" | grep -oE '[^_]+_chr[0-9XYM]+_region_[0-9]+-[0-9]+_hap[0-9]+' || true)
        [ -z "$header" ] && header="$name_noext"
        {
            echo ">$header"
            grep -v "^>" "$f"
        } >> "$output"
    done
}

export -f merge_region

# --- Prepare a list of regions with files ---
tmp_file=$(mktemp)
while IFS= read -r f; do
    name_noext=$(basename "$f" | sed -E 's/(\.vcf)?\.(fa|fasta)$//')
    region=$(echo "$name_noext" | grep -oE 'chr[0-9XYM]+_region_[0-9]+-[0-9]+' || true)
    [ -z "$region" ] && region="unknown_region"
    echo "$region $f" >> "$tmp_file"
done < <(find "$input_dir" -type f \( -name "*.fa" -o -name "*.fasta" \) | sort)

# --- Merge per region ---
cut -d' ' -f1 "$tmp_file" | sort -u | while read -r region; do
    files=($(awk -v r="$region" '$1==r {print $2}' "$tmp_file"))
    merge_region "$region" "$output_dir" "${files[@]}" &
    # Control max parallel jobs
    while (( $(jobs -r | wc -l) >= jobs )); do sleep 0.1; done
done

wait
rm "$tmp_file"

echo "Merge complete."

##########################
## Permissions 
##########################

sleep 200

chmod 644 "$output_dir"/*.fa
