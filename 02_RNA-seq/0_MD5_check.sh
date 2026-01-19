#!/bin/bash

set -euo pipefail

input_file="${input}/${tissue}_sample_name.txt"
output_file="${output}/${tissue}_md5output.txt"

#> "$output_file"

while read -r sample; do
  sample_dir="${path1}/${sample}"

  if [[ ! -d "$sample_dir" ]]; then
    echo "[ERROR] Directory not found: $sample_dir" >> "$output_file"
    continue
  fi

  if [[ ! -f "$sample_dir/MD5.txt" ]]; then
    echo "[ERROR] MD5.txt not found in $sample_dir" >> "$output_file"
    continue
  fi

  echo "### Checking $sample ###" >> "$output_file"
  (
    cd "$sample_dir" || exit
    md5sum -c MD5.txt
  ) >> "$output_file" 2>&1

done < "$input_file"