#!/bin/bash


MD5FILE="MD5"

# 
> "$MD5FILE"


for file in *; do
  if [[ -f "$file" && "$file" != *.gz && "$file" != *.md && "$file" != *.sh && "$file" != "$MD5FILE" ]]; then
  
    md5_before=$(md5 -q "$file")
    echo "Original: $file | MD5: $md5_before" >> "$MD5FILE"
    
    gzip "$file"
    
    md5_after=$(md5 -q "$file.gz")
    echo "Compressed: ${file}.gz | MD5: $md5_after" >> "$MD5FILE"
  fi
done

echo "MD5 checksums recorded in $MD5FILE"