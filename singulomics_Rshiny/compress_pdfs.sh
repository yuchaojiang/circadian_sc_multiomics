#!/usr/bin/env bash

# Usage: ./compress_pdfs.sh /path/to/input_dir /path/to/output_dir quality
# quality: one of screen, ebook, printer, prepress (controls compression level)

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <input_dir> <output_dir> <quality> (screen|ebook|printer|prepress)"
  exit 1
fi

input_dir="$1"
output_dir="$2"
quality="/$3"

mkdir -p "$output_dir"

for file in "$input_dir"/*.pdf; do
  [[ -e "$file" ]] || { echo "No PDF files found in $input_dir"; break; }
  base=$(basename "$file" .pdf)
  out="${output_dir}/${base}_compressed.pdf"

  echo "Compressing $file → $out"
  gs -sDEVICE=pdfwrite \
     -dCompatibilityLevel=1.4 \
     -dPDFSETTINGS="$quality" \
     -dNOPAUSE -dQUIET -dBATCH \
     -sOutputFile="$out" "$file"

  if [[ -f "$out" ]]; then
    echo "  Done."
  else
    echo "  Failed to compress $file"
  fi
done

echo "All PDFs processed."
