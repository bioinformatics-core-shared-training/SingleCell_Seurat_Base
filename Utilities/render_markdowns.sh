#!/usr/bin/env bash

# Usage:
# process all markdowns:
#   bash Utilities/render_markdowns.sh                 
# process specific markdowns by number:
#   bash Utilities/render_markdowns.sh 04 10a 11 

# make sure this is being run from the project directory
if [ ! -f "SingleCell_Seurat_Base.Rproj" ]; then
  echo "Please run this script from the project directory"
  exit 1
fi

# build file list based on provided numbers, or default to all
if [ $# -gt 0 ]; then
  RMD_FILES=()
  for num in "$@"; do
    for f in Markdowns/${num}_*.Rmd; do
      [ -f "$f" ] && RMD_FILES+=("$f")
    done
  done
  if [ ${#RMD_FILES[@]} -eq 0 ]; then
    echo "No matching Rmd files found for: $*"
    exit 1
  fi
else
  RMD_FILES=(Markdowns/*.Rmd)
fi

# render HTML
echo "=== Rendering Markdowns to HTML ==="
for i in "${RMD_FILES[@]}"; do
  echo "Rendering $i..."
  Rscript -e "rmarkdown::render('$i')"
done

# Purl scripts for demos and exercises
echo ""
echo "=== Purling R scripts ==="
for i in "${RMD_FILES[@]}"; do
  output_file="TeachingScripts/$(basename "${i%.Rmd}.R")"
  echo "Purling $i..."
  Rscript -e "knitr::purl('$i', documentation = 0, output = '$output_file')"
done

echo ""
echo "Done."
