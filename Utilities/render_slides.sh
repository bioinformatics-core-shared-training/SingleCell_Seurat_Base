#!/usr/bin/env bash

# This script renders the Rmd slides in the Slides directory to HTML and then
# converts them to PDF using pagedown::chrome_print.
# Requires R, rmarkdown, and pagedown packages to be installed
#
# Usage:
#   bash Utilities/render_slides.sh [options] [deck_numbers...]
#
# Options:
#   --pdf-only       Skip rendering Rmd to HTML, only convert existing HTML to PDF
#   --help           Show this help message
#
# Examples:
#   bash Utilities/render_slides.sh                  # render and convert all slides
#   bash Utilities/render_slides.sh --pdf-only        # convert all HTML to PDF only
#   bash Utilities/render_slides.sh 03 07             # render and convert slides 03 and 07
#   bash Utilities/render_slides.sh --pdf-only 03 07  # convert slides 03 and 07 to PDF only

# make sure this is being run from the project directory
if [ ! -f "SingleCell_Seurat_Base.Rproj" ]; then
  echo "Please run this script from the project directory"
  exit 1
fi

# parse arguments
PDF_ONLY=false
DECK_NUMBERS=()

for arg in "$@"; do
  case "$arg" in
    --pdf-only)
      PDF_ONLY=true
      ;;
    --help)
      sed -n '/#/p' "$0" | head -20
      exit 0
      ;;
    *)
      DECK_NUMBERS+=("$arg")
      ;;
  esac
done

# build file lists based on deck numbers or default to all
if [ ${#DECK_NUMBERS[@]} -gt 0 ]; then
  RMD_FILES=()
  HTML_FILES=()
  for num in "${DECK_NUMBERS[@]}"; do
    for f in Slides/${num}*.Rmd; do
      [ -f "$f" ] && RMD_FILES+=("$(basename "$f")")
    done
    for f in Slides/${num}*.html; do
      [ -f "$f" ] && HTML_FILES+=("$f")
    done
  done
else
  RMD_FILES=($(cd Slides && ls *.Rmd 2>/dev/null))
  HTML_FILES=($(ls Slides/*.html 2>/dev/null))
fi

# render HTML slides from Rmd
if [ "$PDF_ONLY" = false ]; then
  echo "=== Rendering Slides to HTML ==="
  for rmd in "${RMD_FILES[@]}"; do
    echo "Rendering $rmd..."
    conda run -n RNAseq Rscript -e "rmarkdown::render('Slides/$rmd')"
  done
  # refresh HTML file list to include any newly rendered files
  if [ ${#DECK_NUMBERS[@]} -gt 0 ]; then
    HTML_FILES=()
    for num in "${DECK_NUMBERS[@]}"; do
      for f in Slides/${num}*.html; do
        [ -f "$f" ] && HTML_FILES+=("$f")
      done
    done
  else
    HTML_FILES=($(ls Slides/*.html 2>/dev/null))
  fi
fi

# convert HTML slides to PDF using pagedown::chrome_print
echo ""
echo "=== Converting Slides to PDF ==="
HTML_LIST=$(printf "'%s'," "${HTML_FILES[@]}")
HTML_LIST="c(${HTML_LIST%,})"
conda run -n RNAseq Rscript -e "
  files <- ${HTML_LIST}
  for (f in files) {
    out <- sub('html$', 'pdf', f)
    cat('Converting', f, '...\n')
    tryCatch(
      pagedown::chrome_print(f, output = out),
      error = function(e) cat('Error:', conditionMessage(e), '\n')
    )
  }
"

echo ""
echo "Done."
