#!/usr/bin/env bash

# make sure this is being run from the project directory
if [ ! -f "SingleCell_Seurat_Base.Rproj" ]; then
  echo "Please run this script from the project directory"
  exit 1
fi

# render HTML
for i in Markdowns/*.Rmd
do
  Rscript -e "rmarkdown::render('$i')"
done

# Purl scripts for demos and exercises
for i in Markdowns/*.Rmd
do
  output_file="course_files/Demonstrations/$(basename "${i%.Rmd}.R")"
  
  Rscript -e "knitr::purl('$i', documentation = 0, output = '$output_file')"
done
