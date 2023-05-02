#!/bin/bash

# Set the directory containing the Rmarkdown files (current folder)
RMD_DIR="."

# Loop through all Rmarkdown files in the directory
for rmd_file in "${RMD_DIR}"/*.Rmd; do
  # Get the output file name (HTML by default)
  output_file="${rmd_file%.Rmd}.html"
  
  # Run the Rmarkdown file and save the output to the same folder
  Rscript -e "rmarkdown::render('${rmd_file}', output_file='${output_file}')"
done


