#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

rmarkdown::render(
    "src/report_template.Rmd",
    output_file = paste0("../", args[1]),
    params = list(rds_paths = args[2:length(args)])
)
