#!/usr/bin/env Rscripts

## make alignment stat

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2){
  print("Usage: Rscript alignment_stat.R <file|dir> <out_dir>")
  print("")
  print("Option:")
  print("  <file|dir>    The directory of alignment output, or a file save multiple stat.csv, each one in one-line")
  print("  <out_dir>     The directory where to save the results")
  stop("arguments failed")
}

input  <- args[1]
output <- args[2]
template <- NULL
if (file.exists(args[3])) {
  template <- args[3]
}

## check packages
if (! require("devtools")) {
  install.packages("devtools")
}

if (! require(goldclipReport)) {
  devtools::install_github("bakerwm/goldclipReport")
}

if (! require(dplyr)) {
  install.packages("dplyr")
}

library(goldclipReport)

## check arguments
if(dir.exists(input)) {
  stat_files <- list.files(input, "*mapping_stat.csv", all.files = TRUE,
                           full.names = TRUE, ignore.case = TRUE)
  stat_files <- normalizePath(stat_files)
} else if(file.exists(input)) {
 stat_files <- readr::read_lines(input, skip_empty_rows = TRUE)
} else {
  stat_files <- NULL
}

if(length(stat_files) == 0) {
  stop("stat.csv file not detected")
}

if (is.null(template)) {
  print("using default template")
  align_report(stat_files, output, preview = FALSE)
} else {
  print("custom template")
  align_report(stat_files, output, template = template, preview = FALSE)
}
## save
print(paste0("Saving results in ", output))

## EOF
