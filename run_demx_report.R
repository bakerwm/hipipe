#!/usr/bin/env Rscripts
# run demx report
#
# Generate plots

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  print("Usage: Rscript run_demx_report.R <path|results>")
  print("")
  print("Options:")
  print("    path|results*   path to demx directory")
  stop("arguments failed")
}

# count.txt, path_out,
# resPath <- args[1]

resPath <- normalizePath(args[1])

suppressPackageStartupMessages(library(dplyr))
library(goldclipReport)
library(ggplot2)
library(ggrepel)

# run DESseq2 analysis
goldclipReport::demxReport(resPath)

# EOF

