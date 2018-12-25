#!/usr/bin/env Rscripts

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1){
  print("Usage: Rscript qc_report.R <qc.dir>")
  print("")
  print("Option:")
  print("  qc.dir    The directory of fastqc output files")
  stop("arguments failed")
}

qc.dir <- args[1]

if (! dir.exists(qc.dir)) {
  stop("directory not exists")
}
qc.report <- file.path(qc.dir, "report")

if (! require("devtools")) {
  install.packages("devtools")
}

if (! require(goldclipReport)) {
  devtools::install_github("bakerwm/goldclipReport")
}

if (! require(fastqcr)) {
  devtools::install_github("bakerwm/fastqcr")
}

if (! require(dplyr)) {
  install.packages("dplyr")
}

library(goldclipReport)

FastqcReport(qc.dir, qc.report, preview = FALSE)

## save
print(paste0("Saving results in ", qc.report))

## EOF


