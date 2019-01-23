#!/usr/bin/env Rscripts
# run DESeq2 analysis for matrix data
# matrix: featureCounts output
#
# Generate plots

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  print("Usage: Rscript run_deseq2.R <count.txt> <path_plot>")
  print("")
  print("Options:")
  print("  input     The matrix output of featureCounts, or path of DESeq2 analysis")
  print("  organism  The reference genome, [dm3, hg19]")
  print("  output    The directory to save plot file")
  stop("arguments failed")
}

# count.txt, path_out,
input     <- args[1]
organism  <- args[2]
output    <- args[3]

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DESeq2))
library(goldclipReport)
library(ggplot2)
library(ggrepel)

# run DESseq2 analysis
if(dir.exists(input)) {
  print("skip DESeq2 analysis")
  path_de <- input
} else if(file.exists(input)) {
  print("run DESeq2 analysis")
  path_de <- file.path(output, "de_analysis")
  DESeq2_for_featureCounts(input, organism, path_de, pvalue_cutoff = 0.1)
}

# Generate publish quality plots
fs <- list.files(path_de, "transcripts_deseq2.csv$", TRUE, TRUE, TRUE)
fs <- fs[1]
print(fs)

##----------------------------------------------------------------------------##
## save table
fs2 <- gsub(".csv$", ".fix.xls", fs)
if(file.exists(fs2)) {
  print("xls file exists")
} else {
  ## load gene names
  f_name <- paste("genelist", organism, "rda", sep = ".")
  f <- system.file("extdata", f_name, package = "goldclipReport")
  load(f) # genelist

  ## load matrix table
  df <- DESeq2_csv2df(fs)
  df2 <- dplyr::mutate(df, id = as.character(id)) %>%
    dplyr::mutate(id = plyr::mapvalues(id, genelist$gene_id, genelist$gene_name, FALSE)) %>%
    dplyr::filter(! is.na(padj)) %>%
    dplyr::arrange(padj)

  ## save table
  readr::write_delim(df2, fs2, delim = "\t", col_names = TRUE)
}

#-----------------------------------------------------------------------------##
print("generate publishable plots")
print(paste0("found DESeq2 ouptut: ", fs2))
path_pdf <- file.path(output, "report")
# DESeq2_publish_plot(fs2, path_pdf, gene_labels = c("mdg1", "blood", "HeT-A", "Burdock"), save2pdf = TRUE)
tmp <- DESeq2_publish_plot(fs2, path_pdf, save2pdf = TRUE)


# EOF

