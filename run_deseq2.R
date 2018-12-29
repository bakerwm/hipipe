#!/usr/bin/env Rscripts
# run DESeq2 analysis for matrix data
# matrix: featureCounts output
#
# Generate plots

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  print("Usage: Rscript run_deseq2.R <count.txt> <path_plot>")
  print("")
  print("Option:")
  print("  count.txt    The count of genes by featureCounts")
  print("  organism     The reference genome, [dm3, hg19]")
  print("  path_plot    The directory to save plot file")
  stop("arguments failed")
}

# count.txt, path_out,
file_count   <- args[1]
organism     <- args[2]
path_results <- args[3]

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DESeq2))
library(goldclipReport)
library(ggplot2)
library(ggrepel)
# file_count   <- "results/count/count.txt"
# path_results <- "results"

# run DESseq2 analysis
print("run DESeq2 analysis")
path_de <- file.path(path_results, "de_analysis")
DESeq2_for_featureCounts(file_count, organism, path_de, pvalue_cutoff = 0.1)

# Generate plots
fs <- list.files(path_de, "transcripts_deseq2.csv$", TRUE, TRUE, TRUE)
print(fs)
#
# ##----------------------------------------------------------------------------##
# ## save table
f_name = paste("genelist", organism, "rda", sep = ".")
f = system.file("extdata", f_name, package = "goldclipReport")
load(f) # genelist

df <- DESeq2_csv2df(fs)
df2 <- dplyr::mutate(df, id = as.character(id)) %>%
  dplyr::mutate(id = plyr::mapvalues(id, genelist$gene_id, genelist$gene_name, FALSE)) %>%
  dplyr::filter(! is.na(padj)) %>%
  dplyr::arrange(padj)
fs2 <- paste0(fs, ".fix.xls")
readr::write_delim(df2, fs2, delim = "\t", col_names = TRUE)
# ##----------------------------------------------------------------------------##

print("generate publishable plots")
print(paste0("found DESeq2 ouptut: ", fs[1]))
path_pdf <- file.path(path_results, "report")
tmp <- lapply(fs, function(f) {
  DESeq2_publish_plot(f, path_pdf, gene_labels = c("mdg1", "blood", "HeT-A", "Burdock"), save2pdf = TRUE)
})

# EOF

