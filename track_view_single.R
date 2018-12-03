#!/usr/bin/Rscript

##------------------------------------------------------##
## IMPORTANT !!!
##
## import bigWig, only for TE (current)
##
## method
## 1. read bigWig file as GRanges object by rtracklayer::import.bw()
## support filter conditions: chromosome, star, end
##
## 2. normalize IP by input sample, in GRanges object format by
## trackViewer::GRoperator()
##
## 3. Convert GRanges object to data.frame, 1-base solution
##
## 4. Generate coverage plot using ggplot2 with data.frame data
##
## 5. multiple plots in one page
##
## option-1: make plots for single bigWig file
##
## option-2: normalize IP by input, and make plots for single bigWig file
##
## option-3: make plots for dual bigWig files, wt vs mut
##
## option-4: normalize each IP by correspond input.
##
##------------------------------------------------------##

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  print("Usage: track_view.R <ip.bw> <input.bw> <chr.size> <n>")
  print("")
  print("Option:")
  print("ip.bw    the bigWig file for IP sample")
  print("input.bw the bigWig file for Input sample")
  print("pdf_out  the pdf file to save plots")
  print("chr.size the chrom size of chromosomes")
  print("n        the list of TEs to display, default: TRUE")
  stop("argument failed")
}

ip_bw    <- args[1]
input_bw <- args[2]
pdf_out  <- args[3]
fasize   <- args[4]
if (length(args) > 4) {
  n <- unlist(strsplit(args[5], " "))
} else {
  n <- TRUE
}

# packages
if (! require(goldclipReport)) {
  stop("Need to install goldclipReport")
}
library(goldclipReport)

if (! require(rtracklayer)) {
  BiocInstaller::biocLite("rtracklayer")
  library(rtracklayer)
}

if (! require(trackViewer)) {
  BiocInstaller::biocLite("trackViewer")
  library(trackViewer)
}

if (! require(GenomicRanges)) {
  BiocInstaller::biocLite("GenomicRanges")
  library(GenomicRanges)
}

if (! require(dplyr)) {
  install.packages("dplyr")
  library(dplyr)
}

if (! require(ggplot2)) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (! require(ggridges)) {
  install.packages("ggridges")
  library(ggridges)
}

if (! require(tidyr)) {
  install.packages("tidyr")
  library(tidyr)
}

##----------------------------------------------------------------------------##
## List of ChIPseq data sets
df_list <- chipseq_bw_parser(ip_bw, fasize, input_bw, n)

plist <- lapply(df_list, function(d){
  coverage_plot_single(d, fill.color = "orange",
                      exclude.minus.scores = TRUE)
})

plot_n_pages(plist, nrow = 5, ncol = 2, pdf_out = pdf_out)



## EOF
