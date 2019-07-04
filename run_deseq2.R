#!/usr/bin/env Rscripts
# run DESeq2 analysis for matrix data
# matrix: featureCounts output
#
# Generate plots

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  print("Usage: Rscript run_deseq2.R <count.txt> <path_plot>")
  print("")
  print("Options:")
  print("    control_txt*   the output of featureCounts of Control samples")
  print("  treatment_txt*   the output of featureCounts of Treatment samples")
  print("       organism*   the reference genome, [dm3, hg19]")
  print("         output*   the directory to save plot file")
  print("   control_name   the name for control samples")
  print(" treatment_name   the name for treatment samples")
  print("        p_value   the pvalue cutoff, default: 0.05")
  stop("arguments failed")
}


# count.txt, path_out,
countA    <- args[1]
countB    <- args[2]
organism  <- args[3]
output    <- args[4]
nameA     <- args[5]
nameB     <- args[6]
pvalue    <- args[7]


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DESeq2))
library(goldclipReport)
library(ggplot2)
library(ggrepel)

# run DESseq2 analysis
goldclipReport::deseqHub(countA, countB,
                         organism = organism,
                         nameA = nameA,
                         nameB = nameB,
                         outdir = output,
                         pvalue_cutoff = pvalue,
                         readable = FALSE)
path_de <- output

## Generate publish quality figures
cnt_fix <- file.path(output, "transcripts_deseq2.fix.xls")
stopifnot(file.exists(cnt_fix))

#-----------------------------------------------------------------------------##
print("generate publishable plots")
print(paste0("found DESeq2 ouptut: ", cnt_fix))
tmp <- DESeq2_publish_plot(cnt_fix, output, save2pdf = TRUE)

# EOF

