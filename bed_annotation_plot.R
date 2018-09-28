#!/usr/bin/Rscript


args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 2) {
    fn <- args[1]
    fout <- args[2]
} else {
    print("Usage: bed_annotation_plot.R <anno.stat> <out.pdf>")
    stop('At lease two argument')
}

library(goldclipReport)
library(dplyr)
library(ggplot2)
df <- AnnoParser(fn, dm3_basic = TRUE)
p <- BedAnnoBarplot(df, type="percentage", category='peaks')

# define width
# width (inch): 3 - 6
nsmp <- length(unique(df$id))
w <- nsmp + 2
w <- ifelse(w > 6, 6, w)
pdf(fout, width = w, height = 4, paper = 'a4')
print(p)
dev.off()