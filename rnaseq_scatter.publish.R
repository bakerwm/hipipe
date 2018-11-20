

## scatter plot
setwd("~/work/yu_2018/projects/20181004_chipseq_zk/hipipe-chip-output/20181010-scatter-plot")
source("rnaseq_plot_functions.R")
library(goldclipReport)
library(ggplot2)
library(ggrepel)

annotation.labels <- c("HeT-A", "blood", "mdg1", "Burdock")
organism = "dm3"

# attp2 vs dNxf2 KD
x <- "/home/wangming/work/yu_2018/projects/20180720_RNAseq_nxf2/results/RNAseq_rerun/shAttp2_vs_shNxf2/transposon_analysis/de_analysis/RNAseq_shattp2_vs_RNAseq_shNxf2/transcripts_deseq2.csv"
pdf1 <- "RNAseq_attp2_vs_dNxf2.te.pdf"
x.label <- "dNxf2 wt"
y.label <- "dNxf2 KD"
p <- make_deseq2_scatter_plot2(x, organism = "dm3",
                              x.label = x.label, y.label = y.label,
                              show.sig.labels = FALSE,
                              annotation.labels = annotation.labels,
                              pvalue_cutoff = 0.05, save2pdf = NULL)
pdf(pdf1, width = 2.5, height = 2.5, paper = "a4")
print(p)
dev.off()


# CG9754 KD vs dNxf2 KD
pdf2 <- "RNAseq_CG9754_vs_dNxf2.te.pdf"
x.label <- "CG9754 KD"
y.label <- "dNxf2 KD"
x1 <- "/home/wangming/work/yu_2018/projects/20180720_RNAseq_nxf2/results/RNAseq_rerun/dsWhite_vs_dsCG9754/transposon_analysis/de_analysis/RNAseq_dsWhite_vs_RNAseq_dsCG9754/transcripts_deseq2.csv"
x2 <- "/home/wangming/work/yu_2018/projects/20180720_RNAseq_nxf2/results/RNAseq_rerun/shAttp2_vs_shNxf2/transposon_analysis/de_analysis/RNAseq_shattp2_vs_RNAseq_shNxf2/transcripts_deseq2.csv"
df1 <- DESeq2_csv2df(x1)
df2 <- DESeq2_csv2df(x2)
df_norm <- merge(dplyr::select(df1, 1:3), dplyr::select(df2, 1:3), by = "id") %>%
      dplyr::mutate(RNAseq_dsCG9754_norm = RNAseq_dsCG9754 * RNAseq_shattp2 / RNAseq_dsWhite) %>%
  dplyr::select(id, RNAseq_dsCG9754_norm, RNAseq_shNxf2)
names(df_norm) <- c("id", "x", "y")
df_norm$id <- as.character(df_norm$id)
# convert to log10
# replace zero values by 0.001 (rpm)
df_norm$x[df_norm$x == 0] <- 0.001
df_norm$y[df_norm$y == 0] <- 0.001
df_norm <- dplyr::mutate(df_norm, x = log10(x), y = log10(y))

p2 <- publish_plot_de_scatter(df_norm, x.label, y.label, point.size = 1.2,
                             point.color = "grey50", sig.labels = NULL,
                             sig.size = 1.2, sig.color = "red",
                             annotation.labels = annotation.labels,
                             annotation.point.size = 1.2,
                             annotation.point.color = "red",
                             save2pdf = FALSE)
pdf(pdf2, width = 2.5, height = 2.5, paper = "a4")
print(p2)
dev.off()




################################################################################
# dNxf2 het vs dNxf2 mut
setwd("~/work/yu_2018/projects/20181004_chipseq_zk/hipipe-chip-output/20181010-scatter-plot/")
pdf3 <- "RNAseq_dNxf2_wt_vs_dNxf2_mut.te.pdf"
x <- "/home/wangming/work/yu_2018/projects/20181004_chipseq_zk/results/w1118_vs_Nxf2_22aamut_25124/transposon_analysis/de_analysis/RNAseq_W1118_vs_RNAseq_NXF2_22aamut_25124/transcripts_deseq2.csv"
x.label <- "dNxf2 wt"
y.label <- "dNxf2 mut"
annotation.labels <- c("HeT-A", "blood", "mdg1")
p <- make_deseq2_scatter_plot2(x, organism = "dm3",
                               x.label = x.label, y.label = y.label,
                               show.sig.labels = TRUE,
                               annotation.labels = annotation.labels,
                               pvalue_cutoff = 0.05, save2pdf = NULL)

x.title <- expression(paste(italic("W1118 "), "[", log["10"], " rpm]"))
y.title <- expression(paste(italic("dNxf2 mut "), "[", log["10"], " rpm]"))
p <- p + xlab(x.title) + ylab(y.title)

# pdf(pdf3, width = 3, height = 3, paper = "a4")
# print(p)
# dev.off()

## significant TE in nxf2_wt_vs_nxf2_mut
df_nxf2_mut <- DESeq2_csv2df(x)
figF_sig_TE <- dplyr::filter(df_nxf2_mut, padj < 0.05)

pdf4 <- "RNAseq_CG9754_vs_dNxf2.te.pdf"
x.label <- "Panx mut"
y.label <- "dNxf2 mut"
x1 <- "~/work/yu_2018/projects/20181004_chipseq_zk/results/w1118_vs_Nxf2_22aamut_25124/transposon_analysis/de_analysis/RNAseq_W1118_vs_RNAseq_NXF2_22aamut_25124/transcripts_deseq2.csv"
x2 <- "~/work/yu_2018/projects/20181004_chipseq_zk/results/CG9754_het_vs_mut/transposon_analysis/de_analysis/RNAseq_CG9754_het_vs_RNAseq_CG9754_mut/transcripts_deseq2.csv"
df1 <- DESeq2_csv2df(x1)
df2 <- DESeq2_csv2df(x2)
df_norm <- merge(dplyr::select(df1, 1:3), dplyr::select(df2, 1:3), by = "id") %>%
  dplyr::mutate(panx_mut_norm = RNAseq_CG9754_mut * RNAseq_W1118 / RNAseq_CG9754_het) %>%
  dplyr::select(id, panx_mut_norm, RNAseq_NXF2_22aamut_25124)
names(df_norm) <- c("id", "x", "y")
df_norm$id <- as.character(df_norm$id)

# fix NA and missing values
df_norm$x[is.infinite(df_norm$x)] = 0
df_norm$y[is.infinite(df_norm$y)] = 0

# convert to log10
# replace zero values by 0.001 (rpm)
df_norm$x[df_norm$x == 0] <- 0.001
df_norm$y[df_norm$y == 0] <- 0.001
df_norm <- dplyr::mutate(df_norm, x = log10(x), y = log10(y))

annotation.labels = figF_sig_TE$id

p2 <- publish_plot_de_scatter(df_norm, x.label, y.label, point.size = 1.2,
                              point.color = "grey50", sig.labels = NULL,
                              sig.size = 1.2, sig.color = "red",
                              annotation.labels = annotation.labels,
                              annotation.point.size = 1.2,
                              annotation.point.color = "red",
                              show.labels.text = FALSE,
                              save2pdf = FALSE)

x.title <- expression(paste(italic("CG9754 mut "), "[", log["10"], " rpm]"))
y.title <- expression(paste(italic("dNxf2 mut "), "[", log["10"], " rpm]"))
p2 <- p2 + xlab(x.title) + ylab(y.title)

# pdf(pdf4, width = 3, height = 3, paper = "a4")
# print(p2)
# dev.off()

figureLegend <- "Figure. Differential express analysis between indicated samples.
a. significant expressed TEs (p value < 0.05) were labeled in red. (blood is on the right of
HeT-A in plot).
b. significant expressed TEs in dNxf2_mut vs dNxf2_wt were labeled in red in this plot."

library(cowplot)
pdf("extended_figure2fg-TE-scatter.pdf", width = 5.5, height = 6.3, paper = "a4")
pg1 <- cowplot::plot_grid(p, p2, NULL, NULL, align = "hv", ncol = 2, labels = "auto")
pg2 <- cowplot::add_sub(pg1, figureLegend, x = 0, hjust = 0, size = 10)
cowplot::ggdraw(pg2)
dev.off()
















#
# ## fix Zhang Yiqun
# setwd("~/work/yu_2018/projects/20181004_chipseq_zk/results_publish/")
# z1a <- "mut_vs_wt.TE.csv"
# z2a <- "panx_mut_vs_nxf_mut.TE.csv"
# df_z1a <- read.csv(z1a)
# df_z2a <- read.csv(z2a)
#
# # counts
#
#






