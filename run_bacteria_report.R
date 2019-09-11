#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("At least one argument must be supplied (bacteria_dir)", call. = FALSE)
}

bacDir <- args[1]
pdfOut <- file.path(bacDir, "bacteria_content.pdf")

# search kraken2.stat file
statFiles <- list.files(bacDir, "*kraken2.stat",
                        all.files = TRUE,
                        full.names = TRUE)

# files
if (length(statFiles) == 0) {
  stop("No *.kraken2.stat file found", call. = FALSE)
}

# start program
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
library(pheatmap)


makeplot <- function(x, pdfOut, width = 8, height = 8) {
  # csv file
  csvFile <- gsub(".pdf$", ".csv", pdfOut)

  # read file
  # RPM
  dl <- lapply(x, function(i){
    f_name <- gsub(".kraken2.stat", "", basename(i))
    readr::read_delim(i, "\t", col_types = readr::cols()) %>%
      mutate(reads_in_tax = reads_in_tax * 1e6 / reads_total) %>%
      dplyr::select(name, reads_in_tax) %>%
      dplyr::rename(!! f_name := reads_in_tax)
  })

  df <- goldclipReport::merge_list_of_dataframes(dl, by = "name")
  # df <- df %>% dplyr::select(sort(colnames(df)))
  df <- tibble::column_to_rownames(df, "name")
  ma <- log10(df + 1)

  ## save to table
  write.csv(df, csvFile, row.names = TRUE, quote = FALSE)

  ## make plot
  library(RColorBrewer)
  cc = colorRampPalette(rev(brewer.pal(n = 7,
                                       name = "RdYlBu")))
  breaks = seq(0, 6, length.out = 100)

  library(pheatmap)
  p <- pheatmap(ma,
                silent = TRUE,
                cluster_cols = FALSE,
                color = cc(100),
                breaks = breaks,
                border_color = "grey40")
  pdf(pdfOut, width = width, height = height)
  print(p)
  dev.off()
}

makeplot(statFiles, pdfOut, width = 8, height = 8)

## EOF
