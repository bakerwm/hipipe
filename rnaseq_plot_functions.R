


#' make high quality scatter plot
#'
#' DESeq2 analysis, padj < 0.05 for significant
#' count matrix to plot on x-axis and y-axis
#'
#' @param x path to the csv file of DESeq2 output, to parse the padj values
#'   require mean values in col-2, col-3
#'
#' @param annotation.labels a set of characters for annotation
#'
#' @param save2pdf A logical or charactervalue, whether or not save the
#'   plot to PDF file, default: FALSE
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export
#'
make_deseq2_scatter_plot2 <- function(x, organism = "dm3", #x.name, y.name,
                                      x.label = NULL, y.label = NULL,
                                      show.sig.labels = TRUE,
                                      annotation.labels = NULL,
                                      pvalue_cutoff = 0.05, save2pdf = NULL) {
  stopifnot(file.exists(x))
  
  # read DESeq2 csv file
  df_deseq2 <- DESeq2_csv2df(x)
  df_deseq2$id <- as.character(df_deseq2$id)
  df_sig <- dplyr::filter(df_deseq2, padj <= pvalue_cutoff)
  sig.labels <- df_sig$id
  if (! isTRUE(show.sig.labels)) {
    sig.labels <- NULL
  }
  
  # check file size
  if (length(df_deseq2) < 13) {
    stop("csv file columns not enough")
  }
  
  # prepare data table
  df_mean <- dplyr::select(df_deseq2, 1:3)
  smp <- names(df_mean)[-1]
  names(df_mean) <- c("id", "x", "y")
  # convert to log10
  # replace zero values by 0.001 (rpm)
  df_mean$x[df_mean$x == 0] <- 0.001
  df_mean$y[df_mean$y == 0] <- 0.001
  df_mean <- dplyr::mutate(df_mean, x = log10(x), y = log10(y))
  
  # prepare plot
  x.label <- ifelse(is.null(x.label), smp[1], x.label)
  y.label <- ifelse(is.null(y.label), smp[2], y.label)
  p <- publish_plot_de_scatter(df_mean, x.label, y.label, point.size = 1.2,
                               point.color = "grey50", sig.labels = sig.labels,
                               sig.size = 1.2, sig.color = "red",
                               annotation.labels = annotation.labels,
                               annotation.point.size = 1.2,
                               annotation.point.color = "red",
                               save2pdf = FALSE)
  return(p)
}







#' make high quality scatter plot
#'
#' DESeq2 analysis, padj < 0.05 for significant
#' count matrix to plot on x-axis and y-axis
#'
#' @param x path to the csv file of DESeq2 output, to parse the padj values
#'
#' @param f path to the matrix file, save the counts for each record,
#'   featureCounts, or HTSeq
#'
#' @param annotation.labels a set of characters for annotation
#'
#' @param save2pdf A logical or charactervalue, whether or not save the
#'   plot to PDF file, default: FALSE
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export
#'
make_deseq2_scatter_plot <- function(x, f, organism = "dm3", #x.name, y.name,
                                     x.label = NULL, y.label = NULL,
                                     show.sig.labels = TRUE,
                                     annotation.labels = NULL,
                                     pvalue_cutoff = 0.05, save2pdf = NULL) {
  stopifnot(file.exists(x) & file.exists(f))
  
  # read DESeq2 csv file
  df_deseq2 <- DESeq2_csv2df(x)
  df_sig <- dplyr::filter(df_deseq2, padj <= pvalue_cutoff)
  sig.labels <- df_sig$id
  if (! isTRUE(show.sig.labels)) {
    sig.labels <- NULL
  }
  
  # read featureCounts
  df <- featureCountsReader(f, normalizeTo1M = TRUE, fixZero = 1,
                            organism = "dm3")
  if (length(df) < 2) {
    stop("less than 2 BAM files of replicate samples")
  }
  
  # cal mean values for replicates
  smp <- unique(gsub("_rep[1-9]$", "", names(df)))
  if (length(smp) != 2) {
    stop("only accept 2 conditions")
  }
  
  # df means
  df_mean <- data.frame(id = rownames(df),
                        x  = rowMeans(dplyr::select(df, starts_with(smp[1]))),
                        y  = rowMeans(dplyr::select(df, starts_with(smp[2]))),
                        stringsAsFactors = FALSE)
  # convert to log10
  # replace zero values by 0.001 (rpm)
  df_mean$x[df_mean$x == 0] <- 0.001
  df_mean$y[df_mean$y == 0] <- 0.001
  df_mean <- dplyr::mutate(df_mean, x = log10(x), y = log10(y))
  
  
  # prepare plot
  x.label <- ifelse(is.null(x.label), smp[1], x.label)
  y.label <- ifelse(is.null(y.label), smp[2], y.label)
  p <- publish_plot_de_scatter(df_mean, x.label, y.label, point.size = 1.5,
                               point.color = "grey30", sig.labels = sig.labels,
                               sig.size = 1.5, sig.color = "red",
                               annotation.labels = annotation.labels,
                               annotation.point.size = 1.5,
                               annotation.point.color = "red",
                               save2pdf = FALSE)
  return(p)
}














#'
#' @param data the data.frame for plot, should contain "id" "x", "y",
#'   columns, x and y are values plot on x-axis and y-axis. id is the
#'   values to annotate.
#'
#' @param annotation.labels a set of ids for annotation, should match "id" column,
#'   default: NULL
#'
#' @param sig.labels a set of ids for the significant ids, which will be
#'   shown in red, larger size. default: NULL
#'
#' @param x.label label on x-axis,
#'   !!! for current publication, shuold like this: dsCG9754 [log10 rpm]
#'
#' @param y.label label on y-axis, the same as x.label
#'
#' @param save2df A logical or charactervalue, whether or not save the
#'   plot to PDF file, default: FALSE
#'
#'
#' @import ggplot2
#' @import ggrepel
#'
#' @export
#'
publish_plot_de_scatter <- function(data, x.label, y.label,
                                    point.size = 1,
                                    point.color = "grey30",
                                    sig.labels = NULL,
                                    sig.size = 1,
                                    sig.color = "red",
                                    annotation.labels = NULL,
                                    annotation.point.size = 1,
                                    annotation.point.color = "red",
                                    show.labels.text = TRUE,
                                    save2pdf = FALSE) {
  stopifnot(is.data.frame(data))
  stopifnot(all(c("id", "x", "y") %in% names(data)))
  stopifnot(is.character(data$id))
  stopifnot(is.numeric(data$x) & is.numeric(data$y))
  stopifnot(is.character(x.label) & is.character(y.label))
  
  # axis title
  x.title <- bquote(.(x.label)~"["*log["10"]~"rpm]")
  y.title <- bquote(.(y.label)~"["*log["10"]~"rpm]")
  
  # annotation
  if (! is.null(annotation.labels)) {
    df_annotation <- dplyr::filter(data, id %in% annotation.labels)
  }
  
  # significant points
  if (! is.null(sig.labels)) {
    df_sig <- dplyr::filter(data, id %in% sig.labels) # only sig
    data <- dplyr::filter(data, ! id %in% sig.labels) # exclude sig
  } else {
    df_sig <- NULL
  }
  
  # make basic plot
  if ("group" %in% names(df_sig)) {
    p <- ggplot(data, aes(x, y, color = group)) +
      geom_point(size = point.size)
  } else {
    p <- ggplot(data, aes(x, y)) +
      geom_point(color = point.color, size = sig.size)
  }
  # extra lines
  p <- p +
    xlab(x.title) + ylab(y.title) +
    scale_x_continuous(limits = c(0, 5),
                       breaks = seq(0, 5, 1),
                       labels = seq(0, 5, 1),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 5),
                       breaks = seq(0, 5, 1),
                       labels = seq(0, 5, 1),
                       expand = c(0, 0)) +
    theme_linedraw() +
    theme(panel.border = element_rect(color = "black", fill = NA, size = .7),
          plot.title   = element_text(color = "black", hjust = .5, size = 12),
          panel.grid   = element_blank(),
          axis.line    = element_blank(),
          axis.ticks.length = unit(.15, "cm"),
          axis.ticks   = element_line(color = "black", size = .5),
          axis.text    = element_text(color = "black", size = 10),
          axis.title   = element_text(color = "black", size = 12))
  
  # add sig points
  if (is.data.frame(df_sig)) {
    # specify colors
    if ("group" %in% names(df_sig)) {
      p <- p + geom_point(aes(x, y, color = group), data = df_sig, size = sig.size)
    } else {
      p <- p + geom_point(data = df_sig, color = sig.color, size = sig.size)
    }
  }
  
  # add annotation
  if (nrow(df_annotation) > 0) {
    if (isTRUE(show.labels.text)) {
      p <- p +
        geom_text_repel(
          data          = df_annotation,
          label         = df_annotation$id,
          size          = 4,
          box.padding   = .5,
          segment.size  = 0.5,
          segment.color = "grey50",
          direction     = "both")
    }
    p <- p +
      geom_point(data = df_annotation, color = annotation.point.color,
                 size = annotation.point.size)
  }
  
  # add ablines
  p <- p + geom_abline(slope = 1, intercept = 0, color = "black", size = .5) +
    geom_abline(slope = 1, intercept = log10(2), color = "black",
                linetype = 2, size = .5) +
    geom_abline(slope = 1, intercept = -log10(2), color = "black",
                linetype = 2, size = .5)
  
  # whether save to pdf
  if (is.character(save2pdf)) {
    pdf_dir <- dirname(savd2pdf)
    if (! dir.exists(dirname(save2pdf))) {
      warning("directory of pdf file not exists")
    } else {
      pdf(rpt_file, width = 4, height = 4, paper = "a4")
      print(p)
      dev.off()
      print(paste("save plot to file:", save2pdf))
    }
  }
  
  
  return(p)
}

