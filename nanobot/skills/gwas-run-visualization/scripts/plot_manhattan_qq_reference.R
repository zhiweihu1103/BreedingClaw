#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})

options(scipen = 9999)

usage <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript scripts/plot_manhattan_qq_reference.R",
      "--prefix TRAIT_PREFIX",
      "[--sig-file PREFIX.sig.Manhattan.pmap.txt]",
      "[--qq-file PREFIX.all.QQ.pmap.txt]",
      "[--outdir plots]",
      sep = " "
    ),
    "\n"
  )
}

parse_args <- function(args) {
  parsed <- list(
    prefix = NULL,
    sig_file = NULL,
    qq_file = NULL,
    outdir = "."
  )

  index <- 1
  while (index <= length(args)) {
    key <- args[[index]]
    value <- if (index < length(args)) args[[index + 1]] else NULL

    if (key %in% c("-h", "--help")) {
      usage()
      quit(save = "no", status = 0)
    }

    if (is.null(value) || startsWith(value, "--")) {
      stop("Missing value for argument: ", key, call. = FALSE)
    }

    if (key == "--prefix") {
      parsed$prefix <- value
    } else if (key == "--sig-file") {
      parsed$sig_file <- value
    } else if (key == "--qq-file") {
      parsed$qq_file <- value
    } else if (key == "--outdir") {
      parsed$outdir <- value
    } else {
      stop("Unknown argument: ", key, call. = FALSE)
    }

    index <- index + 2
  }

  if (is.null(parsed$prefix)) {
    stop("--prefix is required", call. = FALSE)
  }

  if (is.null(parsed$sig_file)) {
    parsed$sig_file <- paste0(parsed$prefix, ".sig.Manhattan.pmap.txt")
  }

  if (is.null(parsed$qq_file)) {
    parsed$qq_file <- paste0(parsed$prefix, ".all.QQ.pmap.txt")
  }

  parsed
}

require_columns <- function(data, required, label) {
  missing_cols <- setdiff(required, colnames(data))
  if (length(missing_cols) > 0) {
    stop(
      label,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
}

add_optional_columns <- function(data, columns) {
  for (column in columns) {
    if (!column %in% colnames(data)) {
      data[[column]] <- NA
    }
  }
  data
}

write_tsv <- function(data, path) {
  write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

sample_background_points <- function(data, keep_mask, max_background_points = 200000L, order_cols = NULL) {
  keep_mask <- as.logical(keep_mask)
  keep_mask[is.na(keep_mask)] <- FALSE

  keep_df <- data[keep_mask, , drop = FALSE]
  background_df <- data[!keep_mask, , drop = FALSE]

  if (nrow(background_df) > max_background_points) {
    set.seed(123)
    background_df <- background_df[sample.int(nrow(background_df), max_background_points), , drop = FALSE]
  }

  out <- rbind(keep_df, background_df)
  if (!is.null(order_cols) && nrow(out) > 0) {
    ordering <- do.call(order, unname(out[, order_cols, drop = FALSE]))
    out <- out[ordering, , drop = FALSE]
  }
  out
}

build_manhattan_plot <- function(data, chr_info, suggestive_line, bonferroni_line, prefix) {
  palette <- c(
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
    "#ff7f00", "#addd33", "#a65628", "#f781bf", "#999999"
  )

  ggplot(data, aes(x = x, y = log10p)) +
    geom_point(aes(color = factor(CHR)), size = 1.5) +
    scale_color_manual(values = rep_len(palette, length(unique(data$CHR)))) +
    scale_x_continuous(
      breaks = chr_info$marker,
      labels = paste0("Chr", chr_info$CHR),
      expand = c(0.01, 0)
    ) +
    geom_hline(yintercept = suggestive_line, linetype = "dashed", linewidth = 0.6) +
    geom_hline(
      yintercept = bonferroni_line,
      linetype = "dashed",
      linewidth = 0.6,
      color = "red"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      title = paste("GWAS Manhattan Plot -", prefix),
      x = "Chromosome",
      y = expression(-log[10](italic(P)))
    )
}

build_qq_plot <- function(data, prefix) {
  ggplot(data, aes(x = -log10(expected_p), y = -log10(P))) +
    geom_point(size = 1.2, color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      title = paste("Q-Q Plot -", prefix),
      x = expression(Expected ~ -log[10](italic(P))),
      y = expression(Observed ~ -log[10](italic(P)))
    )
}

downsample_qq <- function(data, max_points = 200000, tail_fraction = 0.01) {
  n_total <- nrow(data)
  if (n_total <= max_points) {
    return(data)
  }

  tail_count <- max(1, ceiling(n_total * tail_fraction))
  tail_count <- min(tail_count, n_total, max_points)
  remaining <- max_points - tail_count
  tail_index <- seq_len(tail_count)

  if (remaining <= 0 || tail_count == n_total) {
    return(data[tail_index, , drop = FALSE])
  }

  set.seed(123)
  body_index <- sample((tail_count + 1):n_total, remaining)
  data[sort(c(tail_index, body_index)), , drop = FALSE]
}

save_plot_set <- function(manhattan_plot, qq_plot, prefix, outdir) {
  combined_plot <- plot_grid(manhattan_plot, qq_plot, ncol = 2, rel_widths = c(4, 1.5))

  ggsave(
    filename = file.path(outdir, paste0(prefix, ".manhattan_qq.pdf")),
    plot = combined_plot,
    width = 12,
    height = 6,
    units = "in"
  )
  ggsave(
    filename = file.path(outdir, paste0(prefix, ".manhattan_qq.png")),
    plot = combined_plot,
    width = 12,
    height = 6,
    units = "in",
    dpi = 300
  )
  ggsave(
    filename = file.path(outdir, paste0(prefix, ".manhattan.pdf")),
    plot = manhattan_plot,
    width = 10,
    height = 6,
    units = "in"
  )
  ggsave(
    filename = file.path(outdir, paste0(prefix, ".manhattan.png")),
    plot = manhattan_plot,
    width = 10,
    height = 6,
    units = "in",
    dpi = 300
  )
  ggsave(
    filename = file.path(outdir, paste0(prefix, ".qq.pdf")),
    plot = qq_plot,
    width = 5,
    height = 5,
    units = "in"
  )
  ggsave(
    filename = file.path(outdir, paste0(prefix, ".qq.png")),
    plot = qq_plot,
    width = 5,
    height = 5,
    units = "in",
    dpi = 300
  )
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  if (!file.exists(args$sig_file)) {
    stop("Significant Manhattan input not found: ", args$sig_file, call. = FALSE)
  }
  if (!file.exists(args$qq_file)) {
    stop("QQ input not found: ", args$qq_file, call. = FALSE)
  }

  dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

  cat("[INFO] Trait prefix:", args$prefix, "\n")
  cat("[INFO] Significant Manhattan input:", args$sig_file, "\n")
  cat("[INFO] QQ input:", args$qq_file, "\n")

  manhattan_sig <- read.delim(args$sig_file, header = TRUE, sep = "\t", check.names = FALSE)
  qq_data <- read.delim(args$qq_file, header = TRUE, sep = "\t", check.names = FALSE)

  require_columns(manhattan_sig, c("SNP", "CHR", "BP", "P"), "Significant Manhattan input")
  require_columns(qq_data, c("P"), "QQ input")

  qq_data <- add_optional_columns(qq_data, c("SNP", "CHR", "BP"))

  manhattan_sig$CHR <- as.numeric(manhattan_sig$CHR)
  manhattan_sig$BP <- as.numeric(manhattan_sig$BP)
  manhattan_sig$P <- as.numeric(manhattan_sig$P)
  qq_data$P <- as.numeric(qq_data$P)
  qq_data$CHR <- suppressWarnings(as.numeric(qq_data$CHR))
  qq_data$BP <- suppressWarnings(as.numeric(qq_data$BP))

  valid_nonzero_p <- qq_data$P[!is.na(qq_data$P) & qq_data$P > 0]
  if (length(valid_nonzero_p) == 0) {
    stop("QQ input does not contain any non-zero p-values", call. = FALSE)
  }

  min_nonzero_p <- min(valid_nonzero_p)
  cat("[INFO] Minimum non-zero p-value:", min_nonzero_p, "\n")

  manhattan_sig$P[manhattan_sig$P == 0] <- min_nonzero_p
  qq_data$P[qq_data$P == 0] <- min_nonzero_p

  manhattan_sig <- manhattan_sig[!is.na(manhattan_sig$CHR) & !is.na(manhattan_sig$BP) & !is.na(manhattan_sig$P), ]
  qq_data <- qq_data[!is.na(qq_data$P), ]

  if (all(c("SNP", "CHR", "BP") %in% colnames(qq_data)) && any(!is.na(qq_data$CHR)) && any(!is.na(qq_data$BP))) {
    manhattan_all <- qq_data[!is.na(qq_data$CHR) & !is.na(qq_data$BP), c("SNP", "CHR", "BP", "P"), drop = FALSE]
  } else {
    manhattan_all <- manhattan_sig
  }

  chr_info <- aggregate(BP ~ CHR, manhattan_all, max)
  chr_info <- chr_info[order(chr_info$CHR), , drop = FALSE]
  chr_info$offset <- c(0, cumsum(chr_info$BP)[-nrow(chr_info)])
  chr_info$marker <- chr_info$offset + chr_info$BP / 2

  manhattan_all <- manhattan_all[order(manhattan_all$CHR, manhattan_all$BP), , drop = FALSE]
  manhattan_all$offset <- chr_info$offset[match(manhattan_all$CHR, chr_info$CHR)]
  manhattan_all$x <- manhattan_all$BP + manhattan_all$offset
  manhattan_all$log10p <- -log10(manhattan_all$P)

  suggestive_line <- -log10(1 / nrow(qq_data))
  bonferroni_line <- -log10(0.05 / nrow(qq_data))
  cat("[INFO] Suggestive line:", suggestive_line, "\n")
  cat("[INFO] Bonferroni line:", bonferroni_line, "\n")

  manhattan_suggestive <- subset(manhattan_all, log10p >= suggestive_line & log10p < bonferroni_line)
  manhattan_significant <- subset(manhattan_all, log10p >= bonferroni_line)
  manhattan_plot_data <- sample_background_points(
    data = manhattan_all,
    keep_mask = (manhattan_all$P <= 1e-4) | (manhattan_all$SNP %in% manhattan_sig$SNP),
    max_background_points = 200000L,
    order_cols = c("CHR", "BP")
  )

  write_tsv(
    manhattan_all[, c("SNP", "CHR", "BP", "P", "log10p")],
    file.path(args$outdir, paste0(args$prefix, ".manhattan.all.tsv"))
  )
  write_tsv(
    manhattan_suggestive[, c("SNP", "CHR", "BP", "P", "log10p")],
    file.path(args$outdir, paste0(args$prefix, ".manhattan.suggestive.tsv"))
  )
  write_tsv(
    manhattan_significant[, c("SNP", "CHR", "BP", "P", "log10p")],
    file.path(args$outdir, paste0(args$prefix, ".manhattan.significant.tsv"))
  )

  qq_data <- qq_data[order(qq_data$P), , drop = FALSE]
  qq_data$expected_p <- seq_len(nrow(qq_data)) / nrow(qq_data)
  qq_plot_data <- downsample_qq(qq_data)

  write_tsv(
    qq_plot_data[, c("SNP", "CHR", "BP", "P", "expected_p")],
    file.path(args$outdir, paste0(args$prefix, ".qq.plot.tsv"))
  )

  manhattan_plot <- build_manhattan_plot(
    data = manhattan_plot_data,
    chr_info = chr_info,
    suggestive_line = suggestive_line,
    bonferroni_line = bonferroni_line,
    prefix = args$prefix
  )
  qq_plot <- build_qq_plot(qq_plot_data, args$prefix)

  save_plot_set(manhattan_plot, qq_plot, args$prefix, args$outdir)
  cat("[INFO] Manhattan and QQ plotting completed\n")
}

main()
