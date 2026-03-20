#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) return(y)
  if (length(x) == 1 && (is.na(x) || identical(x, ""))) return(y)
  x
}

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  defaults <- list(
    trait_label = "Trait",
    chrom_sizes = "",
    global_manhattan = "",
    qtl_regions = "",
    local_manhattan = "",
    gene_models = "",
    highlight_genes = "",
    out = "ld_qtl_summary.pdf",
    y_min = "3",
    flank_bp = "200000",
    max_loci = "0",
    genomewide_p = as.character(0.05 / 4033947)
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Bad argument: ", key)
    name <- gsub("-", "_", substring(key, 3))
    if (!(name %in% names(defaults))) stop("Unknown argument: ", key)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      defaults[[name]] <- "true"
      i <- i + 1L
    } else {
      defaults[[name]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  defaults$y_min <- as.numeric(defaults$y_min)
  defaults$flank_bp <- as.integer(defaults$flank_bp)
  defaults$max_loci <- as.integer(defaults$max_loci)
  defaults$genomewide_p <- as.numeric(defaults$genomewide_p)
  defaults
}

read_tsv <- function(path) {
  if (is.null(path) || identical(path, "") || !file.exists(path)) return(data.frame())
  info <- file.info(path)
  if (is.na(info$size) || info$size == 0) return(data.frame())
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )
}

ensure_cols <- function(df, cols, label) {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

alpha_col <- function(col, alpha = 1) {
  grDevices::adjustcolor(col, alpha.f = alpha)
}

sample_background_points <- function(df, keep_mask, max_background_points = 200000L, order_cols = NULL) {
  keep_mask <- as.logical(keep_mask)
  keep_mask[is.na(keep_mask)] <- FALSE

  keep_df <- df[keep_mask, , drop = FALSE]
  background_df <- df[!keep_mask, , drop = FALSE]

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

r2_bin <- function(r2) {
  if (is.na(r2)) return("#d1d1d1")
  if (r2 >= 0.8) return("#b2182b")
  if (r2 >= 0.6) return("#ef8a62")
  if (r2 >= 0.4) return("#fddbc7")
  if (r2 >= 0.2) return("#67a9cf")
  "#d1d1d1"
}

infer_chrom_sizes <- function(global_df, qtl_df) {
  candidates <- data.frame()
  if (nrow(global_df) > 0) {
    tmp <- aggregate(bp ~ chrom, data = global_df, FUN = max)
    colnames(tmp) <- c("chrom", "length_bp")
    candidates <- rbind(candidates, tmp)
  }
  if (nrow(qtl_df) > 0) {
    tmp <- aggregate(panel_end ~ chrom, data = qtl_df, FUN = max)
    colnames(tmp) <- c("chrom", "length_bp")
    candidates <- rbind(candidates, tmp)
  }
  if (nrow(candidates) == 0) {
    stop("Unable to infer chromosome sizes. Provide chrom_sizes or valid GWAS/QTL coordinates.")
  }
  out <- aggregate(length_bp ~ chrom, data = candidates, FUN = max)
  out <- out[order(as.numeric(out$chrom)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

build_offsets <- function(chrom_sizes) {
  chrom_sizes <- chrom_sizes[order(as.numeric(chrom_sizes$chrom)), , drop = FALSE]
  offsets <- numeric(nrow(chrom_sizes))
  names(offsets) <- as.character(chrom_sizes$chrom)
  if (nrow(chrom_sizes) > 1) {
    offsets[2:nrow(chrom_sizes)] <- cumsum(chrom_sizes$length_bp[1:(nrow(chrom_sizes) - 1L)])
  }
  centers <- offsets + chrom_sizes$length_bp / 2
  list(offsets = offsets, centers = centers)
}

normalize_inputs <- function(opts) {
  global_df <- read_tsv(opts$global_manhattan)
  qtl_df <- read_tsv(opts$qtl_regions)
  local_df <- read_tsv(opts$local_manhattan)
  gene_df <- read_tsv(opts$gene_models)
  highlight_df <- read_tsv(opts$highlight_genes)
  chrom_sizes <- read_tsv(opts$chrom_sizes)

  ensure_cols(global_df, c("snp_id", "chrom", "bp"), "global_manhattan")
  ensure_cols(qtl_df, c("locus_id", "chrom", "lead_snp_id", "lead_bp", "qtl_start", "qtl_end"), "qtl_regions")
  ensure_cols(local_df, c("locus_id", "snp_id", "chrom", "bp"), "local_manhattan")

  if (!("mlog10_p" %in% colnames(global_df))) {
    ensure_cols(global_df, c("p_value"), "global_manhattan")
    global_df$mlog10_p <- -log10(pmax(as.numeric(global_df$p_value), .Machine$double.xmin))
  }
  if (!("p_value" %in% colnames(global_df))) {
    global_df$p_value <- 10^(-as.numeric(global_df$mlog10_p))
  }
  if (!("lead_p" %in% colnames(qtl_df))) {
    qtl_df$lead_p <- 10^(-as.numeric(qtl_df$lead_mlog10_p))
  }
  if (!("mlog10_p" %in% colnames(local_df))) {
    ensure_cols(local_df, c("p_value"), "local_manhattan")
    local_df$mlog10_p <- -log10(pmax(as.numeric(local_df$p_value), .Machine$double.xmin))
  }
  if (!("p_value" %in% colnames(local_df))) {
    local_df$p_value <- 10^(-as.numeric(local_df$mlog10_p))
  }
  if (!("r2" %in% colnames(local_df))) local_df$r2 <- NA_real_
  if (!("is_lead" %in% colnames(local_df))) {
    local_df$is_lead <- local_df$snp_id == qtl_df$lead_snp_id[match(local_df$locus_id, qtl_df$locus_id)]
  }
  if (!("panel_start" %in% colnames(qtl_df))) qtl_df$panel_start <- pmax(1L, qtl_df$qtl_start - opts$flank_bp)
  if (!("panel_end" %in% colnames(qtl_df))) qtl_df$panel_end <- qtl_df$qtl_end + opts$flank_bp
  if (!("qtl_type" %in% colnames(qtl_df))) qtl_df$qtl_type <- "qtl"

  global_df$chrom <- as.integer(global_df$chrom)
  global_df$bp <- as.numeric(global_df$bp)
  global_df$mlog10_p <- as.numeric(global_df$mlog10_p)

  qtl_df$chrom <- as.integer(qtl_df$chrom)
  qtl_df$lead_bp <- as.numeric(qtl_df$lead_bp)
  qtl_df$qtl_start <- as.numeric(qtl_df$qtl_start)
  qtl_df$qtl_end <- as.numeric(qtl_df$qtl_end)
  qtl_df$panel_start <- as.numeric(qtl_df$panel_start)
  qtl_df$panel_end <- as.numeric(qtl_df$panel_end)
  qtl_df$lead_p <- as.numeric(qtl_df$lead_p)

  local_df$chrom <- as.integer(local_df$chrom)
  local_df$bp <- as.numeric(local_df$bp)
  local_df$mlog10_p <- as.numeric(local_df$mlog10_p)
  local_df$r2 <- as.numeric(local_df$r2)

  if (nrow(chrom_sizes) == 0) {
    chrom_sizes <- infer_chrom_sizes(global_df, qtl_df)
  } else {
    ensure_cols(chrom_sizes, c("chrom", "length_bp"), "chrom_sizes")
    chrom_sizes$chrom <- as.integer(chrom_sizes$chrom)
    chrom_sizes$length_bp <- as.numeric(chrom_sizes$length_bp)
  }

  if (nrow(gene_df) > 0) {
    ensure_cols(gene_df, c("gene_id", "chrom", "start_bp", "end_bp", "strand"), "gene_models")
    gene_df$chrom <- as.integer(gene_df$chrom)
    gene_df$start_bp <- as.numeric(gene_df$start_bp)
    gene_df$end_bp <- as.numeric(gene_df$end_bp)
  }
  if (nrow(highlight_df) > 0) {
    ensure_cols(highlight_df, c("locus_id", "gene_id"), "highlight_genes")
  }

  list(
    chrom_sizes = chrom_sizes,
    global = global_df,
    qtl = qtl_df,
    local = local_df,
    genes = gene_df,
    highlights = highlight_df
  )
}

layout_lanes <- function(genes, pad = 15000L) {
  if (nrow(genes) == 0) return(integer())
  genes <- genes[order(genes$start_bp, genes$end_bp, genes$gene_id), , drop = FALSE]
  ends <- numeric()
  out <- integer(nrow(genes))
  names(out) <- genes$gene_id
  for (i in seq_len(nrow(genes))) {
    lane <- 0L
    while (lane < length(ends) && genes$start_bp[[i]] <= ends[[lane + 1L]] + pad) {
      lane <- lane + 1L
    }
    if (lane == length(ends)) {
      ends <- c(ends, genes$end_bp[[i]])
    } else {
      ends[[lane + 1L]] <- genes$end_bp[[i]]
    }
    out[[genes$gene_id[[i]]]] <- lane
  }
  out
}

plot_global_panel <- function(trait_label, global_df, loci, offsets, centers, chrom_sizes, y_min, genomewide_p) {
  global_df$x <- global_df$bp + offsets[as.character(global_df$chrom)]
  global_df$c <- ifelse(global_df$chrom %% 2L == 1L, "#4C78A8", "#F58518")
  keep_mask <- global_df$p_value <= 1e-4
  if (nrow(loci) > 0) {
    for (i in seq_len(nrow(loci))) {
      locus <- loci[i, , drop = FALSE]
      keep_mask <- keep_mask |
        (global_df$chrom == locus$chrom &
           global_df$bp >= locus$qtl_start &
           global_df$bp <= locus$qtl_end)
    }
  }
  plot_df <- sample_background_points(
    global_df,
    keep_mask = keep_mask,
    max_background_points = 200000L,
    order_cols = c("chrom", "bp")
  )
  ymax <- max(global_df$mlog10_p, na.rm = TRUE)
  if (nrow(loci) > 0) {
    ymax <- max(ymax, -log10(pmax(min(loci$lead_p, na.rm = TRUE), .Machine$double.xmin)))
  }
  ylim <- c(y_min, max(y_min + 0.1, ymax * 1.08))
  xlim <- c(0, sum(chrom_sizes$length_bp))

  plot(xlim, ylim, type = "n", xaxt = "n", xlab = "", ylab = "-log10(P)", main = "", cex.main = 0.95)
  usr <- par("usr")
  for (i in seq_len(nrow(loci))) {
    locus <- loci[i, , drop = FALSE]
    rect(
      offsets[[as.character(locus$chrom)]] + locus$qtl_start,
      usr[3],
      offsets[[as.character(locus$chrom)]] + locus$qtl_end,
      usr[4],
      col = alpha_col(locus$plot_color, 0.12),
      border = NA
    )
  }
  points(plot_df$x, plot_df$mlog10_p, pch = 16, cex = 0.35, col = alpha_col(plot_df$c, 0.75))
  abline(h = -log10(genomewide_p), col = "#b2182b", lty = 2, lwd = 0.8)
  for (i in seq_len(nrow(loci))) {
    locus <- loci[i, , drop = FALSE]
    lx <- offsets[[as.character(locus$chrom)]] + locus$lead_bp
    ly <- -log10(pmax(locus$lead_p, .Machine$double.xmin))
    points(lx, ly, pch = 21, bg = locus$plot_color, col = "black", cex = 0.9)
    text(lx, ly + 0.15, labels = locus$locus_id, cex = 0.65, pos = 3)
  }
  axis(
    1,
    at = centers[as.character(chrom_sizes$chrom)],
    labels = paste0("Chr", chrom_sizes$chrom),
    cex.axis = 0.75
  )
  title(main = sprintf("%s | LD/QTL summary", trait_label), cex.main = 0.95)
}

plot_local_panel <- function(locus, local_df, y_min) {
  local <- local_df[
    local_df$locus_id == locus$locus_id &
      local_df$bp >= locus$panel_start &
      local_df$bp <= locus$panel_end,
    ,
    drop = FALSE
  ]

  if (nrow(local) == 0) {
    local <- data.frame(
      locus_id = locus$locus_id,
      snp_id = locus$lead_snp_id,
      chrom = locus$chrom,
      bp = locus$lead_bp,
      mlog10_p = -log10(pmax(locus$lead_p, .Machine$double.xmin)),
      r2 = 1,
      stringsAsFactors = FALSE
    )
  }

  local$c <- vapply(local$r2, r2_bin, character(1))
  local_plot <- sample_background_points(
    local,
    keep_mask = local$is_lead | local$r2 >= 0.2 | local$p_value <= 1e-4,
    max_background_points = 20000L,
    order_cols = c("bp", "mlog10_p")
  )
  ymax <- max(local$mlog10_p, -log10(pmax(locus$lead_p, .Machine$double.xmin)), na.rm = TRUE)
  ylim <- c(y_min, max(y_min + 0.1, ymax * 1.08))

  plot(c(locus$panel_start, locus$panel_end), ylim, type = "n", xlab = "", xaxt = "n", ylab = "-log10(P)", main = "")
  rect(locus$qtl_start, ylim[1], locus$qtl_end, ylim[2], col = alpha_col(locus$plot_color, 0.12), border = NA)
  points(local_plot$bp, local_plot$mlog10_p, pch = 16, cex = 0.75, col = alpha_col(local_plot$c, 0.9))
  points(
    locus$lead_bp,
    -log10(pmax(locus$lead_p, .Machine$double.xmin)),
    pch = 23,
    bg = "#111111",
    col = "white",
    cex = 1.1
  )
  title(
    main = sprintf(
      "%s | Chr%s:%s-%s | lead %s (%s)",
      locus$locus_id,
      locus$chrom,
      format(locus$qtl_start, big.mark = ",", scientific = FALSE, trim = TRUE),
      format(locus$qtl_end, big.mark = ",", scientific = FALSE, trim = TRUE),
      locus$lead_snp_id,
      locus$qtl_type
    ),
    adj = 0,
    cex.main = 0.8
  )
}

plot_gene_panel <- function(locus, gene_df, highlight_df) {
  axis_title <- sprintf("Chr%s position (bp)", locus$chrom)
  panel_genes <- gene_df[
    gene_df$chrom == locus$chrom &
      gene_df$start_bp <= locus$panel_end &
      gene_df$end_bp >= locus$panel_start,
    ,
    drop = FALSE
  ]

  if (nrow(panel_genes) == 0) {
    plot(
      c(locus$panel_start, locus$panel_end),
      c(-0.8, 1.2),
      type = "n",
      xlab = "",
      ylab = "",
      yaxt = "n",
      bty = "n"
    )
    rect(locus$qtl_start, -0.8, locus$qtl_end, 1.2, col = alpha_col(locus$plot_color, 0.12), border = NA)
    text(mean(c(locus$panel_start, locus$panel_end)), 0.2, labels = "No gene models provided in this interval", cex = 0.8, col = "#666666")
    mtext(axis_title, side = 1, line = 2.6, cex = 0.9)
    return(invisible(NULL))
  }

  lane_map <- layout_lanes(panel_genes)
  max_lane <- if (length(lane_map) == 0) 0 else max(lane_map)
  plot(
    c(locus$panel_start, locus$panel_end),
    c(-0.8, max_lane + 1.2),
    type = "n",
    xlab = "",
    ylab = "",
    yaxt = "n",
    bty = "n"
  )
  rect(locus$qtl_start, -0.8, locus$qtl_end, max_lane + 1.2, col = alpha_col(locus$plot_color, 0.12), border = NA)

  high_ids <- character()
  if (nrow(highlight_df) > 0) {
    high_ids <- highlight_df$gene_id[highlight_df$locus_id == locus$locus_id]
  }

  op <- par(xpd = NA)
  on.exit(par(op), add = TRUE)

  for (i in seq_len(nrow(panel_genes))) {
    gene <- panel_genes[i, , drop = FALSE]
    lane <- lane_map[[gene$gene_id]]
    y <- max_lane - lane
    is_highlight <- gene$gene_id %in% high_ids
    col <- if (is_highlight) locus$plot_color else "#555555"
    lwd <- if (is_highlight) 2.5 else 1.4
    segments(gene$start_bp, y, gene$end_bp, y, col = col, lwd = lwd, lend = "round")
    if (is_highlight) {
      label_row <- highlight_df[
        highlight_df$locus_id == locus$locus_id & highlight_df$gene_id == gene$gene_id,
        ,
        drop = FALSE
      ]
      label <- if (
        nrow(label_row) > 0 &&
        "label" %in% colnames(label_row) &&
        !is.na(label_row$label[[1]]) &&
        label_row$label[[1]] != ""
      ) {
        label_row$label[[1]]
      } else {
        gene$gene_id
      }
      anchor <- if (gene$strand == "+") gene$start_bp else gene$end_bp
      adj <- if (gene$strand == "+") c(0, 0) else c(1, 0)
      text(anchor, y + 0.2, labels = label, cex = 0.55, pos = 3, col = col, adj = adj)
    }
  }
  mtext(axis_title, side = 1, line = 2.6, cex = 0.9)
}

save_plot <- function(trait_label, global_df, loci, local_df, chrom_sizes, gene_df, highlight_df, out_path, y_min, genomewide_p) {
  if (nrow(loci) == 0) stop("No QTL regions available for plotting.")
  if (nrow(global_df) == 0) stop("global_manhattan is empty.")

  offsets_obj <- build_offsets(chrom_sizes)
  offsets <- offsets_obj$offsets
  centers <- offsets_obj$centers
  heights <- c(2.9, rep(c(2, 1), nrow(loci)))

  open_device <- if (grepl("\\.png$", out_path, ignore.case = TRUE)) {
    function() png(out_path, width = 1800, height = max(900, 420 + nrow(loci) * 520), res = 160)
  } else {
    function() pdf(out_path, width = 14, height = 4.6 + nrow(loci) * 2.7, paper = "special")
  }

  open_device()
  on.exit(dev.off(), add = TRUE)

  layout(matrix(seq_len(length(heights)), ncol = 1), heights = heights)
  par(mar = c(3.5, 4.2, 2.6, 1.2))
  plot_global_panel(trait_label, global_df, loci, offsets, centers, chrom_sizes, y_min, genomewide_p)

  for (i in seq_len(nrow(loci))) {
    locus <- loci[i, , drop = FALSE]
    par(mar = c(1.2, 4.2, 2.2, 1.2))
    plot_local_panel(locus, local_df = local_df, y_min = y_min)
    par(mar = c(4.8, 4.2, 0.4, 1.2))
    plot_gene_panel(locus, gene_df, highlight_df)
  }
}

main <- function() {
  opts <- parse_args()
  inputs <- normalize_inputs(opts)
  chrom_sizes <- inputs$chrom_sizes
  global_df <- inputs$global
  qtl_df <- inputs$qtl
  local_df <- inputs$local
  gene_df <- inputs$genes
  highlight_df <- inputs$highlights

  if (opts$max_loci > 0 && nrow(qtl_df) > opts$max_loci) {
    qtl_df <- qtl_df[seq_len(opts$max_loci), , drop = FALSE]
  }

  locus_colors <- c("#D73027", "#4575B4", "#1A9850", "#984EA3", "#FF7F00")
  qtl_df$plot_color <- vapply(
    seq_len(nrow(qtl_df)),
    function(i) locus_colors[[(i - 1L) %% length(locus_colors) + 1L]],
    character(1)
  )

  save_plot(
    trait_label = opts$trait_label,
    global_df = global_df,
    loci = qtl_df,
    local_df = local_df,
    chrom_sizes = chrom_sizes,
    gene_df = gene_df,
    highlight_df = highlight_df,
    out_path = opts$out,
    y_min = opts$y_min,
    genomewide_p = opts$genomewide_p
  )
}

main()
