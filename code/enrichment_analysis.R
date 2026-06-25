###############################################################################
#  GO (BP, MF, CC) & KEGG Enrichment Analysis for 51 Human Genes
#  without clusterProfiler — uses org.Hs.eg.db, GO.db, KEGGREST, ggplot2
###############################################################################

library(org.Hs.eg.db)
library(GO.db)
library(AnnotationDbi)
library(KEGGREST)
library(ggplot2)
library(DOSE)          # for geneID symbol conversion helper
library(reshape2)
library(parallel)

# ---- 1. Gene list -----------------------------------------------------------
genes <- c(
  "POLR2J", "RPS20", "RPL3", "PTGES3L-AARSD1", "NOP2", "SNRNP27",
  "NXF5", "RPL36", "UTP3", "HNRNPA1L2", "RPS8", "FIP1L1", "RPL10",
  "MRPS28", "POLR2K", "FRG1BP", "RPL8", "EIF4A1", "FUBP1", "POLR2G",
  "RPSA", "POLR2J3", "LCMT2", "DDX19A", "RNASE3", "HNRNPA3", "NCBP2L",
  "SMN1", "RBM4", "RTL3", "EIF2S3B", "CMTR2", "RNASE10", "TRMT12",
  "RPS27L", "PURA", "ZCCHC13", "PNRC2", "RPL37A", "NUDT16", "FAM120A2P",
  "RNASE13", "ZRSR2P1", "ADAT3", "RPS29", "FASTKD5", "POLR2J2", "RPL41",
  "MCTS1", "RPL17", "FMC1-LUC7L2"
)

message("Input genes: ", length(genes))

# Map symbols to Entrez IDs
entrez_map <- tryCatch(
  select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL",
         columns = c("ENTREZID", "GENENAME")),
  error = function(e) NULL
)

# Drop genes that didn't map
entrez_map <- entrez_map[!is.na(entrez_map$ENTREZID) & entrez_map$ENTREZID != "", ]
entrez <- unique(entrez_map$ENTREZID)
message("Mapped to ", length(entrez), " Entrez IDs")

# ---- 2. Universe (all human protein-coding genes with Entrez) ----------------
universe <- keys(org.Hs.eg.db, keytype = "ENTREZID")
message("Universe size: ", length(universe))

# ---- 3. Helper: hypergeometric enrichment for one GO ontology ---------------
#     Returns a data.frame with GO ID, description, ratio, p-value, adjusted p,
#     and the overlapping gene symbols.
do_go_enrichment <- function(gene_vec, universe_vec, ontology = "BP") {

  # Get all GO IDs for this ontology via GO.db's select
  all_go_info <- select(GO.db, keys = ls(GOTERM), columns = c("GOID", "ONTOLOGY", "TERM"))
  go_ids_ont <- all_go_info$GOID[all_go_info$ONTOLOGY == ontology]
  message("  Ontology ", ontology, " — ", length(go_ids_ont), " terms to test")

  # Get gene-to-GO mappings from org.Hs.eg.db
  go2gene <- tryCatch(
    as.list(org.Hs.egGO2ALLEGS),
    error = function(e) NULL
  )
  if (is.null(go2gene)) {
    message("  Could not fetch GO2ALLEGS")
    return(NULL)
  }

  gene_set <- unique(gene_vec)
  n_univ   <- length(unique(universe_vec))
  k        <- length(gene_set)   # size of query in universe

  .enrich_one <- function(goid) {
    # genes annotated to this GO term (must be in universe)
    term_genes <- go2gene[[goid]]
    if (is.null(term_genes)) return(NULL)
    term_genes <- intersect(term_genes, universe_vec)
    M <- length(term_genes)   # term size in universe
    if (M < 2) return(NULL)

    x <- length(intersect(gene_set, term_genes))   # overlap
    if (x < 1) return(NULL)

    # hypergeometric / Fisher's exact test
    pval <- phyper(x - 1, M, n_univ - M, k, lower.tail = FALSE)

    # genes contributing
    overlap_genes <- intersect(gene_set, term_genes)
    overlap_symbols <- entrez_map$SYMBOL[entrez_map$ENTREZID %in% overlap_genes]
    overlap_symbols <- unique(overlap_symbols)

    data.frame(
      GOID        = goid,
      Ontology    = ontology,
      Description = Term(GOTERM[[goid]]),
      GeneRatio   = paste0(x, "/", k),
      BgRatio     = paste0(M, "/", n_univ),
      pvalue      = pval,
      Count       = x,
      geneID      = paste(sort(overlap_symbols), collapse = "/"),
      stringsAsFactors = FALSE
    )
  }

  res_list <- mclapply(go_ids_ont, .enrich_one, mc.cores = 4)
  res <- do.call(rbind, res_list[!sapply(res_list, is.null)])
  if (is.null(res) || nrow(res) == 0) return(NULL)

  # BH correction
  res$p.adjust <- p.adjust(res$pvalue, method = "BH")
  res <- res[order(res$p.adjust), ]
  rownames(res) <- NULL
  res
}

# ---- 4. Run GO enrichment for BP, MF, CC ------------------------------------
message("\n=== GO Enrichment Analysis ===")
go_bp <- do_go_enrichment(entrez, universe, "BP")
go_mf <- do_go_enrichment(entrez, universe, "MF")
go_cc <- do_go_enrichment(entrez, universe, "CC")

# Filter significant (p.adjust < 0.05)
sig_bp <- if (!is.null(go_bp)) go_bp[go_bp$p.adjust < 0.05, ] else NULL
sig_mf <- if (!is.null(go_mf)) go_mf[go_mf$p.adjust < 0.05, ] else NULL
sig_cc <- if (!is.null(go_cc)) go_cc[go_cc$p.adjust < 0.05, ] else NULL

message("Sig BP: ", nrow(sig_bp), "  Sig MF: ", nrow(sig_mf), "  Sig CC: ", nrow(sig_cc))

# ---- 5. KEGG enrichment via KEGGREST ----------------------------------------
message("\n=== KEGG Enrichment Analysis ===")

# Get all human KEGG pathways
hsa_pathways <- tryCatch(
  keggList("pathway", "hsa"),
  error = function(e) NULL
)

if (!is.null(hsa_pathways)) {
  path_ids <- gsub("^path:", "", names(hsa_pathways))
  message("  KEGG pathways: ", length(path_ids))

  # For each pathway, get its gene list and test overlap
  kegg_res_list <- lapply(seq_along(path_ids), function(i) {
    pid <- path_ids[i]
    Sys.sleep(0.1)   # be gentle to KEGG server
    pw <- tryCatch(keggGet(pid), error = function(e) NULL)
    if (is.null(pw) || is.null(pw[[1]]$GENE)) return(NULL)

    pw_genes <- pw[[1]]$GENE
    # KEGG returns lines: EntrezID, name, ... — take every other line starting from 1
    pw_entrez <- pw_genes[seq(1, length(pw_genes), 2)]
    pw_entrez <- intersect(pw_entrez, universe)

    M <- length(pw_entrez)
    if (M < 2) return(NULL)

    x <- length(intersect(entrez, pw_entrez))
    if (x < 1) return(NULL)

    pval <- phyper(x - 1, M, length(universe) - M, length(entrez),
                   lower.tail = FALSE)

    overlap_symbols <- entrez_map$SYMBOL[entrez_map$ENTREZID %in%
                                          intersect(entrez, pw_entrez)]
    overlap_symbols <- unique(overlap_symbols)

    data.frame(
      ID          = pid,
      Description = hsa_pathways[i],
      GeneRatio   = paste0(x, "/", length(entrez)),
      BgRatio     = paste0(M, "/", length(universe)),
      pvalue      = pval,
      Count       = x,
      geneID      = paste(sort(overlap_symbols), collapse = "/"),
      stringsAsFactors = FALSE
    )
  })
  kegg_res <- do.call(rbind, kegg_res_list[!sapply(kegg_res_list, is.null)])
  if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
    kegg_res$p.adjust <- p.adjust(kegg_res$pvalue, method = "BH")
    kegg_res <- kegg_res[order(kegg_res$p.adjust), ]
    rownames(kegg_res) <- NULL
  } else {
    kegg_res <- NULL
    message("  No KEGG enrichment results")
  }
  sig_kegg <- if (!is.null(kegg_res)) kegg_res[kegg_res$p.adjust < 0.05, ] else NULL
  message("Sig KEGG: ", nrow(sig_kegg))
} else {
  kegg_res <- NULL
  sig_kegg <- NULL
  message("  Could not connect to KEGG REST server")
}

# ---- 6. Save full results ---------------------------------------------------
outdir <- "/saturn/zhaoty/evo_project/results/enrichment_analysis"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

save(go_bp, go_mf, go_cc, sig_bp, sig_mf, sig_cc,
     kegg_res, sig_kegg, entrez_map,
     file = file.path(outdir, "enrichment_results.RData"))

write.csv(if (!is.null(go_bp)) go_bp else data.frame(),
          file.path(outdir, "GOBP_full.csv"), row.names = FALSE)
write.csv(if (!is.null(go_mf)) go_mf else data.frame(),
          file.path(outdir, "GOMF_full.csv"), row.names = FALSE)
write.csv(if (!is.null(go_cc)) go_cc else data.frame(),
          file.path(outdir, "GOCC_full.csv"), row.names = FALSE)
write.csv(if (!is.null(kegg_res)) kegg_res else data.frame(),
          file.path(outdir, "KEGG_full.csv"), row.names = FALSE)

# ---- 7. Visualization -------------------------------------------------------
message("\n=== Generating Plots ===")

# Color palette
my_colors <- c(BP = "#E64B35", MF = "#4DBBD5", CC = "#00A087", KEGG = "#3C5488")

# ---------- helper: enrichment dot-plot ----------
dot_plot <- function(df, title, n_top = 20) {
  if (is.null(df) || nrow(df) == 0) {
    message("  No significant terms for ", title)
    return(NULL)
  }
  df <- head(df, n_top)
  # parse GeneRatio
  gr <- strsplit(df$GeneRatio, "/")
  df$Ratio <- sapply(gr, function(x) as.numeric(x[1]) / as.numeric(x[2]))

  df$Description <- factor(df$Description, levels = rev(df$Description))

  p <- ggplot(df, aes(x = Ratio, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "#E41A1C", high = "#377EB8",
                         name = "p.adjust", trans = "log10") +
    scale_size_continuous(range = c(3, 10)) +
    labs(title = title, x = "GeneRatio", y = "") +
    theme_minimal(base_size = 13) +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "right",
          panel.grid.minor = element_blank())
  p
}

# ---------- helper: bar-plot ----------
bar_plot <- function(df, title, n_top = 20) {
  if (is.null(df) || nrow(df) == 0) {
    message("  No significant terms for ", title)
    return(NULL)
  }
  df <- head(df, n_top)
  df$Description <- factor(df$Description, levels = rev(df$Description))

  p <- ggplot(df, aes(x = -log10(p.adjust), y = Description)) +
    geom_bar(stat = "identity", aes(fill = Count)) +
    scale_fill_gradient(low = "#92C5DE", high = "#B2182B", name = "Count") +
    labs(title = title, x = expression(-log[10](p.adjust)), y = "") +
    theme_minimal(base_size = 13) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "right",
          panel.grid.minor = element_blank())
  p
}

# ---------- helper: simplified "cnet" style ----------
# For each term, plot its overlapping genes as a bipartite graph
#
# NOTE: This is a much-simplified version of cnetplot. We generate a
# bipartite layout using igraph if available, or a semi-random layout with base R.
plot_go_network <- function(df, title, max_terms = 10, seed = 42) {
  if (is.null(df) || nrow(df) == 0) {
    message("  No data for network plot: ", title)
    return(invisible(NULL))
  }

  # Check we have a reasonable number of rows
  top_df <- head(df, max_terms)
  if (nrow(top_df) == 0) return(invisible(NULL))

  # Parse geneID column
  term_list <- strsplit(top_df$geneID, "/")
  names(term_list) <- top_df$Description

  # Shorten term descriptions for display
  short_names <- function(nms, maxw = 40) {
    sapply(nms, function(x) {
      if (nchar(x) > maxw) paste0(substr(x, 1, maxw - 3), "...")
      else x
    })
  }

  # Build edge list
  edges <- do.call(rbind, lapply(names(term_list), function(term) {
    genes <- term_list[[term]]
    if (length(genes) == 0) return(NULL)
    data.frame(term = term, gene = genes, stringsAsFactors = FALSE)
  }))
  if (is.null(edges) || nrow(edges) == 0) return(invisible(NULL))

  # Unique terms and genes
  all_terms <- unique(edges$term)
  all_genes <- unique(edges$gene)

  # Simple radial layout: terms on left arc, genes on right arc
  n_terms <- length(all_terms)
  n_genes <- length(all_genes)

  # Positions
  y_t <- seq(1, 0, length.out = n_terms + 2)[-c(1, n_terms + 2)]
  y_g <- seq(1, 0, length.out = n_genes + 2)[-c(1, n_genes + 2)]

  layout_df <- data.frame(
    node   = c(all_terms, all_genes),
    type   = c(rep("Term", n_terms), rep("Gene", n_genes)),
    x      = c(rep(0, n_terms), rep(1, n_genes)),
    y      = c(y_t, y_g),
    stringsAsFactors = FALSE
  )

  # Merge in colors
  term_colors <- scales::hue_pal()(n_terms)
  names(term_colors) <- all_terms
  layout_df$color <- ifelse(layout_df$type == "Term",
                            term_colors[layout_df$node],
                            "#AAAAAA")
  layout_df$label <- ifelse(layout_df$type == "Term",
                            short_names(layout_df$node, 35),
                            layout_df$node)

  # Write to PDF as a grid-based network
  pdf(file.path(outdir, paste0(gsub("[ /]", "_", title), "_network.pdf")),
      width = 14, height = max(6, n_terms * 0.5 + n_genes * 0.2))
  par(mar = c(2, 2, 3, 2))
  plot.new()
  plot.window(xlim = c(-0.15, 1.15), ylim = c(-0.05, 1.05))

  # Edges (semi-transparent)
  for (i in seq_len(nrow(edges))) {
    from <- edges$term[i]
    to   <- edges$gene[i]
    x0 <- layout_df$x[layout_df$node == from]
    y0 <- layout_df$y[layout_df$node == from]
    x1 <- layout_df$x[layout_df$node == to]
    y1 <- layout_df$y[layout_df$node == to]
    segments(x0, y0, x1, y1, col = "#00000015", lwd = 1.2)
  }

  # Nodes
  points(layout_df$x, layout_df$y,
         pch = ifelse(layout_df$type == "Term", 21, 19),
         cex = ifelse(layout_df$type == "Term", 2.5, 1.8),
         col = "#333333",
         bg  = layout_df$color)

  # Labels
  text(layout_df$x - 0.01, layout_df$y,
       labels = layout_df$label,
       pos = ifelse(layout_df$x < 0.5, 2, 4),
       cex = 0.7, xpd = TRUE)
  title(title, cex.main = 1.2, font.main = 2)
  dev.off()
  message("  Network plot saved: ", title)
}

# ---------- Generate all plots ----------

# Dot plots
p <- dot_plot(sig_bp, "GO Biological Process (BP)")
if (!is.null(p)) ggsave(file.path(outdir, "GOBP_dotplot.pdf"), p, width = 10, height = 7)

p <- dot_plot(sig_mf, "GO Molecular Function (MF)")
if (!is.null(p)) ggsave(file.path(outdir, "GOMF_dotplot.pdf"), p, width = 10, height = 7)

p <- dot_plot(sig_cc, "GO Cellular Component (CC)")
if (!is.null(p)) ggsave(file.path(outdir, "GOCC_dotplot.pdf"), p, width = 10, height = 7)

p <- dot_plot(sig_kegg, "KEGG Pathway")
if (!is.null(p)) ggsave(file.path(outdir, "KEGG_dotplot.pdf"), p, width = 10, height = 7)

# Bar plots
p <- bar_plot(sig_bp, "GO Biological Process (BP)")
if (!is.null(p)) ggsave(file.path(outdir, "GOBP_barplot.pdf"), p, width = 10, height = 7)

p <- bar_plot(sig_mf, "GO Molecular Function (MF)")
if (!is.null(p)) ggsave(file.path(outdir, "GOMF_barplot.pdf"), p, width = 10, height = 7)

p <- bar_plot(sig_cc, "GO Cellular Component (CC)")
if (!is.null(p)) ggsave(file.path(outdir, "GOCC_barplot.pdf"), p, width = 10, height = 7)

p <- bar_plot(sig_kegg, "KEGG Pathway")
if (!is.null(p)) ggsave(file.path(outdir, "KEGG_barplot.pdf"), p, width = 10, height = 7)

# Combined multi-panel dot-plot (top 10 from each)
combined <- rbind(
  if (!is.null(sig_bp)) cbind(head(sig_bp, 10), Ont = "BP") else NULL,
  if (!is.null(sig_mf)) cbind(head(sig_mf, 10), Ont = "MF") else NULL,
  if (!is.null(sig_cc)) cbind(head(sig_cc, 10), Ont = "CC") else NULL
)

if (!is.null(combined) && nrow(combined) > 0) {
  gr <- strsplit(combined$GeneRatio, "/")
  combined$Ratio <- sapply(gr, function(x) as.numeric(x[1]) / as.numeric(x[2]))
  combined$Description <- factor(combined$Description,
                                  levels = rev(unique(combined$Description)))

  p <- ggplot(combined, aes(x = Ratio, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "#E41A1C", high = "#377EB8",
                         name = "p.adjust", trans = "log10") +
    scale_size_continuous(range = c(3, 10)) +
    facet_grid(Ont ~ ., scales = "free_y", space = "free_y") +
    labs(title = "GO Enrichment (BP / MF / CC)", x = "GeneRatio", y = "") +
    theme_minimal(base_size = 12) +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
          strip.text.y  = element_text(angle = 0, face = "bold", size = 11),
          legend.position = "right",
          panel.grid.minor = element_blank())
  ggsave(file.path(outdir, "GO_combined_dotplot.pdf"), p, width = 12, height = 10)
}

# Combined multi-panel bar-plot
combined_bar <- rbind(
  if (!is.null(sig_bp)) cbind(head(sig_bp, 10), Ont = "BP") else NULL,
  if (!is.null(sig_mf)) cbind(head(sig_mf, 10), Ont = "MF") else NULL,
  if (!is.null(sig_cc)) cbind(head(sig_cc, 10), Ont = "CC") else NULL
)

if (!is.null(combined_bar) && nrow(combined_bar) > 0) {
  combined_bar$Description <- factor(combined_bar$Description,
                                      levels = rev(unique(combined_bar$Description)))
  p <- ggplot(combined_bar, aes(x = -log10(p.adjust), y = Description)) +
    geom_bar(stat = "identity", aes(fill = Count)) +
    scale_fill_gradient(low = "#92C5DE", high = "#B2182B", name = "Count") +
    facet_grid(Ont ~ ., scales = "free_y", space = "free_y") +
    labs(title = "GO Enrichment (BP / MF / CC)", x = expression(-log[10](p.adjust)), y = "") +
    theme_minimal(base_size = 12) +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold", size = 14),
          strip.text.y  = element_text(angle = 0, face = "bold", size = 11),
          legend.position = "right",
          panel.grid.minor = element_blank())
  ggsave(file.path(outdir, "GO_combined_barplot.pdf"), p, width = 12, height = 10)
}

# Network-style plots (simplified cnet)
plot_go_network(sig_bp, "BP Gene-Term Network", max_terms = 10)
plot_go_network(sig_mf, "MF Gene-Term Network", max_terms = 10)
plot_go_network(sig_cc, "CC Gene-Term Network", max_terms = 10)
if (!is.null(sig_kegg)) plot_go_network(sig_kegg, "KEGG Gene-Pathway Network", max_terms = 10)

# ---- 8. Summary CSV ---------------------------------------------------------
summary_list <- list()
if (!is.null(sig_bp) && nrow(sig_bp) > 0) {
  tmp <- head(sig_bp, 10)
  summary_list[["GO_BP"]] <- data.frame(
    Category = "GO BP", Description = tmp$Description,
    p.adjust = signif(tmp$p.adjust, 3), Count = tmp$Count,
    Genes = tmp$geneID
  )
}
if (!is.null(sig_mf) && nrow(sig_mf) > 0) {
  tmp <- head(sig_mf, 10)
  summary_list[["GO_MF"]] <- data.frame(
    Category = "GO MF", Description = tmp$Description,
    p.adjust = signif(tmp$p.adjust, 3), Count = tmp$Count,
    Genes = tmp$geneID
  )
}
if (!is.null(sig_cc) && nrow(sig_cc) > 0) {
  tmp <- head(sig_cc, 10)
  summary_list[["GO_CC"]] <- data.frame(
    Category = "GO CC", Description = tmp$Description,
    p.adjust = signif(tmp$p.adjust, 3), Count = tmp$Count,
    Genes = tmp$geneID
  )
}
if (!is.null(sig_kegg) && nrow(sig_kegg) > 0) {
  tmp <- head(sig_kegg, 10)
  summary_list[["KEGG"]] <- data.frame(
    Category = "KEGG", Description = tmp$Description,
    p.adjust = signif(tmp$p.adjust, 3), Count = tmp$Count,
    Genes = tmp$geneID
  )
}
if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, summary_list)
  write.csv(summary_df, file.path(outdir, "enrichment_summary.csv"), row.names = FALSE)
}

message("\n=== DONE ===")
message("All results saved to: ", outdir)
