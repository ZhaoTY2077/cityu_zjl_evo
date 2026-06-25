###############################################################################
#  GO (BP, MF, CC) & KEGG Enrichment Analysis
#  for 391 human RBP genes with 100% human-chimp similarity
#  Uses org.Hs.eg.db, GO.db, KEGGREST, ggplot2 (same pipeline as previous)
###############################################################################

library(org.Hs.eg.db)
library(GO.db)
library(AnnotationDbi)
library(KEGGREST)
library(ggplot2)
library(reshape2)
library(parallel)

# ---- 1. Read gene list ------------------------------------------------------
genes <- readLines("/saturn/zhaoty/evo_project/files/chimp_100pct_genes.txt")
genes <- unique(trimws(genes))
genes <- genes[genes != "" & !grepl("^#", genes)]
message("Input genes: ", length(genes))

# Map symbols to Entrez IDs
entrez_map <- tryCatch(
  select(org.Hs.eg.db, keys = genes, keytype = "SYMBOL",
         columns = c("ENTREZID", "GENENAME")),
  error = function(e) NULL
)

entrez_map <- entrez_map[!is.na(entrez_map$ENTREZID) & entrez_map$ENTREZID != "", ]
entrez <- unique(entrez_map$ENTREZID)
message("Mapped to ", length(entrez), " Entrez IDs")

# ---- 2. Universe ------------------------------------------------------------
universe <- keys(org.Hs.eg.db, keytype = "ENTREZID")
message("Universe size: ", length(universe))

# ---- 3. GO enrichment helper (same as before) -------------------------------
do_go_enrichment <- function(gene_vec, universe_vec, ontology = "BP") {
  # Get all GO IDs for this ontology
  all_go_info <- select(GO.db, keys = ls(GOTERM), columns = c("GOID", "ONTOLOGY", "TERM"))
  go_ids_ont <- all_go_info$GOID[all_go_info$ONTOLOGY == ontology]
  message("  Ontology ", ontology, " — ", length(go_ids_ont), " terms to test")

  # Gene-to-GO mappings
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
  k        <- length(gene_set)

  .enrich_one <- function(goid) {
    term_genes <- go2gene[[goid]]
    if (is.null(term_genes)) return(NULL)
    term_genes <- intersect(term_genes, universe_vec)
    M <- length(term_genes)
    if (M < 2) return(NULL)

    x <- length(intersect(gene_set, term_genes))
    if (x < 1) return(NULL)

    pval <- phyper(x - 1, M, n_univ - M, k, lower.tail = FALSE)

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

  res$p.adjust <- p.adjust(res$pvalue, method = "BH")
  res <- res[order(res$p.adjust), ]
  rownames(res) <- NULL
  res
}

# ---- 4. Run GO enrichment ---------------------------------------------------
message("\n=== GO Enrichment Analysis ===")
go_bp <- do_go_enrichment(entrez, universe, "BP")
go_mf <- do_go_enrichment(entrez, universe, "MF")
go_cc <- do_go_enrichment(entrez, universe, "CC")

sig_bp <- if (!is.null(go_bp)) go_bp[go_bp$p.adjust < 0.05, ] else NULL
sig_mf <- if (!is.null(go_mf)) go_mf[go_mf$p.adjust < 0.05, ] else NULL
sig_cc <- if (!is.null(go_cc)) go_cc[go_cc$p.adjust < 0.05, ] else NULL
message("Sig BP: ", nrow(sig_bp), "  Sig MF: ", nrow(sig_mf), "  Sig CC: ", nrow(sig_cc))

# ---- 5. KEGG enrichment -----------------------------------------------------
message("\n=== KEGG Enrichment Analysis ===")
hsa_pathways <- tryCatch(keggList("pathway", "hsa"), error = function(e) NULL)

if (!is.null(hsa_pathways)) {
  path_ids <- gsub("^path:", "", names(hsa_pathways))
  message("  KEGG pathways: ", length(path_ids))

  kegg_res_list <- lapply(seq_along(path_ids), function(i) {
    pid <- path_ids[i]
    Sys.sleep(0.1)
    pw <- tryCatch(keggGet(pid), error = function(e) NULL)
    if (is.null(pw) || is.null(pw[[1]]$GENE)) return(NULL)

    pw_genes <- pw[[1]]$GENE
    pw_entrez <- pw_genes[seq(1, length(pw_genes), 2)]
    pw_entrez <- intersect(pw_entrez, universe)

    M <- length(pw_entrez)
    if (M < 2) return(NULL)

    x <- length(intersect(entrez, pw_entrez))
    if (x < 1) return(NULL)

    pval <- phyper(x - 1, M, length(universe) - M, length(entrez), lower.tail = FALSE)

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
  kegg_res <- NULL; sig_kegg <- NULL
  message("  Could not connect to KEGG REST server")
}

# ---- 6. Save full results ---------------------------------------------------
outdir <- "/saturn/zhaoty/evo_project/enrichment_chimp_results"
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

# dot-plot
dot_plot <- function(df, title, n_top = 20) {
  if (is.null(df) || nrow(df) == 0) { message("  No sig terms: ", title); return(NULL) }
  df <- head(df, n_top)
  gr <- strsplit(df$GeneRatio, "/")
  df$Ratio <- sapply(gr, function(x) as.numeric(x[1]) / as.numeric(x[2]))
  df$Description <- factor(df$Description, levels = rev(df$Description))
  ggplot(df, aes(x = Ratio, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "#E41A1C", high = "#377EB8",
                         name = "p.adjust", trans = "log10") +
    scale_size_continuous(range = c(3, 10)) +
    labs(title = title, x = "GeneRatio", y = "") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "right", panel.grid.minor = element_blank())
}

# bar-plot
bar_plot <- function(df, title, n_top = 20) {
  if (is.null(df) || nrow(df) == 0) { message("  No sig terms: ", title); return(NULL) }
  df <- head(df, n_top)
  df$Description <- factor(df$Description, levels = rev(df$Description))
  ggplot(df, aes(x = -log10(p.adjust), y = Description)) +
    geom_bar(stat = "identity", aes(fill = Count)) +
    scale_fill_gradient(low = "#92C5DE", high = "#B2182B", name = "Count") +
    labs(title = title, x = expression(-log[10](p.adjust)), y = "") +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "right", panel.grid.minor = element_blank())
}

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

# Combined multi-panel dot-plot
combined <- rbind(
  if (!is.null(sig_bp)) cbind(head(sig_bp, 10), Ont = "BP") else NULL,
  if (!is.null(sig_mf)) cbind(head(sig_mf, 10), Ont = "MF") else NULL,
  if (!is.null(sig_cc)) cbind(head(sig_cc, 10), Ont = "CC") else NULL
)
if (!is.null(combined) && nrow(combined) > 0) {
  gr <- strsplit(combined$GeneRatio, "/")
  combined$Ratio <- sapply(gr, function(x) as.numeric(x[1]) / as.numeric(x[2]))
  combined$Description <- factor(combined$Description, levels = rev(unique(combined$Description)))
  p <- ggplot(combined, aes(x = Ratio, y = Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "#E41A1C", high = "#377EB8",
                         name = "p.adjust", trans = "log10") +
    scale_size_continuous(range = c(3, 10)) +
    facet_grid(Ont ~ ., scales = "free_y", space = "free_y") +
    labs(title = "GO Enrichment (BP / MF / CC)", x = "GeneRatio", y = "") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          strip.text.y = element_text(angle = 0, face = "bold", size = 11),
          legend.position = "right", panel.grid.minor = element_blank())
  ggsave(file.path(outdir, "GO_combined_dotplot.pdf"), p, width = 12, height = 10)
}

# Combined bar-plot
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
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          strip.text.y = element_text(angle = 0, face = "bold", size = 11),
          legend.position = "right", panel.grid.minor = element_blank())
  ggsave(file.path(outdir, "GO_combined_barplot.pdf"), p, width = 12, height = 10)
}

# ---- 8. Summary CSV ---------------------------------------------------------
summary_list <- list()
if (!is.null(sig_bp) && nrow(sig_bp) > 0)
  summary_list[["GO_BP"]] <- data.frame(Category = "GO BP",
    Description = head(sig_bp,10)$Description,
    p.adjust = signif(head(sig_bp,10)$p.adjust, 3),
    Count = head(sig_bp,10)$Count, Genes = head(sig_bp,10)$geneID)
if (!is.null(sig_mf) && nrow(sig_mf) > 0)
  summary_list[["GO_MF"]] <- data.frame(Category = "GO MF",
    Description = head(sig_mf,10)$Description,
    p.adjust = signif(head(sig_mf,10)$p.adjust, 3),
    Count = head(sig_mf,10)$Count, Genes = head(sig_mf,10)$geneID)
if (!is.null(sig_cc) && nrow(sig_cc) > 0)
  summary_list[["GO_CC"]] <- data.frame(Category = "GO CC",
    Description = head(sig_cc,10)$Description,
    p.adjust = signif(head(sig_cc,10)$p.adjust, 3),
    Count = head(sig_cc,10)$Count, Genes = head(sig_cc,10)$geneID)
if (!is.null(sig_kegg) && nrow(sig_kegg) > 0)
  summary_list[["KEGG"]] <- data.frame(Category = "KEGG",
    Description = head(sig_kegg,10)$Description,
    p.adjust = signif(head(sig_kegg,10)$p.adjust, 3),
    Count = head(sig_kegg,10)$Count, Genes = head(sig_kegg,10)$geneID)

if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, summary_list)
  write.csv(summary_df, file.path(outdir, "enrichment_summary.csv"), row.names = FALSE)
}

message("\n=== DONE ===")
message("All results saved to: ", outdir)
