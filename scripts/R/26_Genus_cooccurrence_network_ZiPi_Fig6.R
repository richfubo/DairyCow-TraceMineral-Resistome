##############################################################################
# Co-occurrence Network Analysis & Zi-Pi Topological Assessment
#
# Description:
#   1. Constructs phyloseq object (ps_filter) from raw feature tables.
#   2. Subsets samples for target time points (D30, D60, D90).
#   3. Calculates within-group Spearman co-occurrence networks (|r| >= 0.7, p < 0.05).
#   4. Computes topological Zi-Pi scores and generates no-label composite plots.
##############################################################################

rm(list = ls())
options(warn = -1)

# ============================================================================
# 0. Load Required Packages
# ============================================================================
packages <- c("tidyverse", "phyloseq", "igraph", "ggraph", "patchwork", "psych")

for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

# ============================================================================
# 1. Relative Paths Setup
# ============================================================================
data_dir <- "data"
out_dir  <- "output"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

metadata_file <- file.path(data_dir, "metadata.csv")
feature_file  <- file.path(data_dir, "feature-table.tsv")
taxonomy_file <- file.path(data_dir, "taxonomy.tsv")
ps_save_file  <- file.path(data_dir, "ps_filter.rds")

target_times   <- c("D30", "D60", "D90")
groups_to_plot <- c("ITM", "OTM1", "OTM2")

# Check required raw input files
required_files <- c(metadata_file, feature_file, taxonomy_file)
if (!all(file.exists(required_files))) {
  stop("Some required input files are missing in 'data/' directory. Please verify file paths.")
}

# ============================================================================
# 2. Step 1: Build & Save phyloseq Object (ps_filter)
# ============================================================================
cat("\n========================================\n")
cat("Step 1/2: Processing raw files & building ps_filter\n")
cat("========================================\n")

# 1) Load Metadata
meta <- read.csv(metadata_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

if (!"ID" %in% colnames(meta)) meta$ID <- rownames(meta)
meta$Time  <- as.character(meta$Time)
meta$Group <- as.character(meta$Group)
meta$Group <- dplyr::recode(meta$Group, "A" = "ITM", "B" = "OTM1", "C" = "OTM2")

# 2) Load Feature Table (handle BIOM headers if present)
first_line <- readLines(feature_file, n = 1, warn = FALSE)
skip_n <- ifelse(grepl("^# Constructed from biom file", first_line), 1, 0)

feature <- read.delim(
  feature_file, sep = "\t", header = TRUE, skip = skip_n,
  row.names = 1, check.names = FALSE, comment.char = "", quote = ""
)

# 3) Load Taxonomy Table
taxonomy <- read.delim(
  taxonomy_file, sep = "\t", header = TRUE,
  row.names = 1, check.names = FALSE, comment.char = "", quote = ""
)

# 4) Construct Phyloseq Object
ps_raw <- phyloseq(
  otu_table(as.matrix(feature), taxa_are_rows = TRUE),
  tax_table(as.matrix(taxonomy)),
  sample_data(meta)
)

# 5) Filter for Target Times (D30, D60, D90) & Groups, removing 0-abundance taxa
keep_samples <- rownames(meta)[meta$Time %in% target_times & meta$Group %in% groups_to_plot]
ps_sub <- prune_samples(keep_samples, ps_raw)
ps_filter <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)

# Save processed RDS
saveRDS(ps_filter, ps_save_file)
cat("✅ ps_filter.rds created and saved to 'data/ps_filter.rds'\n")

cat("\nRetained Sample Distribution (D30 + D60 + D90):\n")
print(table(sample_data(ps_filter)$Group, sample_data(ps_filter)$Time))

# ============================================================================
# 3. Core Network & Zi-Pi Computation Function
# ============================================================================
get_network_and_plots_clean <- function(ps_input, group_name, color_low, color_high, group_color) {
  cat(paste0("\n>>> Analyzing Group: ", group_name, "...\n"))
  
  meta_sub <- as(sample_data(ps_input), "data.frame")
  target_samples <- rownames(meta_sub)[meta_sub$Group == group_name]
  
  if (length(target_samples) < 3) return(NULL)
  
  sub_ps <- prune_samples(target_samples, ps_input)
  sub_ps <- prune_taxa(taxa_sums(sub_ps) > 0, sub_ps)
  
  otu_dat <- t(as(otu_table(sub_ps), "matrix"))
  
  # Calculate Spearman Correlation
  cor_res <- tryCatch(
    psych::corr.test(otu_dat, method = "spearman", adjust = "none", ci = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(cor_res)) return(NULL)
  
  r_mat <- cor_res$r
  p_mat <- cor_res$p
  
  # Filter by Thresholds (|r| >= 0.7, p < 0.05)
  r_mat[abs(r_mat) < 0.7] <- 0
  r_mat[p_mat >= 0.05]    <- 0
  diag(r_mat)             <- 0
  
  ig <- graph_from_adjacency_matrix(r_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
  ig <- delete.vertices(ig, which(degree(ig) == 0))
  
  if (vcount(ig) == 0) return(NULL)
  
  # --- Plot Co-occurrence Network ---
  E(ig)$weight_abs <- abs(E(ig)$weight)
  V(ig)$degree     <- degree(ig)
  
  p_net <- ggraph(ig, layout = "fr", weights = weight_abs) +
    geom_edge_link(color = "grey85", width = 0.4, alpha = 0.6) +
    geom_node_point(aes(fill = degree, size = degree), shape = 21, color = "white", stroke = 0.5) +
    scale_fill_gradient(low = color_low, high = color_high, guide = "none") +
    scale_size(range = c(2, 6), guide = "none") +
    theme_void() +
    labs(
      title    = group_name,
      subtitle = paste0("Nodes: ", vcount(ig), " | Edges: ", ecount(ig))
    ) +
    theme(
      plot.title    = element_text(face = "bold", size = 16, hjust = 0.5, color = group_color),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray50")
    )
  
  # --- Calculate Zi-Pi Topological Metrics ---
  ig_pos <- delete.edges(ig, E(ig)[weight < 0])
  
  if (ecount(ig_pos) > 0) {
    mod_res <- tryCatch(cluster_fast_greedy(ig_pos), error = function(e) NULL)
    module_membership <- if(!is.null(mod_res)) membership(mod_res) else rep(1, vcount(ig))
  } else {
    module_membership <- rep(1, vcount(ig))
  }
  
  names(module_membership) <- names(V(ig))
  res_df <- data.frame(Node = names(V(ig)), Zi = NA, Pi = NA)
  
  for (i in 1:nrow(res_df)) {
    node             <- res_df$Node[i]
    mod              <- module_membership[node]
    neighbors_nodes  <- neighbors(ig, node)
    ki               <- length(neighbors_nodes)
    neighbors_in_mod <- neighbors_nodes[names(neighbors_nodes) %in% names(which(module_membership == mod))]
    kis              <- length(neighbors_in_mod)
    mod_nodes        <- names(which(module_membership == mod))
    k_s              <- degree(ig, v = mod_nodes)
    mean_ks          <- mean(k_s)
    sd_ks            <- sd(k_s)
    
    if (is.na(sd_ks) || sd_ks == 0) res_df$Zi[i] <- 0 else res_df$Zi[i] <- (kis - mean_ks) / sd_ks
    
    sum_term <- 0
    for (m in unique(module_membership)) {
      kis_m <- length(neighbors_nodes[names(neighbors_nodes) %in% names(which(module_membership == m))])
      if (ki > 0) sum_term <- sum_term + (kis_m / ki)^2
    }
    res_df$Pi[i] <- 1 - sum_term
  }
  
  res_df$Role <- case_when(
    res_df$Zi > 2.5  & res_df$Pi > 0.62  ~ "Network Hubs",
    res_df$Zi > 2.5  & res_df$Pi <= 0.62 ~ "Module Hubs",
    res_df$Zi <= 2.5 & res_df$Pi > 0.62  ~ "Connectors",
    TRUE                                  ~ "Peripherals"
  )
  
  p_zipi <- ggplot(res_df, aes(x = Pi, y = Zi)) +
    annotate("rect", xmin = 0,    xmax = 0.62, ymin = -Inf, ymax = 2.5, fill = "gray95", alpha = 0.5) +
    annotate("rect", xmin = 0,    xmax = 0.62, ymin = 2.5,  ymax = Inf, fill = "#FFF3E0", alpha = 0.4) +
    annotate("rect", xmin = 0.62, xmax = 1,    ymin = -Inf, ymax = 2.5, fill = "#E1F5FE", alpha = 0.4) +
    annotate("rect", xmin = 0.62, xmax = 1,    ymin = 2.5,  ymax = Inf, fill = "#E8F5E9", alpha = 0.4) +
    geom_point(aes(color = Role == "Peripherals"), size = 2.5, alpha = 0.7) +
    scale_color_manual(values = c("TRUE" = "gray70", "FALSE" = group_color), guide = "none") +
    geom_vline(xintercept = 0.62, linetype = "dashed", color = "gray50", linewidth = 0.3) +
    geom_hline(yintercept = 2.5,  linetype = "dashed", color = "gray50", linewidth = 0.3) +
    labs(x = "Among-module connectivity (Pi)", y = "Within-module connectivity (Zi)") +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank())
  
  return(list(net = p_net, zipi = p_zipi))
}

# ============================================================================
# 4. Step 2: Batch Analysis & Composite Plotting
# ============================================================================
cat("\n========================================\n")
cat("Step 2/2: Computing networks and Zi-Pi plots\n")
cat("========================================\n")

res_ITM  <- get_network_and_plots_clean(ps_filter, "ITM",  "#E0F2F1", "#00695C", "#8ECFC9")
res_OTM1 <- get_network_and_plots_clean(ps_filter, "OTM1", "#FFF3E0", "#E65100", "#FFBE7A")
res_OTM2 <- get_network_and_plots_clean(ps_filter, "OTM2", "#FFEBEE", "#C62828", "#FA7F6F")

if (!is.null(res_ITM) & !is.null(res_OTM1) & !is.null(res_OTM2)) {
  row1 <- res_ITM$net + res_OTM1$net + res_OTM2$net
  row2 <- res_ITM$zipi + res_OTM1$zipi + res_OTM2$zipi
  
  final_plot <- row1 / row2 + plot_annotation(tag_levels = 'A')
  
  out_pdf <- file.path(out_dir, "Figure6_Clean_Composite.pdf")
  out_png <- file.path(out_dir, "Figure6_Clean_Composite.png")
  
  ggsave(out_pdf, final_plot, width = 16, height = 10)
  ggsave(out_png, final_plot, width = 16, height = 10, dpi = 600)
  
  cat("\n========================================\n")
  cat("Analysis successfully completed!\n")
  cat("Outputs generated:\n")
  cat("1. RDS object : ", ps_save_file, "\n")
  cat("2. PDF Figure : ", out_pdf, "\n")
  cat("3. PNG Figure : ", out_png, "\n")
  cat("========================================\n")
}