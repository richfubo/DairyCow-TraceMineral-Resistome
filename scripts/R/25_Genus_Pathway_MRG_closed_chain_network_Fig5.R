############################################################
# Fig. 5D Network Plot Analysis
# Genus – pathway category – MRG subtype closed-chain network
#
# D30 + D60 + D90
#
# Spearman correlation
# FDR < 0.05
# |rho| >= 0.5
#
# Nodes:
#    Genus = orange (#F4A62A)
#    Pathway category = green (#2CA25F)
#    MRG subtype = blue (#2B8CBE)
#
# Edges:
#    Positive = red (#D6604D)
#    Negative = blue (#4393C3)
#    Genus–Pathway / Pathway–MRG = solid
#    Direct Genus–MRG = dashed
############################################################

rm(list = ls())
options(warn = -1)


## =========================================================
## 0. Packages
## =========================================================

packages <- c(
  "tidyverse",
  "igraph",
  "ggrepel",
  "scales"
)

for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}


## =========================================================
## 1. File paths (Relative Paths)
## =========================================================

# 数据输入目录与输出目录设置
data_dir <- "data"
out_dir  <- "output"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

metadata_file     <- file.path(data_dir, "metadata.csv")
feature_file      <- file.path(data_dir, "feature-table.tsv")
taxonomy_file     <- file.path(data_dir, "taxonomy.tsv")
pathway_file      <- file.path(data_dir, "pathabundance_relab_unstratified.tsv")
mrg_file          <- file.path(data_dir, "MRG_merged_subtype_abundance.tsv")
pathway_stat_file <- file.path(data_dir, "stat_MetaCyc_pathways_all_comparisons_D30D60D90.csv")

required_files <- c(
  metadata_file,
  feature_file,
  taxonomy_file,
  pathway_file,
  mrg_file,
  pathway_stat_file
)

cat("\nChecking input files:\n")
print(data.frame(
  File = required_files,
  Exists = file.exists(required_files)
))

if (!all(file.exists(required_files))) {
  stop("Some input files were not found in the 'data/' directory. Please check file paths.")
}


## =========================================================
## 2. Parameters
## =========================================================

target_times <- c("D30", "D60", "D90")

groups_to_plot <- c(
  "ITM",
  "OTM1",
  "OTM2"
)

fdr_cut <- 0.05
rho_cut <- 0.50

# 筛选时对 rho 保留 2 位小数
screen_rounded_rho <- TRUE

# 特征数量控制
top_n_genera <- 25
top_n_mrg    <- 80

pseudo_count <- 1e-8

set.seed(123)


## =========================================================
## 3. Metadata
## =========================================================

meta <- read.csv(
  metadata_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!all(c("ID", "Time", "Group") %in% colnames(meta))) {
  stop("metadata.csv must contain ID, Time, and Group.")
}

meta <- meta %>%
  mutate(
    SampleID = as.character(ID),
    Time     = as.character(Time),
    Group    = as.character(Group)
  )

# 兼容 A/B/C 到 ITM/OTM1/OTM2 的映射
meta$Group <- recode(
  meta$Group,
  "A" = "ITM",
  "B" = "OTM1",
  "C" = "OTM2"
)

meta <- meta %>%
  filter(
    Time %in% target_times,
    Group %in% groups_to_plot
  )

meta$Group <- factor(
  meta$Group,
  levels = groups_to_plot
)

meta$Time <- factor(
  meta$Time,
  levels = c("D30", "D60", "D90")
)

rownames(meta) <- meta$SampleID

cat("\nSamples retained:\n")
print(table(meta$Group, meta$Time))


## =========================================================
## 4. Genus abundance
## =========================================================

first_line <- readLines(
  feature_file,
  n = 1,
  warn = FALSE
)

skip_n <- ifelse(
  grepl("^# Constructed from biom file", first_line),
  1,
  0
)

feature <- read.delim(
  feature_file,
  sep = "\t",
  header = TRUE,
  skip = skip_n,
  check.names = FALSE,
  comment.char = "",
  quote = "",
  stringsAsFactors = FALSE
)

colnames(feature)[1] <- "FeatureID"

feature <- feature %>%
  filter(
    !is.na(FeatureID),
    FeatureID != "",
    !str_detect(FeatureID, "^#")
  ) %>%
  mutate(
    FeatureID = as.character(FeatureID)
  )

sample_cols <- setdiff(
  colnames(feature),
  "FeatureID"
)

feature[sample_cols] <- lapply(
  feature[sample_cols],
  function(x) as.numeric(as.character(x))
)

# 基于全量特征计算分母
sample_total <- colSums(
  feature[, sample_cols, drop = FALSE],
  na.rm = TRUE
)


## ---------- taxonomy ----------

taxonomy <- read.delim(
  taxonomy_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  comment.char = "",
  quote = "",
  stringsAsFactors = FALSE
)

colnames(taxonomy) <- gsub("^Feature\\.ID$", "FeatureID", colnames(taxonomy))
colnames(taxonomy) <- gsub("^Feature ID$", "FeatureID", colnames(taxonomy))

if (!"FeatureID" %in% colnames(taxonomy)) {
  colnames(taxonomy)[1] <- "FeatureID"
}

if (!"Taxon" %in% colnames(taxonomy)) {
  stop("taxonomy.tsv must contain a Taxon column.")
}

taxonomy <- taxonomy %>%
  mutate(
    FeatureID = as.character(FeatureID),
    Taxon     = as.character(Taxon)
  )


## ---------- extract genus ----------

extract_genus <- function(taxon) {
  taxon <- as.character(taxon)
  parts <- unlist(strsplit(taxon, ";"))
  parts <- str_trim(parts)

  genus <- parts[str_detect(parts, "^g__|^D_5__")]

  if (length(genus) == 0) {
    return("Unclassified")
  }

  genus <- genus[1]
  genus <- str_replace(genus, "^g__", "")
  genus <- str_replace(genus, "^D_5__", "")
  genus <- str_trim(genus)

  if (
    is.na(genus) ||
    genus == "" ||
    genus == "__" ||
    str_detect(tolower(genus), "uncultured|unclassified|unknown|metagenome")
  ) {
    return("Unclassified")
  }

  genus
}

clean_genus_name <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "_[0-9]+$", "")
  x <- str_replace(x, "_[A-E]$", "")
  x
}

taxonomy <- taxonomy %>%
  mutate(
    Genus       = sapply(Taxon, extract_genus),
    Genus_clean = clean_genus_name(Genus)
  )


## ---------- remove non-interpretable genus names ----------

bad_genus_pattern <- paste(
  c(
    "^CAG", "^UBA", "^UCG", "^RUG", "^SFM", "^PeH", "^RF", "^NK",
    "^GCA", "^GCF", "^DTU", "^QAM", "^QEM", "^QHM", "^QYM", "^QZM",
    "^QGM", "^uncultured", "^Unclassified", "^unknown", "^metagenome",
    "^Ambiguous_taxa", "^Incertae_Sedis", "^Family_XIII",
    "^Lachnospiraceae_[A-Z0-9]", "^Ruminococcaceae_[A-Z0-9]",
    "^Eubacterium_", "^Blautia_[A-Z]_", "^Clostridium_[A-Z]",
    "^Bacteroides_[A-Z]", "^Prevotella_[A-Z]", "^Cryptobacteroides$",
    "^Paramuribaculum$"
  ),
  collapse = "|"
)

taxonomy <- taxonomy %>%
  mutate(
    Is_interpretable =
      !is.na(Genus_clean) &
      Genus_clean != "" &
      Genus_clean != "Unclassified" &
      !str_detect(Genus, regex(bad_genus_pattern, ignore_case = TRUE))
  )


## ---------- combine feature table and taxonomy ----------

feature_tax <- feature %>%
  left_join(
    taxonomy %>% select(FeatureID, Genus_clean, Is_interpretable),
    by = "FeatureID"
  ) %>%
  mutate(
    Genus_clean      = ifelse(is.na(Genus_clean), "Unclassified", Genus_clean),
    Is_interpretable = ifelse(is.na(Is_interpretable), FALSE, Is_interpretable)
  )


## ---------- aggregate ASVs to genus ----------

genus_table <- feature_tax %>%
  filter(Is_interpretable) %>%
  group_by(Genus_clean) %>%
  summarise(
    across(all_of(sample_cols), sum, na.rm = TRUE),
    .groups = "drop"
  )

genus_mat <- genus_table %>%
  column_to_rownames("Genus_clean") %>%
  as.matrix()

mode(genus_mat) <- "numeric"
genus_mat <- t(genus_mat)


## ---------- convert to relative abundance ----------

genus_mat <- sweep(
  genus_mat,
  1,
  sample_total[rownames(genus_mat)],
  "/"
)

genus_mat[!is.finite(genus_mat)] <- 0


## ---------- keep target samples & top genera ----------

common_genus_samples <- intersect(rownames(genus_mat), rownames(meta))

genus_post <- genus_mat[common_genus_samples, , drop = FALSE]
meta_post  <- meta[common_genus_samples, , drop = FALSE]

genus_mean  <- colMeans(genus_post, na.rm = TRUE)
keep_genera <- names(sort(genus_mean, decreasing = TRUE))[1:min(top_n_genera, length(genus_mean))]

genus_post <- genus_post[, keep_genera, drop = FALSE]

cat("\nTop genera retained:\n")
print(data.frame(
  Genus = keep_genera,
  Mean_abundance = genus_mean[keep_genera]
))


## =========================================================
## 5. HUMAnN MetaCyc pathway categories
## =========================================================

pathway_raw <- read.delim(
  pathway_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  comment.char = "",
  quote = "",
  stringsAsFactors = FALSE
)

colnames(pathway_raw)[1] <- "Pathway"

pathway_raw <- pathway_raw %>%
  filter(
    !is.na(Pathway),
    Pathway != "",
    !str_detect(Pathway, "\\|"),
    !str_detect(Pathway, "UNMAPPED|UNINTEGRATED")
  )

pathway_sample_cols <- setdiff(colnames(pathway_raw), "Pathway")
pathway_raw[pathway_sample_cols] <- lapply(
  pathway_raw[pathway_sample_cols],
  function(x) as.numeric(as.character(x))
)


## ---------- pathway descriptions ----------

path_stat <- read.csv(
  pathway_stat_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if ("Description" %in% colnames(path_stat)) {
  pathway_description <- path_stat %>%
    select(Pathway, Description) %>%
    distinct(Pathway, .keep_all = TRUE)

  pathway_raw <- pathway_raw %>%
    left_join(pathway_description, by = "Pathway") %>%
    mutate(
      Category_text = ifelse(
        is.na(Description) | Description == "",
        Pathway,
        Description
      )
    )
} else {
  pathway_raw$Category_text <- pathway_raw$Pathway
}


## ---------- classify MetaCyc pathways ----------

classify_pathway_category <- function(x) {
  y <- tolower(x)
  case_when(
    str_detect(
      y,
      paste0(
        "glycol|gluconeogen|pentose|starch|glycogen|",
        "galactose|sucrose|fructose|mannose|lactose|",
        "fermentation|pyruvate|butanoate|lactate|",
        "acetate|bifidobacterium shunt|carbohydrate"
      )
    ) ~ "Carbohydrate metabolism / fermentation",

    str_detect(
      y,
      paste0(
        "amino|lysine|methionine|cysteine|arginine|",
        "ornithine|glutamate|glutamine|alanine|",
        "valine|leucine|isoleucine|histidine|",
        "tryptophan|phenylalanine|tyrosine|",
        "polyamine|putrescine|spermidine"
      )
    ) ~ "Amino acid / polyamine metabolism",

    str_detect(
      y,
      paste0(
        "purine|pyrimidine|nucleotide|nucleoside|",
        "adenosine|guanosine|cytidine|uridine|",
        "thymidine|ump|cmp|amp|gmp|dtmp|",
        "ribonucl|deoxyribonucl"
      )
    ) ~ "Nucleotide metabolism",

    str_detect(
      y,
      paste0(
        "cofactor|nad|nadh|nadp|menaquin|quinol|",
        "heme|tetrapyrrole|biotin|thiamin|folate|",
        "riboflavin|pantothenate|cobalamin|vitamin|",
        "redox|queuosine"
      )
    ) ~ "Cofactor / redox metabolism",

    str_detect(
      y,
      paste0(
        "peptidoglycan|cell wall|cell envelope|",
        "membrane|lipid|lipopolysaccharide|o-antigen|",
        "teichoic|capsule|udp-n-acetyl|fatty acid|",
        "phospholipid"
      )
    ) ~ "Cell envelope / membrane",

    str_detect(
      y,
      paste0(
        "tca|tricarboxylic|glyoxylate|respiration|",
        "electron transfer|oxidative phosphorylation|",
        "atp|energy|aerobic|anaerobic"
      )
    ) ~ "Energy metabolism",

    TRUE ~ "Other"
  )
}

keep_categories <- c(
  "Carbohydrate metabolism / fermentation",
  "Amino acid / polyamine metabolism",
  "Nucleotide metabolism",
  "Cofactor / redox metabolism",
  "Cell envelope / membrane",
  "Energy metabolism"
)

pathway_raw <- pathway_raw %>%
  mutate(Pathway_category = classify_pathway_category(Category_text))

cat("\nNumber of pathways in each category:\n")
print(
  pathway_raw %>%
    filter(Pathway_category %in% keep_categories) %>%
    count(Pathway_category)
)


## ---------- category abundance = sum of pathways ----------

pathway_category_table <- pathway_raw %>%
  filter(Pathway_category %in% keep_categories) %>%
  group_by(Pathway_category) %>%
  summarise(
    across(all_of(pathway_sample_cols), sum, na.rm = TRUE),
    .groups = "drop"
  )

pathway_mat <- pathway_category_table %>%
  column_to_rownames("Pathway_category") %>%
  as.matrix()

mode(pathway_mat) <- "numeric"
pathway_post <- t(pathway_mat)

rownames(pathway_post) <- gsub("_Abundance$", "", rownames(pathway_post))
rownames(pathway_post) <- gsub("_merged.*$", "", rownames(pathway_post))
rownames(pathway_post) <- gsub("\\.tsv$", "", rownames(pathway_post))

cat("\nPathway categories:\n")
print(colnames(pathway_post))


## =========================================================
## 6. MRG subtype abundance
## =========================================================

mrg <- read.delim(
  mrg_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  comment.char = "",
  quote = "",
  stringsAsFactors = FALSE
)

colnames(mrg)[1] <- "Gene"
mrg <- mrg %>% mutate(Gene = as.character(Gene))

mrg_sample_cols <- setdiff(colnames(mrg), "Gene")
mrg[mrg_sample_cols] <- lapply(
  mrg[mrg_sample_cols],
  function(x) as.numeric(as.character(x))
)

mrg_mat <- mrg %>%
  group_by(Gene) %>%
  summarise(
    across(all_of(mrg_sample_cols), sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

mode(mrg_mat) <- "numeric"
mrg_post_all <- t(mrg_mat)

rownames(mrg_post_all) <- gsub("_Abundance$", "", rownames(mrg_post_all))
rownames(mrg_post_all) <- gsub("_merged.*$", "", rownames(mrg_post_all))
rownames(mrg_post_all) <- gsub("\\.tsv$", "", rownames(mrg_post_all))

mrg_samples <- intersect(rownames(mrg_post_all), rownames(meta))
mrg_tmp     <- mrg_post_all[mrg_samples, , drop = FALSE]

mrg_mean <- colMeans(mrg_tmp, na.rm = TRUE)
keep_mrg <- names(sort(mrg_mean, decreasing = TRUE))[1:min(top_n_mrg, length(mrg_mean))]

mrg_post <- mrg_post_all[, keep_mrg, drop = FALSE]

cat("\nNumber of MRG subtypes retained:", ncol(mrg_post), "\n")


## =========================================================
## 7. Match samples & Log10 Transformation
## =========================================================

common_samples <- Reduce(
  intersect,
  list(
    rownames(genus_post),
    rownames(pathway_post),
    rownames(mrg_post),
    rownames(meta)
  )
)

cat("\nCommon samples:", length(common_samples), "\n")

if (length(common_samples) < 20) {
  stop("Too few matched samples. Check sample names.")
}

genus_final   <- genus_post[common_samples, , drop = FALSE]
pathway_final <- pathway_post[common_samples, , drop = FALSE]
mrg_final     <- mrg_post[common_samples, , drop = FALSE]
meta_final    <- meta[common_samples, , drop = FALSE]

cat("\nFinal sample distribution:\n")
print(table(meta_final$Group, meta_final$Time))

genus_cor   <- log10(genus_final + pseudo_count)
pathway_cor <- log10(pathway_final + pseudo_count)
mrg_cor     <- log10(mrg_final + pseudo_count)


## =========================================================
## 8. Spearman correlation function
## =========================================================

cor_by_group <- function(mat_x, mat_y, metadata, x_type, y_type) {
  res_list <- list()

  for (grp in groups_to_plot) {
    idx <- as.character(metadata$Group) == grp
    x_g <- mat_x[idx, , drop = FALSE]
    y_g <- mat_y[idx, , drop = FALSE]

    for (x_name in colnames(x_g)) {
      for (y_name in colnames(y_g)) {
        x <- x_g[, x_name]
        y <- y_g[, y_name]

        use   <- complete.cases(x, y)
        x2    <- x[use]
        y2    <- y[use]
        n_use <- length(x2)

        if (n_use < 6 || length(unique(x2)) < 2 || length(unique(y2)) < 2) {
          rho  <- NA_real_
          pval <- NA_real_
        } else {
          ct <- suppressWarnings(
            cor.test(x2, y2, method = "spearman", exact = FALSE)
          )
          rho  <- unname(ct$estimate)
          pval <- ct$p.value
        }

        res_list[[length(res_list) + 1]] <- data.frame(
          Group     = grp,
          From      = x_name,
          To        = y_name,
          From_type = x_type,
          To_type   = y_type,
          Rho       = rho,
          P_value   = pval,
          N         = n_use,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  bind_rows(res_list)
}


## =========================================================
## 9. Calculate Correlations & BH Correction
## =========================================================

cat("\nCalculating Genus–Pathway correlations...\n")
cor_GP <- cor_by_group(genus_cor, pathway_cor, meta_final, "Genus", "Pathway category")

cat("Calculating Pathway–MRG correlations...\n")
cor_PM <- cor_by_group(pathway_cor, mrg_cor, meta_final, "Pathway category", "MRG subtype")

cat("Calculating Genus–MRG correlations...\n")
cor_GM <- cor_by_group(genus_cor, mrg_cor, meta_final, "Genus", "MRG subtype")


# BH correction
cor_GP <- cor_GP %>%
  group_by(Group) %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  ungroup()

cor_PM <- cor_PM %>%
  group_by(Group, To) %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  ungroup()

cor_GM <- cor_GM %>%
  group_by(Group) %>%
  mutate(FDR = p.adjust(P_value, method = "BH")) %>%
  ungroup()


## =========================================================
## 10. Filter Significant Edges & Save CSVs
## =========================================================

add_significance <- function(x) {
  x %>%
    mutate(
      Rho_for_screening = if (screen_rounded_rho) round(abs(Rho), 2) else abs(Rho),
      Significant       = !is.na(FDR) & FDR < fdr_cut & Rho_for_screening >= rho_cut
    )
}

cor_GP <- add_significance(cor_GP)
cor_PM <- add_significance(cor_PM)
cor_GM <- add_significance(cor_GM)

sig_GP <- cor_GP %>% filter(Significant)
sig_PM <- cor_PM %>% filter(Significant)
sig_GM <- cor_GM %>% filter(Significant)

cat("\nSignificant GP edges:", nrow(sig_GP), "\n")
cat("Significant PM edges:", nrow(sig_PM), "\n")
cat("Significant GM edges:", nrow(sig_GM), "\n")

write.csv(cor_GP, file.path(out_dir, "Fig5D_Genus_Pathway_correlations.csv"), row.names = FALSE)
write.csv(cor_PM, file.path(out_dir, "Fig5D_Pathway_MRG_correlations.csv"), row.names = FALSE)
write.csv(cor_GM, file.path(out_dir, "Fig5D_Genus_MRG_correlations.csv"), row.names = FALSE)


## =========================================================
## 11. Construct Closed Chains
## =========================================================

chains <- sig_GP %>%
  rename(
    Genus            = From,
    Pathway_category = To,
    Rho_GP           = Rho,
    FDR_GP           = FDR
  ) %>%
  select(Group, Genus, Pathway_category, Rho_GP, FDR_GP) %>%
  inner_join(
    sig_PM %>%
      rename(
        Pathway_category = From,
        MRG_subtype      = To,
        Rho_PM           = Rho,
        FDR_PM           = FDR
      ) %>%
      select(Group, Pathway_category, MRG_subtype, Rho_PM, FDR_PM),
    by = c("Group", "Pathway_category")
  ) %>%
  inner_join(
    sig_GM %>%
      rename(
        Genus       = From,
        MRG_subtype = To,
        Rho_GM      = Rho,
        FDR_GM      = FDR
      ) %>%
      select(Group, Genus, MRG_subtype, Rho_GM, FDR_GM),
    by = c("Group", "Genus", "MRG_subtype")
  ) %>%
  mutate(
    Mean_abs_rho = (abs(Rho_GP) + abs(Rho_PM) + abs(Rho_GM)) / 3,
    Max_FDR      = pmax(FDR_GP, FDR_PM, FDR_GM)
  ) %>%
  distinct(Group, Genus, Pathway_category, MRG_subtype, .keep_all = TRUE)

write.csv(
  chains,
  file.path(out_dir, "Genus_PathwayCategory_MRGsubtype_closed_chains.csv"),
  row.names = FALSE
)

cat("\nNumber of closed chains:", nrow(chains), "\n")

if (nrow(chains) == 0) {
  stop("No closed chains detected.")
}


## =========================================================
## 12. Prepare Network Nodes and Edges
## =========================================================

edges_GP <- chains %>%
  transmute(
    Group, From = Genus, To = Pathway_category,
    From_type = "Genus", To_type = "Pathway category",
    Edge_class = "Genus–Pathway", Rho = Rho_GP, FDR = FDR_GP
  )

edges_PM <- chains %>%
  transmute(
    Group, From = Pathway_category, To = MRG_subtype,
    From_type = "Pathway category", To_type = "MRG subtype",
    Edge_class = "Pathway–MRG", Rho = Rho_PM, FDR = FDR_PM
  )

edges_GM <- chains %>%
  transmute(
    Group, From = Genus, To = MRG_subtype,
    From_type = "Genus", To_type = "MRG subtype",
    Edge_class = "Genus–MRG", Rho = Rho_GM, FDR = FDR_GM
  )

edges <- bind_rows(edges_GP, edges_PM, edges_GM) %>%
  distinct(Group, From, To, From_type, To_type, .keep_all = TRUE) %>%
  mutate(
    Correlation  = ifelse(Rho > 0, "Positive correlation", "Negative correlation"),
    Relationship = ifelse(Edge_class == "Genus–MRG", "Direct Genus–MRG", "Genus–Pathway or Pathway–MRG"),
    AbsRho       = abs(Rho),
    From_key     = paste(Group, From_type, From, sep = "___"),
    To_key       = paste(Group, To_type, To, sep = "___")
  )

nodes <- bind_rows(
  chains %>% transmute(Group, Name = Genus, Node_type = "Genus"),
  chains %>% transmute(Group, Name = Pathway_category, Node_type = "Pathway category"),
  chains %>% transmute(Group, Name = MRG_subtype, Node_type = "MRG subtype")
) %>%
  distinct() %>%
  mutate(Node_key = paste(Group, Node_type, Name, sep = "___"))


## ---------- clean labels ----------

clean_mrg_label <- function(x) {
  x <- str_replace_all(x, "__", " ")
  str_replace_all(x, "_", " ")
}

clean_pathway_label <- function(x) {
  ifelse(nchar(x) > 36, paste0(substr(x, 1, 33), "..."), x)
}

nodes <- nodes %>%
  mutate(
    Label = case_when(
      Node_type == "MRG subtype"      ~ clean_mrg_label(Name),
      Node_type == "Pathway category" ~ clean_pathway_label(Name),
      TRUE                            ~ Name
    )
  )


## ---------- node degree ----------

degree_df <- bind_rows(
  edges %>% transmute(Group, Node_key = From_key),
  edges %>% transmute(Group, Node_key = To_key)
) %>%
  count(Group, Node_key, name = "Degree")

nodes <- nodes %>%
  left_join(degree_df, by = c("Group", "Node_key")) %>%
  mutate(Degree = ifelse(is.na(Degree), 1, Degree))


## =========================================================
## 13. Network Layout Generator
## =========================================================

scale_coord <- function(x) {
  if (length(x) <= 1 || is.na(sd(x)) || sd(x) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

make_group_layout <- function(group_name, nodes, edges) {
  nodes_g <- nodes %>% filter(Group == group_name)
  edges_g <- edges %>% filter(Group == group_name)

  if (nrow(nodes_g) == 0 || nrow(edges_g) == 0) return(NULL)

  vertices <- nodes_g %>%
    transmute(name = Node_key, Label = Label, Node_type = Node_type, Degree = Degree)

  edge_df <- edges_g %>%
    transmute(
      from = From_key, to = To_key, Correlation = Correlation,
      Relationship = Relationship, Rho = Rho, FDR = FDR, AbsRho = AbsRho
    )

  g <- graph_from_data_frame(d = edge_df, vertices = vertices, directed = FALSE)

  set.seed(123)
  lay <- layout_with_fr(g, weights = E(g)$AbsRho, niter = 2500)

  node_pos <- as.data.frame(lay)
  colnames(node_pos) <- c("x", "y")
  node_pos$name      <- V(g)$name
  node_pos$Label     <- V(g)$Label
  node_pos$Node_type <- V(g)$Node_type
  node_pos$Degree    <- V(g)$Degree
  node_pos$Group     <- group_name

  node_pos <- node_pos %>%
    mutate(x = scale_coord(x), y = scale_coord(y))

  edge_pos <- as_data_frame(g, what = "edges") %>%
    left_join(node_pos %>% select(from = name, x_from = x, y_from = y), by = "from") %>%
    left_join(node_pos %>% select(to = name, x_to = x, y_to = y), by = "to") %>%
    mutate(Group = group_name)

  list(nodes = node_pos, edges = edge_pos)
}

layout_list <- lapply(groups_to_plot, make_group_layout, nodes = nodes, edges = edges)
layout_list <- layout_list[!sapply(layout_list, is.null)]

node_plot <- bind_rows(lapply(layout_list, function(x) x$nodes))
edge_plot <- bind_rows(lapply(layout_list, function(x) x$edges))

node_plot$Group <- factor(node_plot$Group, levels = groups_to_plot)
edge_plot$Group <- factor(edge_plot$Group, levels = groups_to_plot)


## =========================================================
## 14. Colors & Linetypes Configuration
## =========================================================

node_colors <- c(
  "Genus"            = "#F4A62A", # orange
  "Pathway category" = "#2CA25F", # green
  "MRG subtype"      = "#2B8CBE"  # blue
)

edge_colors <- c(
  "Positive correlation" = "#D6604D",
  "Negative correlation" = "#4393C3"
)

# 直接 Genus–MRG 关联用短虚线 "22"
edge_linetypes <- c(
  "Genus–Pathway or Pathway–MRG" = "solid",
  "Direct Genus–MRG"             = "22"
)


## =========================================================
## 15. Fig. 5D Network Plotting
## =========================================================

p_network <- ggplot() +

  ## -------------------------------------------------------
  ## Edges
  ## -------------------------------------------------------
  geom_curve(
    data = edge_plot,
    aes(
      x        = x_from,
      y        = y_from,
      xend     = x_to,
      yend     = y_to,
      color    = Correlation,
      linewidth = AbsRho,
      linetype = Relationship
    ),
    curvature = 0.12,
    alpha     = 0.82,
    lineend   = "round"
  ) +

  ## -------------------------------------------------------
  ## Nodes
  ## -------------------------------------------------------
  geom_point(
    data = node_plot,
    aes(
      x    = x,
      y    = y,
      fill = Node_type,
      size = Degree
    ),
    shape  = 21,
    color  = "white",
    stroke = 0.7,
    alpha  = 1
  ) +

  ## -------------------------------------------------------
  ## Labels
  ## -------------------------------------------------------
  geom_text_repel(
    data = node_plot,
    aes(
      x     = x,
      y     = y,
      label = Label
    ),
    family            = "Arial",
    size              = 3.2,
    color             = "grey20",
    max.overlaps      = Inf,
    box.padding       = 0.22,
    point.padding     = 0.20,
    segment.color     = "grey70",
    segment.linewidth = 0.25,
    min.segment.length = 0
  ) +

  ## -------------------------------------------------------
  ## Facet Wrap
  ## -------------------------------------------------------
  facet_wrap(~ Group, nrow = 1) +

  ## -------------------------------------------------------
  ## Scales
  ## -------------------------------------------------------
  scale_fill_manual(values = node_colors, name = "Node type") +
  scale_color_manual(values = edge_colors, name = "Correlation") +
  scale_linetype_manual(values = edge_linetypes, name = "Relationship") +

  scale_linewidth_continuous(
    range  = c(0.35, 1.45),
    limits = c(0.5, 1),
    name   = "|Spearman rho|"
  ) +

  scale_size_continuous(
    range = c(5, 10),
    name  = "Node degree"
  ) +

  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(size = 6, shape = 21, color = "white")
    ),
    color = guide_legend(
      order = 2,
      override.aes = list(linewidth = 1.1, linetype = "solid")
    ),
    linetype = guide_legend(
      order = 3,
      override.aes = list(color = "grey30", linewidth = 0.9)
    ),
    linewidth = guide_legend(order = 4),
    size      = guide_legend(order = 5)
  ) +

  coord_equal() +

  theme_void(base_size = 13, base_family = "Arial") +

  theme(
    text             = element_text(family = "Arial", color = "black"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text       = element_text(face = "bold", size = 15, color = "black", family = "Arial"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key        = element_rect(fill = "white", color = NA),
    legend.title      = element_text(face = "bold", size = 10.5, color = "black"),
    legend.text       = element_text(size = 9, color = "black"),
    legend.position   = "right",
    plot.title        = element_text(face = "bold", hjust = 0.5, size = 15, color = "black"),
    plot.subtitle     = element_text(hjust = 0.5, size = 11, color = "black"),
    panel.border      = element_rect(color = "grey65", fill = NA, linewidth = 0.45),
    plot.margin       = margin(10, 10, 10, 10)
  ) +

  labs(
    title    = "Network associations among genera, pathway categories, and MRG subtypes",
    subtitle = "Closed chains based on within-group Spearman correlations; D30 + D60 + D90"
  )

print(p_network)


## =========================================================
## 16. Save Figure Output (Relative Paths)
## =========================================================

ggsave(
  file.path(out_dir, "Fig5D_closed_chain_network.pdf"),
  p_network,
  width  = 17,
  height = 6.5,
  device = cairo_pdf,
  bg     = "white"
)

ggsave(
  file.path(out_dir, "Fig5D_closed_chain_network.png"),
  p_network,
  width  = 17,
  height = 6.5,
  dpi    = 600,
  bg     = "white"
)

cat("\nAnalysis completed successfully. Output files saved to 'output/' folder.\n")