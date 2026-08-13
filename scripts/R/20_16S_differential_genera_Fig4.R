########################################################
## 16S differential genus analysis
##
## Figure 4D
##
## Post-baseline period:
##   D30 + D60 + D90
##
## Input:
##   metadata/metadata.csv
##
##   processed_data/16S/
##     feature-table.tsv
##     taxonomy.tsv
##
## Analysis:
##   1. Aggregate ASVs/features to genus level
##   2. Convert genus abundance to CPM
##   3. Pool D30, D60, and D90 samples
##   4. Pairwise comparisons:
##        ITM vs OTM1
##        ITM vs OTM2
##        OTM1 vs OTM2
##   5. Wilcoxon rank-sum test (unpaired)
##   6. Benjamini-Hochberg FDR correction
##      separately within each pairwise comparison
##
## Fold change:
##   log2FC = mean(group2) / mean(group1)
##
## Display:
##   Only classified/interpretable genera with
##   FDR < 0.05 are shown.
##
##   If more than 35 significant genera occur in one
##   comparison, retain the 35 with the largest |log2FC|.
##
## Output:
##   results/Fig4_16S_differential_genera/
########################################################


# ======================================================
# 1. Load packages
# ======================================================

library(tidyverse)
library(phyloseq)
library(patchwork)


# ======================================================
# 2. Project paths
#
# Run this script from the repository root:
# DairyCow-TraceMineral-Resistome/
# ======================================================

PROJECT_DIR <- "."


feature_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "16S",
  "feature-table.tsv"
)


taxonomy_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "16S",
  "taxonomy.tsv"
)


meta_path <- file.path(
  PROJECT_DIR,
  "metadata",
  "metadata.csv"
)


output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig4_16S_differential_genera"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Post-baseline sampling times
use_times <- c(
  "D30",
  "D60",
  "D90"
)


# Maximum number of significant genera displayed
# within each pairwise comparison
top_n <- 35


# ======================================================
# 3. Read metadata
# ======================================================

meta <- read.csv(
  meta_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


required_meta_cols <- c(
  "ID",
  "Time",
  "Group"
)


missing_meta_cols <- setdiff(
  required_meta_cols,
  colnames(meta)
)


if (length(missing_meta_cols) > 0) {

  stop(
    paste0(
      "Missing metadata column(s): ",
      paste(
        missing_meta_cols,
        collapse = ", "
      )
    )
  )
}


meta <- meta %>%

  mutate(

    ID = as.character(ID),

    Time = as.character(Time),

    Group = as.character(Group)
  )


# ======================================================
# 4. Check metadata
# ======================================================

expected_groups <- c(
  "ITM",
  "OTM1",
  "OTM2"
)


expected_times <- c(
  "D0",
  "D30",
  "D60",
  "D90"
)


unexpected_groups <- setdiff(
  unique(meta$Group),
  expected_groups
)


if (length(unexpected_groups) > 0) {

  stop(
    paste0(
      "Unexpected Group value(s): ",
      paste(
        unexpected_groups,
        collapse = ", "
      )
    )
  )
}


unexpected_times <- setdiff(
  unique(meta$Time),
  expected_times
)


if (length(unexpected_times) > 0) {

  stop(
    paste0(
      "Unexpected Time value(s): ",
      paste(
        unexpected_times,
        collapse = ", "
      )
    )
  )
}


if (anyDuplicated(meta$ID) > 0) {

  stop(
    "Duplicated sequencing sample IDs were found in metadata.csv."
  )
}


meta$Group <- factor(
  meta$Group,
  levels = expected_groups
)


meta$Time <- factor(
  meta$Time,
  levels = expected_times
)


rownames(meta) <- meta$ID


cat(
  "\n========================================\n"
)

cat(
  "16S differential genus analysis\n"
)

cat(
  "========================================\n"
)


cat(
  "\nSamples by Group × Time:\n"
)


print(
  table(
    meta$Group,
    meta$Time
  )
)


# ======================================================
# 5. Read QIIME2-exported feature table
#
# QIIME2 BIOM export may contain:
# # Constructed from biom file
#
# as the first line.
# ======================================================

first_feature_line <- readLines(
  feature_path,
  n = 1,
  warn = FALSE
)


feature_skip <- ifelse(

  grepl(
    "^# Constructed from biom file",
    first_feature_line
  ),

  1,

  0
)


otu_df <- read.delim(

  feature_path,

  header = TRUE,

  sep = "\t",

  skip = feature_skip,

  check.names = FALSE,

  comment.char = "",

  quote = "",

  stringsAsFactors = FALSE
)


colnames(otu_df)[1] <- "FeatureID"


otu_df$FeatureID <- as.character(
  otu_df$FeatureID
)


otu_sample_cols <- setdiff(
  colnames(otu_df),
  "FeatureID"
)


otu_df[
  otu_sample_cols
] <- lapply(

  otu_df[
    otu_sample_cols
  ],

  function(x) {

    x <- as.numeric(
      as.character(x)
    )


    x[
      is.na(x)
    ] <- 0


    x[
      x < 0
    ] <- 0


    x
  }
)


otu <- as.matrix(
  otu_df[
    ,
    otu_sample_cols,
    drop = FALSE
  ]
)


rownames(otu) <- otu_df$FeatureID


# ======================================================
# 6. Read taxonomy
# ======================================================

tax <- read.delim(

  taxonomy_path,

  header = TRUE,

  sep = "\t",

  fill = TRUE,

  check.names = FALSE,

  comment.char = "",

  quote = "",

  stringsAsFactors = FALSE
)


# Standardize first column
colnames(tax) <- gsub(
  "^Feature\\.ID$",
  "FeatureID",
  colnames(tax)
)


colnames(tax) <- gsub(
  "^Feature ID$",
  "FeatureID",
  colnames(tax)
)


colnames(tax) <- gsub(
  "^#OTU ID$",
  "FeatureID",
  colnames(tax)
)


if (!"FeatureID" %in% colnames(tax)) {

  colnames(tax)[1] <- "FeatureID"
}


if (!"Taxon" %in% colnames(tax)) {

  stop(
    "taxonomy.tsv must contain a Taxon column."
  )
}


tax_clean_df <- tax %>%

  select(
    FeatureID,
    Taxon
  ) %>%

  mutate(
    FeatureID = as.character(FeatureID),
    Taxon = as.character(Taxon)
  )


# ======================================================
# 7. Parse taxonomy
#
# Retain the taxonomy-cleaning logic used in the
# original analysis.
# ======================================================

clean_tax_fun <- function(x) {


  x <- gsub(
    "[a-z]__",
    "",
    x
  )


  x <- gsub(
    " ",
    "",
    x
  )


  parts <- stringr::str_split(
    x,
    ";",
    simplify = TRUE
  )


  if (
    ncol(parts) < 7
  ) {

    parts <- cbind(

      parts,

      matrix(
        NA,
        nrow = nrow(parts),
        ncol = 7 - ncol(parts)
      )
    )
  }


  if (
    ncol(parts) > 7
  ) {

    parts <- parts[
      ,
      1:7,
      drop = FALSE
    ]
  }


  return(
    parts
  )
}


tax_mat <- clean_tax_fun(
  tax_clean_df$Taxon
)


rownames(tax_mat) <- tax_clean_df$FeatureID


colnames(tax_mat) <- c(

  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)


tax_mat[
  is.na(tax_mat) |
    tax_mat == ""
] <- "Unclassified"


# ======================================================
# 8. Match feature table, taxonomy, and metadata
# ======================================================

common_features <- intersect(
  rownames(otu),
  rownames(tax_mat)
)


cat(
  "\nFeature-table features:",
  nrow(otu),
  "\n"
)


cat(
  "Taxonomy features:",
  nrow(tax_mat),
  "\n"
)


cat(
  "Matched features:",
  length(common_features),
  "\n"
)


if (length(common_features) == 0) {

  stop(
    "No common feature IDs were found between feature table and taxonomy."
  )
}


otu <- otu[
  common_features,
  ,
  drop = FALSE
]


tax_mat <- tax_mat[
  common_features,
  ,
  drop = FALSE
]


common_samples <- intersect(
  colnames(otu),
  meta$ID
)


cat(
  "\nFeature-table samples:",
  ncol(otu),
  "\n"
)


cat(
  "Metadata samples:",
  nrow(meta),
  "\n"
)


cat(
  "Matched samples:",
  length(common_samples),
  "\n"
)


if (length(common_samples) != 120) {

  warning(
    paste0(
      "Expected 120 matched samples, but ",
      length(common_samples),
      " were found."
    )
  )
}


missing_in_feature <- setdiff(
  meta$ID,
  colnames(otu)
)


if (length(missing_in_feature) > 0) {

  cat(
    "\nSamples in metadata but missing from feature table:\n"
  )

  print(
    missing_in_feature
  )
}


if (length(common_samples) < 5) {

  stop(
    "Too few common samples between feature table and metadata."
  )
}


otu <- otu[
  ,
  common_samples,
  drop = FALSE
]


meta <- meta[
  common_samples,
  ,
  drop = FALSE
]


stopifnot(
  all(
    colnames(otu) ==
      rownames(meta)
  )
)


# ======================================================
# 9. Build phyloseq object
# ======================================================

ps <- phyloseq(

  otu_table(
    otu,
    taxa_are_rows = TRUE
  ),

  tax_table(
    tax_mat
  ),

  sample_data(
    meta
  )
)


# ======================================================
# 10. Retain post-baseline samples
#
# D30 + D60 + D90
# ======================================================

cat(
  "\nGenerating D30 + D60 + D90 genus-level CPM dataset...\n"
)


post_samples <- sample_names(ps)[

  as.character(
    sample_data(ps)$Time
  ) %in%
    use_times
]


ps_post <- prune_samples(
  post_samples,
  ps
)


ps_post <- prune_taxa(
  taxa_sums(ps_post) > 0,
  ps_post
)


cat(
  "\nPost-baseline samples:",
  nsamples(ps_post),
  "\n"
)


cat(
  "Post-baseline features:",
  ntaxa(ps_post),
  "\n"
)


cat(
  "\nPost-baseline samples by Group × Time:\n"
)


post_meta_check <- as(
  sample_data(ps_post),
  "data.frame"
)


print(
  table(
    post_meta_check$Group,
    post_meta_check$Time
  )
)


# ======================================================
# 11. Aggregate to genus level
#
# NArm = FALSE retains taxa without a classified genus.
# These are included in the statistical testing universe
# and are filtered from display later.
# ======================================================

ps_genus <- tax_glom(

  ps_post,

  taxrank = "Genus",

  NArm = FALSE
)


ps_genus <- prune_taxa(
  taxa_sums(ps_genus) > 0,
  ps_genus
)


cat(
  "\nGenus-level features after tax_glom:",
  ntaxa(ps_genus),
  "\n"
)


# ======================================================
# 12. Convert genus abundance to CPM
#
# CPM = 1e6 × genus abundance / total sample abundance
# ======================================================

ps_post_cpm <- transform_sample_counts(

  ps_genus,

  function(x) {

    if (
      sum(x) == 0
    ) {

      return(
        x
      )
    }


    1e6 *
      x /
      sum(x)
  }
)


# ======================================================
# 13. Functions for genus filtering and display names
# ======================================================

is_unknown_genus <- function(x) {


  x0 <- as.character(x)


  x0[
    is.na(x0)
  ] <- "Unclassified"


  grepl(

    paste0(

      "^Unclassified$|",

      "uncultured|ambiguous|metagenome|",

      "^UBA[0-9]|^UBA$|",

      "^CAG[-_0-9]|^CAG$|",

      "^RUG[0-9]|^RUG$|",

      "^UMGS[0-9]|^UMGS$|",

      "^WR[A-Z0-9]|",

      "^QAKW[0-9]*|",

      "^SFLA[0-9]*|",

      "^SFMI[0-9]*|",

      "^OLB[0-9]*|",

      "^YIM[-_0-9]|",

      "^PeH[0-9]*|",

      "^WQUU[0-9]*|",

      "^QAKW$|^WRAY$|^WRAI$|^WRMH$"
    ),

    x0,

    ignore.case = TRUE
  )
}


clean_genus_display <- function(x) {


  x <- as.character(x)


  # Remove database numeric suffix:
  # Bifidobacterium_388775 -> Bifidobacterium
  x <- gsub(
    "_\\d+$",
    "",
    x
  )


  # Convert underscores to spaces
  x <- gsub(
    "_",
    " ",
    x
  )


  # Remove repeated spaces
  x <- gsub(
    "\\s+",
    " ",
    x
  )


  x <- trimws(
    x
  )


  return(
    x
  )
}


# ======================================================
# 14. Pairwise differential-analysis function
#
# Important:
#
# BH/FDR correction is performed across ALL genus-level
# features tested within the pairwise comparison.
#
# Classification filtering occurs only AFTER FDR
# calculation, reproducing the original analysis.
# ======================================================

run_analysis <- function(
    ps_obj,
    group1,
    group2,
    color1,
    color2,
    top_n = 35
) {


  comparison_name <- paste0(
    group1,
    "_vs_",
    group2
  )


  comparison_title <- paste0(
    group1,
    " vs ",
    group2
  )


  cat(
    "\n>>> Analyzing: ",
    comparison_title,
    " | D30 + D60 + D90\n",
    sep = ""
  )


  meta_sub <- as(
    sample_data(ps_obj),
    "data.frame"
  )


  target_samples <- rownames(meta_sub)[

    as.character(
      meta_sub$Group
    ) %in%
      c(
        group1,
        group2
      )
  ]


  if (length(target_samples) == 0) {

    warning(
      paste0(
        "No samples were found for ",
        comparison_title
      )
    )

    return(
      NULL
    )
  }


  sub_ps <- prune_samples(
    target_samples,
    ps_obj
  )


  sub_ps <- prune_taxa(
    taxa_sums(sub_ps) > 0,
    sub_ps
  )


  otu_dat <- as(
    otu_table(sub_ps),
    "matrix"
  )


  if (!taxa_are_rows(sub_ps)) {

    otu_dat <- t(
      otu_dat
    )
  }


  tax_dat <- as.data.frame(
    tax_table(sub_ps)
  ) %>%

    rownames_to_column(
      "FeatureID"
    )


  sample_meta <- as(
    sample_data(sub_ps),
    "data.frame"
  )


  # Ensure sample order agrees with the abundance matrix
  sample_meta <- sample_meta[
    colnames(otu_dat),
    ,
    drop = FALSE
  ]


  grp <- as.character(
    sample_meta$Group
  )


  # Check group sample numbers
  n_group1 <- sum(
    grp == group1
  )


  n_group2 <- sum(
    grp == group2
  )


  cat(
    "    ",
    group1,
    ": ",
    n_group1,
    " samples; ",
    group2,
    ": ",
    n_group2,
    " samples\n",
    sep = ""
  )


  if (
    n_group1 < 2 ||
    n_group2 < 2
  ) {

    warning(
      paste0(
        "Insufficient samples for ",
        comparison_title
      )
    )

    return(
      NULL
    )
  }


  # ----------------------------------------------------
  # Test every genus-level feature
  # ----------------------------------------------------

  result_list <- vector(
    "list",
    nrow(otu_dat)
  )


  for (
    i in seq_len(
      nrow(otu_dat)
    )
  ) {


    vals <- as.numeric(
      otu_dat[
        i,
        ]
    )


    vals_group1 <- vals[
      grp == group1
    ]


    vals_group2 <- vals[
      grp == group2
    ]


    mean_group1 <- mean(
      vals_group1,
      na.rm = TRUE
    )


    mean_group2 <- mean(
      vals_group2,
      na.rm = TRUE
    )


    # log2FC = group2 / group1
    log2fc <- log2(

      (
        mean_group2 +
          1e-6
      ) /

        (
          mean_group1 +
            1e-6
        )
    )


    pval <- tryCatch(

      wilcox.test(

        vals_group2,

        vals_group1,

        paired = FALSE

      )$p.value,

      error = function(e) {

        NA_real_
      }
    )


    result_list[
      [
        i
      ]
    ] <- tibble(

      FeatureID =
        rownames(
          otu_dat
        )[
          i
        ],

      Log2FC =
        log2fc,

      Pvalue =
        pval,

      Mean_group1 =
        mean_group1,

      Mean_group2 =
        mean_group2,

      N_group1 =
        length(
          vals_group1
        ),

      N_group2 =
        length(
          vals_group2
        )
    )
  }


  # ----------------------------------------------------
  # BH correction
  #
  # Performed before genus classification filtering.
  # ----------------------------------------------------

  res_df <- bind_rows(
    result_list
  ) %>%

    mutate(

      FDR = p.adjust(
        Pvalue,
        method = "BH"
      ),

      Group1 =
        group1,

      Group2 =
        group2,

      Comparison =
        comparison_name
    )


  # ----------------------------------------------------
  # Add taxonomy
  # ----------------------------------------------------

  all_stats <- res_df %>%

    left_join(
      tax_dat,
      by = "FeatureID"
    )


  # ----------------------------------------------------
  # Retain significant classified genera
  # ----------------------------------------------------

  significant_classified <- all_stats %>%

    filter(
      !is.na(Genus)
    ) %>%

    filter(
      !is_unknown_genus(
        Genus
      )
    ) %>%

    mutate(

      Genus_display =
        clean_genus_display(
          Genus
        ),

      Enriched = ifelse(
        Log2FC > 0,
        group2,
        group1
      )
    ) %>%

    filter(
      !is.na(FDR),
      FDR < 0.05
    )


  # ----------------------------------------------------
  # Remove duplicated DISPLAY names
  #
  # This retains the original analysis behavior:
  # if multiple original genus labels collapse to the
  # same cleaned display name, retain the result with
  # the lowest FDR and then the largest |log2FC|.
  # ----------------------------------------------------

  significant_classified <- significant_classified %>%

    group_by(
      Genus_display
    ) %>%

    arrange(
      FDR,
      desc(
        abs(
          Log2FC
        )
      ),
      .by_group = TRUE
    ) %>%

    slice_head(
      n = 1
    ) %>%

    ungroup()


  # ----------------------------------------------------
  # No significant classified genera
  # ----------------------------------------------------

  if (
    nrow(
      significant_classified
    ) == 0
  ) {


    cat(
      "    No significant classified genera at FDR < 0.05.\n"
    )


    return(

      list(

        plot =
          NULL,

        all_stats =
          all_stats,

        significant =
          significant_classified,

        plotted =
          significant_classified
      )
    )
  }


  # ----------------------------------------------------
  # Select genera displayed in figure
  #
  # Largest |log2FC| among significant genera.
  # ----------------------------------------------------

  plot_res <- significant_classified %>%

    arrange(
      desc(
        abs(
          Log2FC
        )
      )
    ) %>%

    slice_head(
      n = top_n
    )


  # ----------------------------------------------------
  # Plot
  # ----------------------------------------------------

  p <- ggplot(

    plot_res,

    aes(

      x =
        reorder(
          Genus_display,
          Log2FC
        ),

      y =
        Log2FC,

      fill =
        Enriched
    )

  ) +


    geom_col(

      width =
        0.7,

      alpha =
        0.9
    ) +


    geom_hline(

      yintercept =
        0,

      linewidth =
        0.35,

      color =
        "black"
    ) +


    coord_flip() +


    scale_fill_manual(

      values =
        setNames(

          c(
            color1,
            color2
          ),

          c(
            group1,
            group2
          )
        )
    ) +


    labs(

      title =
        comparison_title,

      x =
        NULL,

      y =
        "Log2 fold change",

      fill =
        NULL
    ) +


    theme_bw(
      base_size = 11
    ) +


    theme(

      plot.title =
        element_text(
          face = "bold",
          hjust = 0.5
        ),

      axis.text.y =
        element_text(
          face = "italic",
          color = "black"
        ),

      axis.text.x =
        element_text(
          color = "black"
        ),

      axis.title =
        element_text(
          face = "bold"
        ),

      legend.position =
        "none",

      panel.grid.major.y =
        element_blank(),

      panel.grid.minor =
        element_blank()
    )


  cat(
    "    Significant classified genera: ",
    nrow(
      significant_classified
    ),
    "\n",
    sep = ""
  )


  return(

    list(

      plot =
        p,

      all_stats =
        all_stats,

      significant =
        significant_classified,

      plotted =
        plot_res
    )
  )
}


# ======================================================
# 15. Group colors
# ======================================================

group_cols <- c(

  ITM =
    "#8ECFC9",

  OTM1 =
    "#FFBE7A",

  OTM2 =
    "#FA7F6F"
)


# ======================================================
# 16. Run pairwise comparisons
# ======================================================

res_ITM_vs_OTM1 <- run_analysis(

  ps_obj =
    ps_post_cpm,

  group1 =
    "ITM",

  group2 =
    "OTM1",

  color1 =
    group_cols[
      "ITM"
    ],

  color2 =
    group_cols[
      "OTM1"
    ],

  top_n =
    top_n
)


res_ITM_vs_OTM2 <- run_analysis(

  ps_obj =
    ps_post_cpm,

  group1 =
    "ITM",

  group2 =
    "OTM2",

  color1 =
    group_cols[
      "ITM"
    ],

  color2 =
    group_cols[
      "OTM2"
    ],

  top_n =
    top_n
)


res_OTM1_vs_OTM2 <- run_analysis(

  ps_obj =
    ps_post_cpm,

  group1 =
    "OTM1",

  group2 =
    "OTM2",

  color1 =
    group_cols[
      "OTM1"
    ],

  color2 =
    group_cols[
      "OTM2"
    ],

  top_n =
    top_n
)


result_list <- list(

  ITM_vs_OTM1 =
    res_ITM_vs_OTM1,

  ITM_vs_OTM2 =
    res_ITM_vs_OTM2,

  OTM1_vs_OTM2 =
    res_OTM1_vs_OTM2
)


# ======================================================
# 17. Combine complete statistical results
# ======================================================

all_stats_table <- bind_rows(

  lapply(

    result_list,

    function(x) {

      if (is.null(x)) {

        return(
          NULL
        )
      }


      x$all_stats
    }
  )
)


write.csv(

  all_stats_table,

  file.path(
    output_dir,
    "stat_D30D60D90_genus_Wilcoxon_BH_all.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 18. Combine significant classified genera
# ======================================================

significant_table <- bind_rows(

  lapply(

    result_list,

    function(x) {

      if (is.null(x)) {

        return(
          NULL
        )
      }


      x$significant
    }
  )
)


write.csv(

  significant_table,

  file.path(
    output_dir,
    "stat_D30D60D90_genus_Wilcoxon_BH_FDR_lt_0.05_classified.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 19. Combine exact genera displayed in the figure
# ======================================================

plot_table <- bind_rows(

  lapply(

    result_list,

    function(x) {

      if (is.null(x)) {

        return(
          NULL
        )
      }


      x$plotted
    }
  )
)


write.csv(

  plot_table,

  file.path(
    output_dir,
    "plot_D30D60D90_differential_genera_data.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 20. Print significant-genus summary
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "Significant classified genera summary\n"
)

cat(
  "========================================\n"
)


if (
  nrow(
    significant_table
  ) > 0
) {


  sig_summary <- significant_table %>%

    count(
      Comparison,
      name = "n_significant_classified_genera"
    )


  print(
    sig_summary
  )


  write.csv(

    sig_summary,

    file.path(
      output_dir,
      "stat_D30D60D90_significant_genus_summary.csv"
    ),

    row.names = FALSE
  )


} else {


  cat(
    "No significant classified genera detected in any comparison.\n"
  )
}


# ======================================================
# 21. Collect non-empty plots
# ======================================================

plot_list <- list()


for (
  comparison_name in names(
    result_list
  )
) {


  result <- result_list[
    [
      comparison_name
    ]
  ]


  if (
    !is.null(result) &&
    !is.null(result$plot)
  ) {

    plot_list[
      [
        length(
          plot_list
        ) +
          1
      ]
    ] <- result$plot
  }
}


# ======================================================
# 22. Combine plots
#
# Only comparisons containing at least one significant
# classified genus are displayed.
# ======================================================

if (
  length(
    plot_list
  ) == 0
) {


  cat(
    "\nNo figure generated because no comparison contained ",
    "significant classified genera at FDR < 0.05.\n",
    sep = ""
  )


  final_plot <- NULL


} else {


  final_plot <- patchwork::wrap_plots(

    plot_list,

    nrow =
      1
  ) +


    patchwork::plot_annotation(

      title =
        "Differential genera during the post-baseline period",

      subtitle =
        paste0(
          "D30 + D60 + D90; Wilcoxon rank-sum test with ",
          "Benjamini-Hochberg correction; FDR < 0.05"
        ),

      theme =
        theme(

          plot.title =
            element_text(
              size = 16,
              face = "bold",
              hjust = 0.5
            ),

          plot.subtitle =
            element_text(
              size = 11,
              hjust = 0.5
            )
        )
    )


  print(
    final_plot
  )


  # ====================================================
  # 23. Save combined figure
  # ====================================================

  figure_width <- max(

    6,

    5 *
      length(
        plot_list
      )
  )


  ggsave(

    filename = file.path(
      output_dir,
      "Figure_16S_D30D60D90_differential_genera.pdf"
    ),

    plot =
      final_plot,

    width =
      figure_width,

    height =
      6
  )


  ggsave(

    filename = file.path(
      output_dir,
      "Figure_16S_D30D60D90_differential_genera.png"
    ),

    plot =
      final_plot,

    width =
      figure_width,

    height =
      6,

    dpi =
      600
  )
}


# ======================================================
# 24. Finish
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "16S differential genus analysis completed.\n"
)

cat(
  "Samples used: D30 + D60 + D90.\n"
)

cat(
  "Pairwise comparisons:\n"
)

cat(
  "  ITM vs OTM1\n"
)

cat(
  "  ITM vs OTM2\n"
)

cat(
  "  OTM1 vs OTM2\n"
)

cat(
  "Wilcoxon rank-sum test (unpaired).\n"
)

cat(
  "BH/FDR correction applied within each comparison.\n"
)

cat(
  "Feature table:\n"
)

cat(
  feature_path,
  "\n"
)

cat(
  "Taxonomy:\n"
)

cat(
  taxonomy_path,
  "\n"
)

cat(
  "Output directory:\n"
)

cat(
  output_dir,
  "\n"
)

cat(
  "========================================\n"
)