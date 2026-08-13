########################################################
## HUMAnN MetaCyc pathway beta-diversity analysis
##
## MetaCyc pathway PCoA + pairwise PERMANOVA
##
## Input:
##   metadata/metadata.csv
##
##   processed_data/HUMAnN/
##     pathabundance_relab_unstratified.tsv
##
## Analysis:
##   Time points:
##     D0, D30, D60, D90
##
##   For each time point:
##     1. Bray-Curtis dissimilarity
##     2. PCoA
##     3. Pairwise PERMANOVA:
##          ITM vs OTM1
##          ITM vs OTM2
##          OTM1 vs OTM2
##
## PERMANOVA:
##   999 permutations
##
## Multiple-testing correction:
##   Benjamini-Hochberg (BH/FDR)
##   applied across the three pairwise comparisons
##   separately within each time point
##
## Pathway filtering:
##   Mean relative abundance >= 1e-6
##   Prevalence >= 5 samples
##
## Plot:
##   Raw pairwise PERMANOVA P values are displayed,
##   consistent with the original analysis.
##
## Output:
##   results/Fig5_HUMAnN_pathway_PCoA/
########################################################


# ======================================================
# 1. Load packages
# ======================================================

library(tidyverse)
library(vegan)
library(patchwork)


# ======================================================
# 2. Project paths
#
# Run this script from the repository root:
# DairyCow-TraceMineral-Resistome/
# ======================================================

PROJECT_DIR <- "."


pathway_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "HUMAnN",
  "pathabundance_relab_unstratified.tsv"
)


meta_path <- file.path(
  PROJECT_DIR,
  "metadata",
  "metadata.csv"
)


output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig5_HUMAnN_pathway_PCoA"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ======================================================
# 3. Analysis settings
# ======================================================

time_points <- c(
  "D0",
  "D30",
  "D60",
  "D90"
)


group_levels <- c(
  "ITM",
  "OTM1",
  "OTM2"
)


group_pairs <- list(

  c(
    "ITM",
    "OTM1"
  ),

  c(
    "ITM",
    "OTM2"
  ),

  c(
    "OTM1",
    "OTM2"
  )
)


group_colors <- c(

  ITM =
    "#8ECFC9",

  OTM1 =
    "#FFBE7A",

  OTM2 =
    "#FA7F6F"
)


# ======================================================
# 4. Read metadata
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
# 5. Check metadata
# ======================================================

unexpected_groups <- setdiff(
  unique(meta$Group),
  group_levels
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
  time_points
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
  levels = group_levels
)


meta$Time <- factor(
  meta$Time,
  levels = time_points
)


rownames(meta) <- meta$ID


cat(
  "\n========================================\n"
)

cat(
  "HUMAnN MetaCyc pathway PCoA analysis\n"
)

cat(
  "========================================\n"
)


cat(
  "\nMetadata samples:",
  nrow(meta),
  "\n"
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
# 6. Read HUMAnN pathway table
#
# Expected:
#   first column = pathway
#   remaining columns = samples
#
# Input should be the unstratified relative-abundance
# pathway table generated from HUMAnN output.
# ======================================================

pathway <- read.delim(

  pathway_path,

  sep = "\t",

  header = TRUE,

  check.names = FALSE,

  comment.char = "",

  quote = "",

  stringsAsFactors = FALSE
)


colnames(pathway)[1] <- "Pathway"


pathway$Pathway <- as.character(
  pathway$Pathway
)


# ======================================================
# 7. Remove HUMAnN special rows
# ======================================================

pathway <- pathway %>%

  filter(
    !Pathway %in%
      c(
        "UNMAPPED",
        "UNINTEGRATED",
        "UNGROUPED"
      )
  )


# ======================================================
# 8. Convert sample columns to numeric
# ======================================================

sample_cols <- setdiff(
  colnames(pathway),
  "Pathway"
)


pathway[
  sample_cols
] <- lapply(

  pathway[
    sample_cols
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


# ======================================================
# 9. Merge duplicated pathway IDs
# ======================================================

pathway <- pathway %>%

  group_by(
    Pathway
  ) %>%

  summarise(

    across(
      all_of(
        sample_cols
      ),
      ~ sum(
        .x,
        na.rm = TRUE
      )
    ),

    .groups =
      "drop"
  )


cat(
  "\nUnique MetaCyc pathways before filtering:",
  nrow(pathway),
  "\n"
)


# ======================================================
# 10. Construct sample × pathway matrix
# ======================================================

pathway_mat <- pathway %>%

  column_to_rownames(
    "Pathway"
  ) %>%

  as.matrix()


storage.mode(
  pathway_mat
) <- "numeric"


# pathway × sample -> sample × pathway
pathway_mat <- t(
  pathway_mat
)


# ======================================================
# 11. Clean HUMAnN sample names
#
# Retain the cleaning rules from the original script.
# ======================================================

original_sample_names <- rownames(
  pathway_mat
)


clean_sample_names <- original_sample_names


clean_sample_names <- gsub(
  "_Abundance$",
  "",
  clean_sample_names
)


clean_sample_names <- gsub(
  "_merged.*$",
  "",
  clean_sample_names
)


clean_sample_names <- gsub(
  "\\.tsv$",
  "",
  clean_sample_names
)


rownames(pathway_mat) <- clean_sample_names


if (
  anyDuplicated(
    rownames(
      pathway_mat
    )
  ) > 0
) {

  stop(
    paste0(
      "Duplicated sample IDs were generated after ",
      "cleaning HUMAnN sample names."
    )
  )
}


cat(
  "\nHUMAnN pathway table:\n"
)


cat(
  "Samples:",
  nrow(pathway_mat),
  "\n"
)


cat(
  "Pathways:",
  ncol(pathway_mat),
  "\n"
)


# ======================================================
# 12. Export sample-name mapping
# ======================================================

sample_name_mapping <- tibble(

  Original_sample_name =
    original_sample_names,

  ID =
    clean_sample_names
)


write.csv(

  sample_name_mapping,

  file.path(
    output_dir,
    "HUMAnN_sample_name_mapping.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 13. Match pathway data and metadata
# ======================================================

common_samples <- intersect(
  rownames(
    pathway_mat
  ),
  meta$ID
)


cat(
  "\nMatched samples:",
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


missing_in_pathway <- setdiff(
  meta$ID,
  rownames(
    pathway_mat
  )
)


if (length(missing_in_pathway) > 0) {

  cat(
    "\nSamples in metadata but missing from HUMAnN pathway table:\n"
  )

  print(
    missing_in_pathway
  )
}


missing_in_metadata <- setdiff(
  rownames(
    pathway_mat
  ),
  meta$ID
)


if (length(missing_in_metadata) > 0) {

  cat(
    "\nSamples in HUMAnN table but missing from metadata:\n"
  )

  print(
    missing_in_metadata
  )
}


if (length(common_samples) < 5) {

  stop(
    paste0(
      "Too few matched samples. ",
      "Please check HUMAnN sample names and metadata IDs."
    )
  )
}


pathway_mat <- pathway_mat[
  common_samples,
  ,
  drop = FALSE
]


meta <- meta[
  match(
    common_samples,
    meta$ID
  ),
  ,
  drop = FALSE
]


stopifnot(
  all(
    rownames(
      pathway_mat
    ) ==
      meta$ID
  )
)


# ======================================================
# 14. Remove globally absent pathways
# ======================================================

pathway_mat <- pathway_mat[
  ,
  colSums(
    pathway_mat,
    na.rm = TRUE
  ) > 0,
  drop = FALSE
]


cat(
  "\nPathways after removing all-zero features:",
  ncol(pathway_mat),
  "\n"
)


# ======================================================
# 15. Global pathway filtering
#
# Retain original filtering criteria:
#
# mean relative abundance >= 1e-6
# prevalence >= 5 samples
# ======================================================

mean_abundance <- colMeans(
  pathway_mat,
  na.rm = TRUE
)


prevalence <- colSums(
  pathway_mat > 0,
  na.rm = TRUE
)


pathway_filter_table <- tibble(

  Pathway =
    colnames(
      pathway_mat
    ),

  Mean_relative_abundance =
    mean_abundance,

  Prevalence =
    prevalence,

  Retained =
    mean_abundance >= 1e-6 &
    prevalence >= 5
)


write.csv(

  pathway_filter_table,

  file.path(
    output_dir,
    "stat_MetaCyc_pathway_filtering.csv"
  ),

  row.names = FALSE
)


keep_pathways <- pathway_filter_table %>%

  filter(
    Retained
  ) %>%

  pull(
    Pathway
  )


pathway_filtered <- pathway_mat[
  ,
  keep_pathways,
  drop = FALSE
]


cat(
  "\nPathways retained after filtering:",
  ncol(pathway_filtered),
  "\n"
)


if (
  ncol(
    pathway_filtered
  ) < 2
) {

  stop(
    "Too few MetaCyc pathways remained after filtering."
  )
}


# ======================================================
# 16. Check samples after filtering
#
# Samples with zero total pathway abundance cannot be
# used for Bray-Curtis distance analysis.
# ======================================================

sample_keep <- rowSums(
  pathway_filtered,
  na.rm = TRUE
) > 0


zero_samples <- rownames(
  pathway_filtered
)[
  !sample_keep
]


if (length(zero_samples) > 0) {

  cat(
    "\nSamples removed because total pathway abundance was zero:\n"
  )


  print(
    zero_samples
  )


  write.csv(

    tibble(
      ID = zero_samples,
      Reason = "Zero total filtered MetaCyc pathway abundance"
    ),

    file.path(
      output_dir,
      "stat_MetaCyc_excluded_zero_samples.csv"
    ),

    row.names = FALSE
  )
}


pathway_filtered <- pathway_filtered[
  sample_keep,
  ,
  drop = FALSE
]


meta <- meta[
  sample_keep,
  ,
  drop = FALSE
]


cat(
  "\nFinal samples by Group × Time:\n"
)


print(
  table(
    meta$Group,
    meta$Time
  )
)


# ======================================================
# 17. Pairwise PERMANOVA function
#
# Three pairwise comparisons are run independently.
#
# Dissimilarity:
#   Bray-Curtis
#
# Permutations:
#   999
# ======================================================

pairwise_permanova <- function(
    distance_matrix,
    metadata
) {


  result_list <- list()


  for (
    pair in group_pairs
  ) {


    group1 <- pair[1]

    group2 <- pair[2]


    samples_use <- metadata$ID[
      as.character(
        metadata$Group
      ) %in%
        c(
          group1,
          group2
        )
    ]


    meta_sub <- metadata[
      match(
        samples_use,
        metadata$ID
      ),
      ,
      drop = FALSE
    ]


    meta_sub$Group <- droplevels(
      meta_sub$Group
    )


    if (
      length(
        unique(
          meta_sub$Group
        )
      ) < 2
    ) {

      next
    }


    dist_sub <- as.dist(

      distance_matrix[
        samples_use,
        samples_use,
        drop = FALSE
      ]
    )


    set.seed(
      123
    )


    adonis_res <- vegan::adonis2(

      dist_sub ~ Group,

      data =
        meta_sub,

      permutations =
        999
    )


    result_list[
      [
        paste(
          group1,
          group2,
          sep = "_vs_"
        )
      ]
    ] <- tibble(

      Group1 =
        group1,

      Group2 =
        group2,

      Comparison =
        paste(
          group1,
          "vs",
          group2
        ),

      F =
        adonis_res$F[
          1
        ],

      R2 =
        adonis_res$R2[
          1
        ],

      P =
        adonis_res$`Pr(>F)`[
          1
        ],

      N_group1 =
        sum(
          meta_sub$Group ==
            group1
        ),

      N_group2 =
        sum(
          meta_sub$Group ==
            group2
        ),

      Permutations =
        999
    )
  }


  bind_rows(
    result_list
  )
}


# ======================================================
# 18. P-value formatting
# ======================================================

format_p <- function(p) {


  if (is.na(p)) {

    return(
      "NA"
    )
  }


  if (p < 0.001) {

    return(
      "< 0.001"
    )
  }


  sprintf(
    "%.3f",
    p
  )
}


# ======================================================
# 19. Initialize outputs
# ======================================================

plot_list <- list()


permanova_list <- list()


coordinates_list <- list()


variance_list <- list()


# ======================================================
# 20. Analyze each time point separately
# ======================================================

for (
  time_point in time_points
) {


  cat(
    "\n========================================\n"
  )


  cat(
    "Processing:",
    time_point,
    "\n"
  )


  cat(
    "========================================\n"
  )


  # ----------------------------------------------------
  # 20.1 Select samples
  # ----------------------------------------------------

  samples_time <- meta$ID[
    as.character(
      meta$Time
    ) ==
      time_point
  ]


  meta_time <- meta[
    match(
      samples_time,
      meta$ID
    ),
    ,
    drop = FALSE
  ]


  meta_time$Group <- factor(
    meta_time$Group,
    levels = group_levels
  )


  pathway_time <- pathway_filtered[
    samples_time,
    ,
    drop = FALSE
  ]


  cat(
    "\nSamples:\n"
  )


  print(
    table(
      meta_time$Group
    )
  )


  # ----------------------------------------------------
  # 20.2 Remove pathways absent at this time point
  # ----------------------------------------------------

  pathway_time <- pathway_time[
    ,
    colSums(
      pathway_time,
      na.rm = TRUE
    ) > 0,
    drop = FALSE
  ]


  cat(
    "Pathways:",
    ncol(pathway_time),
    "\n"
  )


  if (
    nrow(
      pathway_time
    ) < 3
  ) {

    warning(
      paste0(
        time_point,
        ": too few samples; time point skipped."
      )
    )

    next
  }


  if (
    ncol(
      pathway_time
    ) < 2
  ) {

    warning(
      paste0(
        time_point,
        ": too few pathways; time point skipped."
      )
    )

    next
  }


  # ----------------------------------------------------
  # 20.3 Bray-Curtis
  # ----------------------------------------------------

  bray_time <- vegan::vegdist(

    pathway_time,

    method =
      "bray"
  )


  bray_matrix <- as.matrix(
    bray_time
  )


  # ----------------------------------------------------
  # 20.4 PCoA
  #
  # Retain cmdscale(), consistent with original analysis.
  # ----------------------------------------------------

  pcoa_time <- cmdscale(

    bray_time,

    eig =
      TRUE,

    k =
      2
  )


  eigenvalues <- pcoa_time$eig


  positive_eigenvalues <- eigenvalues[
    eigenvalues > 0
  ]


  if (
    length(
      positive_eigenvalues
    ) < 2
  ) {

    warning(
      paste0(
        time_point,
        ": fewer than two positive PCoA eigenvalues."
      )
    )
  }


  total_positive_eigenvalue <- sum(
    positive_eigenvalues
  )


  pc1_variance <- 100 *
    eigenvalues[1] /
    total_positive_eigenvalue


  pc2_variance <- 100 *
    eigenvalues[2] /
    total_positive_eigenvalue


  cat(
    "PCoA1:",
    round(
      pc1_variance,
      2
    ),
    "%\n"
  )


  cat(
    "PCoA2:",
    round(
      pc2_variance,
      2
    ),
    "%\n"
  )


  # ----------------------------------------------------
  # 20.5 PCoA coordinates
  # ----------------------------------------------------

  pcoa_df <- tibble(

    ID =
      rownames(
        pcoa_time$points
      ),

    PCoA1 =
      pcoa_time$points[
        ,
        1
      ],

    PCoA2 =
      pcoa_time$points[
        ,
        2
      ]
  ) %>%

    left_join(

      meta_time %>%

        select(
          ID,
          Time,
          Group
        ),

      by =
        "ID"
    )


  pcoa_df$Group <- factor(
    pcoa_df$Group,
    levels = group_levels
  )


  coordinates_list[
    [
      time_point
    ]
  ] <- pcoa_df %>%

    mutate(

      PCoA1_variance_percent =
        pc1_variance,

      PCoA2_variance_percent =
        pc2_variance
    )


  variance_list[
    [
      time_point
    ]
  ] <- tibble(

    Time =
      time_point,

    PCoA1_variance_percent =
      pc1_variance,

    PCoA2_variance_percent =
      pc2_variance,

    N_samples =
      nrow(
        pathway_time
      ),

    N_pathways =
      ncol(
        pathway_time
      )
  )


  # ----------------------------------------------------
  # 20.6 Pairwise PERMANOVA
  # ----------------------------------------------------

  pairwise_results <- pairwise_permanova(

    distance_matrix =
      bray_matrix,

    metadata =
      meta_time
  )


  pairwise_results$Time <- time_point


  permanova_list[
    [
      time_point
    ]
  ] <- pairwise_results


  # ----------------------------------------------------
  # 20.7 Extract raw P values for figure annotation
  # ----------------------------------------------------

  p_ITM_OTM1 <- pairwise_results$P[
    pairwise_results$Comparison ==
      "ITM vs OTM1"
  ]


  p_ITM_OTM2 <- pairwise_results$P[
    pairwise_results$Comparison ==
      "ITM vs OTM2"
  ]


  p_OTM1_OTM2 <- pairwise_results$P[
    pairwise_results$Comparison ==
      "OTM1 vs OTM2"
  ]


  title_text <- paste0(

    time_point,

    "\n",

    "ITM vs OTM1: P ",
    ifelse(
      p_ITM_OTM1 < 0.001,
      "< 0.001",
      paste0(
        "= ",
        format_p(
          p_ITM_OTM1
        )
      )
    ),

    " | ITM vs OTM2: P ",
    ifelse(
      p_ITM_OTM2 < 0.001,
      "< 0.001",
      paste0(
        "= ",
        format_p(
          p_ITM_OTM2
        )
      )
    ),

    " | OTM1 vs OTM2: P ",
    ifelse(
      p_OTM1_OTM2 < 0.001,
      "< 0.001",
      paste0(
        "= ",
        format_p(
          p_OTM1_OTM2
        )
      )
    )
  )


  # ----------------------------------------------------
  # 20.8 Plot
  # ----------------------------------------------------

  p <- ggplot(

    pcoa_df,

    aes(

      x =
        PCoA1,

      y =
        PCoA2,

      color =
        Group
    )

  ) +


    geom_point(

      size =
        3.2,

      alpha =
        0.90
    ) +


    stat_ellipse(

      aes(
        group =
          Group
      ),

      type =
        "t",

      linetype =
        2,

      linewidth =
        0.6,

      show.legend =
        FALSE
    ) +


    scale_color_manual(

      values =
        group_colors,

      limits =
        group_levels,

      drop =
        FALSE
    ) +


    labs(

      x =
        paste0(
          "PCoA1 (",
          round(
            pc1_variance,
            2
          ),
          "%)"
        ),

      y =
        paste0(
          "PCoA2 (",
          round(
            pc2_variance,
            2
          ),
          "%)"
        ),

      color =
        "Group",

      title =
        title_text
    ) +


    theme_bw(
      base_size = 12
    ) +


    theme(

      panel.grid =
        element_blank(),


      plot.title =
        element_text(
          size = 8,
          face = "bold"
        ),


      axis.title =
        element_text(
          face = "bold"
        ),


      axis.text =
        element_text(
          color = "black"
        ),


      legend.position =
        "right"
    )


  print(
    p
  )


  plot_list[
    [
      time_point
    ]
  ] <- p


  # ----------------------------------------------------
  # 20.9 Save individual time-point figures
  # ----------------------------------------------------

  ggsave(

    filename = file.path(
      output_dir,
      paste0(
        "Figure_MetaCyc_pathway_PCoA_",
        time_point,
        "_pairwise_PERMANOVA.pdf"
      )
    ),

    plot =
      p,

    width =
      6.5,

    height =
      5
  )


  ggsave(

    filename = file.path(
      output_dir,
      paste0(
        "Figure_MetaCyc_pathway_PCoA_",
        time_point,
        "_pairwise_PERMANOVA.png"
      )
    ),

    plot =
      p,

    width =
      6.5,

    height =
      5,

    dpi =
      600
  )
}


# ======================================================
# 21. Combine pairwise PERMANOVA results
#
# BH correction is performed separately within
# each time point across the three comparisons.
# ======================================================

permanova_df <- bind_rows(
  permanova_list
) %>%

  group_by(
    Time
  ) %>%

  mutate(

    FDR =
      p.adjust(
        P,
        method = "BH"
      )
  ) %>%

  ungroup() %>%

  mutate(

    Time = factor(
      Time,
      levels = time_points
    ),

    Comparison = factor(

      Comparison,

      levels = c(
        "ITM vs OTM1",
        "ITM vs OTM2",
        "OTM1 vs OTM2"
      )
    )
  ) %>%

  arrange(
    Time,
    Comparison
  )


# ======================================================
# 22. Export pairwise PERMANOVA results
# ======================================================

write.csv(

  permanova_df,

  file.path(
    output_dir,
    "stat_MetaCyc_pathway_pairwise_PERMANOVA_by_timepoint.csv"
  ),

  row.names = FALSE
)


cat(
  "\n========================================\n"
)

cat(
  "Pairwise PERMANOVA results\n"
)

cat(
  "========================================\n"
)


print(
  permanova_df
)


# ======================================================
# 23. Export PCoA coordinates
# ======================================================

coordinates_df <- bind_rows(
  coordinates_list
)


write.csv(

  coordinates_df,

  file.path(
    output_dir,
    "stat_MetaCyc_pathway_PCoA_coordinates.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 24. Export PCoA explained variation
# ======================================================

variance_df <- bind_rows(
  variance_list
)


write.csv(

  variance_df,

  file.path(
    output_dir,
    "stat_MetaCyc_pathway_PCoA_explained_variance.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 25. Combined 2 × 2 figure
# ======================================================

required_plots <- c(
  "D0",
  "D30",
  "D60",
  "D90"
)


missing_plots <- setdiff(
  required_plots,
  names(
    plot_list
  )
)


if (length(missing_plots) > 0) {

  warning(
    paste0(
      "Combined figure was not generated because ",
      "the following time point(s) were missing: ",
      paste(
        missing_plots,
        collapse = ", "
      )
    )
  )


  combined_plot <- NULL


} else {


  combined_plot <- (

    plot_list[
      [
        "D0"
      ]
    ] +

      plot_list[
        [
          "D30"
        ]
      ]

  ) / (

    plot_list[
      [
        "D60"
      ]
    ] +

      plot_list[
        [
          "D90"
        ]
      ]

  ) +


    patchwork::plot_layout(

      guides =
        "collect"
    ) &


    theme(

      legend.position =
        "right"
    )


  print(
    combined_plot
  )


  # ====================================================
  # 26. Save combined figure
  # ====================================================

  ggsave(

    filename = file.path(
      output_dir,
      "Figure_MetaCyc_pathway_PCoA_D0_D30_D60_D90_combined.pdf"
    ),

    plot =
      combined_plot,

    width =
      13,

    height =
      10
  )


  ggsave(

    filename = file.path(
      output_dir,
      "Figure_MetaCyc_pathway_PCoA_D0_D30_D60_D90_combined.png"
    ),

    plot =
      combined_plot,

    width =
      13,

    height =
      10,

    dpi =
      600
  )
}


# ======================================================
# 27. Finish
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "HUMAnN MetaCyc pathway PCoA analysis completed.\n"
)

cat(
  "Time points:\n"
)

cat(
  "  D0\n"
)

cat(
  "  D30\n"
)

cat(
  "  D60\n"
)

cat(
  "  D90\n"
)

cat(
  "Pairwise PERMANOVA:\n"
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
  "999 permutations.\n"
)

cat(
  "BH/FDR correction applied within each time point.\n"
)

cat(
  "Raw P values are displayed in the PCoA panels.\n"
)

cat(
  "Pathway input:\n"
)

cat(
  pathway_path,
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