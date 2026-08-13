########################################################
## HUMAnN MetaCyc pathway differential analysis
##
## Figure 5
##
## Post-baseline period:
##   D30 + D60 + D90
##
## Input:
##   metadata/metadata.csv
##
##   processed_data/HUMAnN/
##     pathabundance_relab_unstratified.tsv
##
## Pairwise comparisons:
##   OTM1 vs ITM
##   OTM2 vs ITM
##   OTM2 vs OTM1
##
## Statistical analysis:
##   Wilcoxon rank-sum test (unpaired)
##
## Multiple-testing correction:
##   Benjamini-Hochberg (BH/FDR)
##   applied separately within each pairwise comparison
##
## Pathway filtering:
##   Mean relative abundance >= 1e-6
##   Prevalence >= 5 post-baseline samples
##
## Fold change:
##   log2FC based on group mean relative abundance
##
## Volcano plot:
##   x-axis = log2 fold change
##   y-axis = -log10(FDR)
##   significance threshold = FDR < 0.05
##   no pathway text labels
##
## Output:
##   results/Fig5_HUMAnN_pathway_volcano/
########################################################


# ======================================================
# 1. Load packages
# ======================================================

library(tidyverse)
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
  "Fig5_HUMAnN_pathway_volcano"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ======================================================
# 3. Analysis settings
# ======================================================

use_times <- c(
  "D30",
  "D60",
  "D90"
)


group_levels <- c(
  "ITM",
  "OTM1",
  "OTM2"
)


comparisons <- list(

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


expected_times <- c(
  "D0",
  "D30",
  "D60",
  "D90"
)


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
  levels = group_levels
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
  "HUMAnN MetaCyc differential pathway analysis\n"
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
# 6. Read HUMAnN pathway table
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


# First column normally contains pathway IDs/names
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
# 8. Convert abundance columns to numeric
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

    .groups = "drop"
  )


cat(
  "\nUnique pathways:",
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
# Retained from the original workflow.
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


# Export name mapping
write.csv(

  tibble(

    Original_sample_name =
      original_sample_names,

    ID =
      clean_sample_names
  ),

  file.path(
    output_dir,
    "HUMAnN_sample_name_mapping.csv"
  ),

  row.names = FALSE
)


cat(
  "\nPathway table before sample matching:\n"
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
# 12. Match pathway data with metadata
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
    "\nSamples in metadata but missing from pathway table:\n"
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
    "\nSamples in pathway table but missing from metadata:\n"
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
# 13. Retain D30 + D60 + D90
# ======================================================

post_samples <- meta$ID[
  as.character(
    meta$Time
  ) %in%
    use_times
]


pathway_post <- pathway_mat[
  post_samples,
  ,
  drop = FALSE
]


meta_post <- meta[
  match(
    post_samples,
    meta$ID
  ),
  ,
  drop = FALSE
]


meta_post$Group <- droplevels(
  meta_post$Group
)


meta_post$Time <- droplevels(
  meta_post$Time
)


cat(
  "\n========================================\n"
)

cat(
  "Post-baseline samples\n"
)

cat(
  "========================================\n"
)


print(
  table(
    meta_post$Group,
    meta_post$Time
  )
)


cat(
  "\nTotal D30-D90 samples:",
  nrow(meta_post),
  "\n"
)


if (nrow(meta_post) != 90) {

  warning(
    paste0(
      "Expected 90 D30-D90 samples, but ",
      nrow(meta_post),
      " were retained."
    )
  )
}


# ======================================================
# 14. Remove all-zero pathways
# ======================================================

pathway_post <- pathway_post[
  ,
  colSums(
    pathway_post,
    na.rm = TRUE
  ) > 0,
  drop = FALSE
]


# ======================================================
# 15. Filter low-abundance / low-prevalence pathways
#
# Filtering is based on D30 + D60 + D90 samples,
# consistent with the original differential analysis.
#
# Criteria:
#   mean abundance >= 1e-6
#   prevalence >= 5 samples
# ======================================================

mean_abundance <- colMeans(
  pathway_post,
  na.rm = TRUE
)


prevalence <- colSums(
  pathway_post > 0,
  na.rm = TRUE
)


pathway_filter_table <- tibble(

  Pathway =
    colnames(
      pathway_post
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
    "stat_D30D60D90_MetaCyc_pathway_filtering.csv"
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


pathway_post_f <- pathway_post[
  ,
  keep_pathways,
  drop = FALSE
]


cat(
  "\nPathways retained after filtering:",
  ncol(pathway_post_f),
  "\n"
)


if (ncol(pathway_post_f) == 0) {

  stop(
    "No pathways remained after abundance/prevalence filtering."
  )
}


# ======================================================
# 16. Export filtered pathway matrix
# ======================================================

write.table(

  pathway_post_f,

  file.path(
    output_dir,
    "derived_MetaCyc_pathway_D30D60D90_filtered_matrix.tsv"
  ),

  sep = "\t",

  quote = FALSE,

  col.names = NA
)


# ======================================================
# 17. Parse pathway names
#
# Example:
# PWY-XXXX: pathway description
# ======================================================

parse_pathway_name <- function(x) {


  tibble(

    Pathway =
      x,

    Pathway_ID =
      ifelse(

        grepl(
          ":",
          x
        ),

        sub(
          ":.*$",
          "",
          x
        ),

        x
      ),

    Description =
      ifelse(

        grepl(
          ":",
          x
        ),

        sub(
          "^[^:]+:\\s*",
          "",
          x
        ),

        x
      )
  )
}


# ======================================================
# 18. Differential pathway test function
#
# group2 vs group1:
#
# Log2FC > 0:
#   higher in group2
#
# Log2FC < 0:
#   higher in group1
#
# D30, D60, and D90 are pooled.
#
# Wilcoxon rank-sum tests are intentionally unpaired,
# consistent with the original analysis.
#
# BH correction is performed separately within each
# pairwise comparison across all retained pathways.
# ======================================================

diff_pathway_test <- function(
    mat,
    metadata,
    group1,
    group2,
    pseudocount = 1e-10
) {


  samples_use <- metadata$ID[
    as.character(
      metadata$Group
    ) %in%
      c(
        group1,
        group2
      )
  ]


  mat_sub <- mat[
    samples_use,
    ,
    drop = FALSE
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


  n_group1 <- sum(
    meta_sub$Group == group1
  )


  n_group2 <- sum(
    meta_sub$Group == group2
  )


  cat(
    group1,
    ":",
    n_group1,
    "samples |",
    group2,
    ":",
    n_group2,
    "samples\n"
  )


  if (
    n_group1 < 2 ||
    n_group2 < 2
  ) {

    stop(
      paste0(
        "Insufficient samples for ",
        group2,
        " vs ",
        group1
      )
    )
  }


  pathways <- colnames(
    mat_sub
  )


  result_list <- lapply(

    pathways,

    function(pathway_name) {


      x1 <- mat_sub[
        meta_sub$Group == group1,
        pathway_name
      ]


      x2 <- mat_sub[
        meta_sub$Group == group2,
        pathway_name
      ]


      mean1 <- mean(
        x1,
        na.rm = TRUE
      )


      mean2 <- mean(
        x2,
        na.rm = TRUE
      )


      median1 <- median(
        x1,
        na.rm = TRUE
      )


      median2 <- median(
        x2,
        na.rm = TRUE
      )


      prevalence1 <- sum(
        x1 > 0,
        na.rm = TRUE
      )


      prevalence2 <- sum(
        x2 > 0,
        na.rm = TRUE
      )


      p_value <- tryCatch(

        wilcox.test(

          x2,

          x1,

          paired = FALSE,

          exact = FALSE

        )$p.value,

        error = function(e) {

          NA_real_
        }
      )


      tibble(

        Pathway =
          pathway_name,

        Group_ref =
          group1,

        Group_test =
          group2,

        Comparison =
          paste(
            group2,
            "vs",
            group1
          ),

        N_ref =
          length(
            x1
          ),

        N_test =
          length(
            x2
          ),

        Mean_ref =
          mean1,

        Mean_test =
          mean2,

        Median_ref =
          median1,

        Median_test =
          median2,

        Prevalence_ref =
          prevalence1,

        Prevalence_test =
          prevalence2,

        Mean_difference =
          mean2 -
          mean1,

        Log2FC =
          log2(
            (
              mean2 +
                pseudocount
            ) /
              (
                mean1 +
                  pseudocount
              )
          ),

        P_value =
          p_value
      )
    }
  )


  result_df <- bind_rows(
    result_list
  ) %>%

    mutate(

      # BH correction within this comparison
      FDR =
        p.adjust(
          P_value,
          method = "BH"
        )
    )


  # Add pathway ID / description
  result_df <- result_df %>%

    left_join(

      parse_pathway_name(
        result_df$Pathway
      ),

      by = "Pathway"
    ) %>%

    mutate(

      Direction = case_when(

        FDR < 0.05 &
          Log2FC > 0 ~

          paste0(
            "Enriched in ",
            group2
          ),

        FDR < 0.05 &
          Log2FC < 0 ~

          paste0(
            "Enriched in ",
            group1
          ),

        TRUE ~

          "Not significant"
      ),


      NegLog10P =
        -log10(
          pmax(
            P_value,
            1e-300
          )
        ),


      NegLog10FDR =
        -log10(
          pmax(
            FDR,
            1e-300
          )
        )
    ) %>%

    arrange(
      FDR,
      desc(
        abs(
          Log2FC
        )
      )
    )


  return(
    result_df
  )
}


# ======================================================
# 19. Run pairwise comparisons
# ======================================================

diff_results <- list()


for (
  comparison in comparisons
) {


  reference_group <- comparison[1]

  test_group <- comparison[2]


  comparison_name <- paste(
    test_group,
    "vs",
    reference_group
  )


  cat(
    "\n========================================\n"
  )


  cat(
    "Running comparison:",
    comparison_name,
    "\n"
  )


  cat(
    "========================================\n"
  )


  result <- diff_pathway_test(

    mat =
      pathway_post_f,

    metadata =
      meta_post,

    group1 =
      reference_group,

    group2 =
      test_group
  )


  diff_results[
    [
      comparison_name
    ]
  ] <- result


  write.csv(

    result,

    file.path(
      output_dir,
      paste0(
        "stat_MetaCyc_pathways_",
        test_group,
        "_vs_",
        reference_group,
        "_D30D60D90.csv"
      )
    ),

    row.names = FALSE
  )


  cat(
    "FDR < 0.05:",
    sum(
      result$FDR < 0.05,
      na.rm = TRUE
    ),
    "\n"
  )


  cat(
    "Raw P < 0.05:",
    sum(
      result$P_value < 0.05,
      na.rm = TRUE
    ),
    "\n"
  )
}


# ======================================================
# 20. Combine all differential results
# ======================================================

all_diff <- bind_rows(
  diff_results
)


write.csv(

  all_diff,

  file.path(
    output_dir,
    "stat_MetaCyc_pathways_all_comparisons_D30D60D90.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 21. Export FDR-significant pathways
# ======================================================

significant_diff <- all_diff %>%

  filter(
    FDR < 0.05
  ) %>%

  arrange(
    Comparison,
    FDR,
    desc(
      abs(
        Log2FC
      )
    )
  )


write.csv(

  significant_diff,

  file.path(
    output_dir,
    "stat_MetaCyc_pathways_FDR_lt_0.05_D30D60D90.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 22. Differential pathway summary
# ======================================================

summary_table <- all_diff %>%

  group_by(
    Comparison
  ) %>%

  summarise(

    Total_pathways_tested =
      n(),

    Significant_FDR_0.05 =
      sum(
        FDR < 0.05,
        na.rm = TRUE
      ),

    Nominal_P_0.05 =
      sum(
        P_value < 0.05,
        na.rm = TRUE
      ),

    Enriched_in_ref_FDR =
      sum(
        FDR < 0.05 &
          Log2FC < 0,
        na.rm = TRUE
      ),

    Enriched_in_test_FDR =
      sum(
        FDR < 0.05 &
          Log2FC > 0,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  )


write.csv(

  summary_table,

  file.path(
    output_dir,
    "stat_MetaCyc_pathway_differential_summary_D30D60D90.csv"
  ),

  row.names = FALSE
)


cat(
  "\n========================================\n"
)

cat(
  "Differential pathway summary\n"
)

cat(
  "========================================\n"
)


print(
  summary_table
)


# ======================================================
# 23. Volcano colors
# ======================================================

volcano_colors <- c(

  "Enriched in ITM" =
    "#8ECFC9",

  "Enriched in OTM1" =
    "#FFBE7A",

  "Enriched in OTM2" =
    "#FA7F6F",

  "Not significant" =
    "grey80"
)


# ======================================================
# 24. Volcano plot function
#
# Statistical significance, point color, y-axis, and
# horizontal threshold all consistently use FDR.
# ======================================================

plot_volcano <- function(
    result,
    comparison_name,
    xlim_range = c(
      -4.5,
      4.5
    )
) {


  result_plot <- result %>%

    mutate(

      Volcano_group = factor(

        Direction,

        levels = c(
          "Enriched in ITM",
          "Enriched in OTM1",
          "Enriched in OTM2",
          "Not significant"
        )
      )
    )


  p <- ggplot(

    result_plot,

    aes(
      x = Log2FC,
      y = NegLog10FDR
    )

  ) +


    geom_point(

      aes(
        color =
          Volcano_group
      ),

      alpha =
        0.80,

      size =
        1.8
    ) +


    # log2FC = 0
    geom_vline(

      xintercept =
        0,

      linetype =
        2,

      linewidth =
        0.4
    ) +


    # FDR = 0.05
    geom_hline(

      yintercept =
        -log10(
          0.05
        ),

      linetype =
        2,

      linewidth =
        0.4
    ) +


    scale_color_manual(

      values =
        volcano_colors,

      drop =
        FALSE
    ) +


    coord_cartesian(

      xlim =
        xlim_range
    ) +


    labs(

      title =
        comparison_name,

      x =
        "log2 fold change",

      y =
        "\u2212log10(FDR)",

      color =
        NULL
    ) +


    theme_bw(
      base_size = 12
    ) +


    theme(

      panel.grid =
        element_blank(),


      axis.text =
        element_text(
          color = "black"
        ),


      axis.title =
        element_text(
          face = "bold"
        ),


      plot.title =
        element_text(
          face = "bold",
          size = 12,
          hjust = 0.5
        ),


      legend.title =
        element_blank(),


      legend.position =
        "bottom"
    )


  return(
    p
  )
}


# ======================================================
# 25. Generate volcano plots
# ======================================================

volcano_plots <- list()


volcano_plots[
  [
    "OTM1 vs ITM"
  ]
] <- plot_volcano(

  diff_results[
    [
      "OTM1 vs ITM"
    ]
  ],

  "OTM1 vs ITM"
)


volcano_plots[
  [
    "OTM2 vs ITM"
  ]
] <- plot_volcano(

  diff_results[
    [
      "OTM2 vs ITM"
    ]
  ],

  "OTM2 vs ITM"
)


volcano_plots[
  [
    "OTM2 vs OTM1"
  ]
] <- plot_volcano(

  diff_results[
    [
      "OTM2 vs OTM1"
    ]
  ],

  "OTM2 vs OTM1"
)


# ======================================================
# 26. Combine three panels in one row
# ======================================================

combined_volcano <- (

  volcano_plots[
    [
      "OTM1 vs ITM"
    ]
  ] |

  volcano_plots[
    [
      "OTM2 vs ITM"
    ]
  ] |

  volcano_plots[
    [
      "OTM2 vs OTM1"
    ]
  ]

) +


  patchwork::plot_layout(

    guides =
      "collect"
  ) &


  theme(

    legend.position =
      "bottom"
  )


print(
  combined_volcano
)


# ======================================================
# 27. Save combined volcano figure
# ======================================================

ggsave(

  filename = file.path(
    output_dir,
    "Figure_MetaCyc_pathway_D30D60D90_volcano_combined.pdf"
  ),

  plot =
    combined_volcano,

  width =
    16,

  height =
    5.5
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_MetaCyc_pathway_D30D60D90_volcano_combined.png"
  ),

  plot =
    combined_volcano,

  width =
    16,

  height =
    5.5,

  dpi =
    600
)


# ======================================================
# 28. Save individual volcano plots
# ======================================================

for (
  comparison_name in names(
    volcano_plots
  )
) {


  safe_name <- gsub(
    " ",
    "_",
    comparison_name
  )


  ggsave(

    filename = file.path(
      output_dir,
      paste0(
        "Figure_MetaCyc_volcano_",
        safe_name,
        "_D30D60D90.pdf"
      )
    ),

    plot =
      volcano_plots[
        [
          comparison_name
        ]
      ],

    width =
      6,

    height =
      5
  )


  ggsave(

    filename = file.path(
      output_dir,
      paste0(
        "Figure_MetaCyc_volcano_",
        safe_name,
        "_D30D60D90.png"
      )
    ),

    plot =
      volcano_plots[
        [
          comparison_name
        ]
      ],

    width =
      6,

    height =
      5,

    dpi =
      600
  )
}


# ======================================================
# 29. Finish
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "HUMAnN MetaCyc differential pathway analysis completed.\n"
)

cat(
  "Samples used: D30 + D60 + D90.\n"
)

cat(
  "Pairwise comparisons:\n"
)

cat(
  "  OTM1 vs ITM\n"
)

cat(
  "  OTM2 vs ITM\n"
)

cat(
  "  OTM2 vs OTM1\n"
)

cat(
  "Wilcoxon rank-sum test (unpaired).\n"
)

cat(
  "BH/FDR correction applied separately within each comparison.\n"
)

cat(
  "Volcano y-axis and significance threshold use FDR.\n"
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