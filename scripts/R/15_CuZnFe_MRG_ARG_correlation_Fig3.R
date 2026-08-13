########################################################
## Group-specific associations between Cu/Zn/Fe-related
## MRG abundance and ARG type abundance
##
## Figure 3
##
## Input:
##   metadata/metadata.csv
##
##   processed_data/ARG/
##     normalized_cell.type_matrix.tsv
##
##   processed_data/MRG/
##     MRG_merged_subtype_abundance.tsv
##
## Data processing:
##   MRG merged subtypes are aggregated to their
##   corresponding metal-resistance type/class.
##
##   Example:
##     Copper__subtype1
##     Copper__subtype2
##          ↓
##        Copper
##
## Analysis:
##   Post-baseline samples only:
##   D30, D60, and D90
##
##   Spearman correlations are calculated separately
##   within ITM, OTM1, and OTM2.
##
## Metal categories:
##   Copper
##   Zinc
##   Iron
##
## Multiple-testing correction:
##   Benjamini-Hochberg (BH/FDR)
##   applied within each Group × Metal combination
##
## Plot:
##   x-axis = Copper / Zinc / Iron
##   y-axis = selected ARG types
##   facets = ITM / OTM1 / OTM2
##   color = Spearman rho
##   size = -log10(FDR)
##   stars = BH/FDR significance
##
## Top ARG types:
##   Top 10 within each Group × Metal combination
##   ranked by FDR, then absolute Spearman rho
##
## Output:
##   results/Fig3_MRG_ARG_correlation/
########################################################


# ======================================================
# 1. Load packages
# ======================================================

library(tidyverse)
library(stringr)


# ======================================================
# 2. Project paths
#
# Run this script from the repository root:
# DairyCow-TraceMineral-Resistome/
# ======================================================

PROJECT_DIR <- "."


meta_path <- file.path(
  PROJECT_DIR,
  "metadata",
  "metadata.csv"
)


arg_type_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "ARG",
  "normalized_cell.type_matrix.tsv"
)


mrg_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "MRG",
  "MRG_merged_subtype_abundance.tsv"
)


output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig3_MRG_ARG_correlation"
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


# Number of ARG types selected within each
# Group × Metal combination
top_n <- 10


# Ranking criterion
rank_by <- "FDR"

# Alternative:
# rank_by <- "rawP"


# ======================================================
# 3. Read metadata
# ======================================================

meta_df <- read.csv(
  meta_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


required_cols <- c(
  "ID",
  "Group",
  "Time"
)


missing_cols <- setdiff(
  required_cols,
  colnames(meta_df)
)


if (length(missing_cols) > 0) {

  stop(
    paste0(
      "Missing metadata columns: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  )
}


# ======================================================
# 4. Prepare metadata
# ======================================================

meta_df <- meta_df %>%

  mutate(

    ID = as.character(ID),

    Group = as.character(Group),

    Time = as.character(Time)
  )


expected_groups <- c(
  "ITM",
  "OTM1",
  "OTM2"
)


unexpected_groups <- setdiff(
  unique(meta_df$Group),
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


meta_df$Group <- factor(
  meta_df$Group,
  levels = c(
    "ITM",
    "OTM1",
    "OTM2"
  )
)


meta_df$Time <- factor(
  meta_df$Time,
  levels = c(
    "D0",
    "D30",
    "D60",
    "D90"
  )
)


# Keep post-baseline samples
meta_df <- meta_df %>%

  filter(
    Time %in% use_times
  ) %>%

  droplevels()


if (anyDuplicated(meta_df$ID) > 0) {

  stop(
    "Duplicated sequencing sample IDs were found."
  )
}


# ======================================================
# 5. Metadata checks
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "MRG-ARG group-specific correlation analysis\n"
)

cat(
  "========================================\n"
)


cat(
  "Sampling times included:",
  paste(
    use_times,
    collapse = " + "
  ),
  "\n"
)


cat(
  "\nSamples by Group × Time:\n"
)


print(
  table(
    meta_df$Group,
    meta_df$Time
  )
)


cat(
  "\nTotal post-baseline samples:",
  nrow(meta_df),
  "\n"
)


# ======================================================
# 6. Read ARG type matrix
# ======================================================

arg_df <- read.delim(
  arg_type_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ======================================================
# 7. Read MRG merged-subtype matrix
# ======================================================

mrg_df <- read.delim(
  mrg_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ======================================================
# 8. Function to construct sample × feature matrix
# ======================================================

build_sample_matrix <- function(
    feature_df,
    matrix_name
) {


  colnames(feature_df)[1] <- "feature"


  feature_df$feature <- as.character(
    feature_df$feature
  )


  if (anyDuplicated(feature_df$feature) > 0) {

    warning(
      paste0(
        "Duplicated feature names detected in ",
        matrix_name,
        ". make.unique() will be applied."
      )
    )
  }


  x <- feature_df[
    ,
    -1,
    drop = FALSE
  ]


  x[] <- lapply(

    x,

    function(v) {

      as.numeric(
        as.character(v)
      )
    }
  )


  mat <- t(
    as.matrix(x)
  )


  rownames(mat) <- colnames(
    feature_df
  )[-1]


  colnames(mat) <- make.unique(
    feature_df$feature
  )


  # Missing values -> zero
  mat[
    is.na(mat)
  ] <- 0


  # Negative values -> zero
  mat[
    mat < 0
  ] <- 0


  # Remove globally absent features
  mat <- mat[
    ,
    colSums(
      mat,
      na.rm = TRUE
    ) > 0,
    drop = FALSE
  ]


  return(
    mat
  )
}


# ======================================================
# 9. Build ARG type matrix
# ======================================================

arg_mat <- build_sample_matrix(
  feature_df = arg_df,
  matrix_name = "ARG type matrix"
)


cat(
  "\nARG type matrix: ",
  nrow(arg_mat),
  " samples × ",
  ncol(arg_mat),
  " ARG types\n",
  sep = ""
)


# ======================================================
# 10. Prepare MRG merged-subtype abundance matrix
#
# MRG merged subtype names are expected to contain
# their resistance type before "__".
#
# Examples:
#
# Copper__xxx -> Copper
# Zinc__xxx   -> Zinc
# Iron__xxx   -> Iron
# ======================================================

colnames(mrg_df)[1] <- "feature"


mrg_df$feature <- as.character(
  mrg_df$feature
)


mrg_sample_cols <- colnames(
  mrg_df
)[-1]


# Convert abundance columns to numeric
mrg_df[
  mrg_sample_cols
] <- lapply(

  mrg_df[
    mrg_sample_cols
  ],

  function(v) {

    v <- as.numeric(
      as.character(v)
    )


    v[
      is.na(v)
    ] <- 0


    v[
      v < 0
    ] <- 0


    return(
      v
    )
  }
)


# ======================================================
# 11. Extract MRG type/class
# ======================================================

mrg_df <- mrg_df %>%

  mutate(

    MRG_type = if_else(

      str_detect(
        feature,
        fixed("__")
      ),

      str_replace(
        feature,
        "__.*$",
        ""
      ),

      feature
    )
  )


cat(
  "\n========================================\n"
)

cat(
  "MRG merged-subtype -> type aggregation\n"
)

cat(
  "========================================\n"
)


cat(
  "Number of MRG merged subtypes:",
  nrow(mrg_df),
  "\n"
)


cat(
  "Number of MRG types/classes:",
  length(
    unique(
      mrg_df$MRG_type
    )
  ),
  "\n"
)


cat(
  "\nMRG types detected:\n"
)


print(
  sort(
    unique(
      mrg_df$MRG_type
    )
  )
)


# ======================================================
# 12. Export MRG subtype-to-type mapping
# ======================================================

mrg_mapping_df <- mrg_df %>%

  select(
    feature,
    MRG_type
  ) %>%

  distinct()


write.csv(

  mrg_mapping_df,

  file.path(
    output_dir,
    "MRG_merged_subtype_to_type_mapping.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 13. Aggregate MRG merged subtypes to MRG type
# ======================================================

mrg_type_df <- mrg_df %>%

  select(
    MRG_type,
    all_of(
      mrg_sample_cols
    )
  ) %>%

  group_by(
    MRG_type
  ) %>%

  summarise(

    across(
      all_of(
        mrg_sample_cols
      ),
      ~ sum(
        .x,
        na.rm = TRUE
      )
    ),

    .groups = "drop"
  )


# Export derived type matrix
write.table(

  mrg_type_df,

  file.path(
    output_dir,
    "derived_MRG_type_abundance.tsv"
  ),

  sep = "\t",

  quote = FALSE,

  row.names = FALSE
)


# ======================================================
# 14. Construct sample × MRG type matrix
# ======================================================

mrg_x <- mrg_type_df[
  ,
  -1,
  drop = FALSE
]


mrg_mat <- t(
  as.matrix(
    mrg_x
  )
)


rownames(mrg_mat) <- colnames(
  mrg_type_df
)[-1]


colnames(mrg_mat) <- make.unique(
  as.character(
    mrg_type_df$MRG_type
  )
)


mrg_mat[
  is.na(mrg_mat)
] <- 0


mrg_mat[
  mrg_mat < 0
] <- 0


cat(
  "\nMRG type matrix: ",
  nrow(mrg_mat),
  " samples × ",
  ncol(mrg_mat),
  " MRG types\n",
  sep = ""
)


# ======================================================
# 15. Clean ARG and MRG type names
# ======================================================

clean_type_name <- function(x) {

  x %>%

    as.character() %>%

    str_replace_all(
      "_",
      " "
    ) %>%

    str_replace_all(
      "-",
      " "
    ) %>%

    str_replace_all(
      "\\s+",
      " "
    ) %>%

    str_trim()
}


colnames(arg_mat) <- clean_type_name(
  colnames(arg_mat)
)


colnames(mrg_mat) <- clean_type_name(
  colnames(mrg_mat)
)


colnames(arg_mat) <- make.unique(
  colnames(arg_mat)
)


colnames(mrg_mat) <- make.unique(
  colnames(mrg_mat)
)


# ======================================================
# 16. Check required metal categories
# ======================================================

required_metals <- c(
  "Copper",
  "Zinc",
  "Iron"
)


missing_metals <- setdiff(
  required_metals,
  colnames(mrg_mat)
)


if (length(missing_metals) > 0) {

  stop(
    paste0(
      "Required MRG type(s) not found after aggregation: ",
      paste(
        missing_metals,
        collapse = ", "
      )
    )
  )
}


cat(
  "\nMRG metal categories used:\n"
)


print(
  required_metals
)


# ======================================================
# 17. Match ARG, MRG, and metadata samples
# ======================================================

common_samples <- Reduce(

  intersect,

  list(

    rownames(
      arg_mat
    ),

    rownames(
      mrg_mat
    ),

    meta_df$ID
  )
)


cat(
  "\nMatched D30-D90 samples:",
  length(common_samples),
  "\n"
)


if (length(common_samples) != 90) {

  warning(
    paste0(
      "Expected 90 post-baseline samples, but ",
      length(common_samples),
      " were matched."
    )
  )
}


if (length(common_samples) < 5) {

  stop(
    paste0(
      "Too few common samples. ",
      "Please check ARG, MRG, and metadata sample IDs."
    )
  )
}


arg_mat <- arg_mat[
  common_samples,
  ,
  drop = FALSE
]


mrg_mat <- mrg_mat[
  common_samples,
  ,
  drop = FALSE
]


meta_df <- meta_df[
  match(
    common_samples,
    meta_df$ID
  ),
  ,
  drop = FALSE
]


stopifnot(
  all(
    rownames(arg_mat) ==
      rownames(mrg_mat)
  )
)


stopifnot(
  all(
    rownames(arg_mat) ==
      meta_df$ID
  )
)


cat(
  "\nSamples after alignment:\n"
)


print(
  table(
    meta_df$Group,
    meta_df$Time
  )
)


# ======================================================
# 18. Build Copper / Zinc / Iron MRG abundance table
# ======================================================

metal_df <- tibble(

  ID =
    common_samples,

  Copper =
    as.numeric(
      mrg_mat[
        ,
        "Copper"
      ]
    ),

  Zinc =
    as.numeric(
      mrg_mat[
        ,
        "Zinc"
      ]
    ),

  Iron =
    as.numeric(
      mrg_mat[
        ,
        "Iron"
      ]
    )

) %>%

  left_join(

    meta_df %>%

      select(
        ID,
        Group,
        Time
      ),

    by = "ID"
  )


# Export exact metal abundances used
write.csv(

  metal_df,

  file.path(
    output_dir,
    "D30D60D90_CuZnFe_MRG_abundance_used.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 19. Convert ARG type matrix to long format
# ======================================================

arg_long <- as.data.frame(
  arg_mat
) %>%

  tibble::rownames_to_column(
    "ID"
  ) %>%

  pivot_longer(

    cols = -ID,

    names_to =
      "ARG_type",

    values_to =
      "ARG_abundance"
  )


# ======================================================
# 20. Build correlation input table
# ======================================================

corr_input <- arg_long %>%

  left_join(
    metal_df,
    by = "ID"
  ) %>%

  filter(
    !is.na(Group),
    !is.na(Time)
  ) %>%

  pivot_longer(

    cols = c(
      Copper,
      Zinc,
      Iron
    ),

    names_to =
      "Metal",

    values_to =
      "MRG_abundance"
  ) %>%

  mutate(

    Group = factor(
      Group,
      levels = c(
        "ITM",
        "OTM1",
        "OTM2"
      )
    ),

    Metal = factor(
      Metal,
      levels = c(
        "Copper",
        "Zinc",
        "Iron"
      )
    ),

    Time = factor(
      Time,
      levels = c(
        "D30",
        "D60",
        "D90"
      )
    )
  )


# ======================================================
# 21. Safe Spearman correlation function
# ======================================================

safe_spearman <- function(
    x,
    y
) {


  keep <- complete.cases(
    x,
    y
  )


  x <- x[
    keep
  ]


  y <- y[
    keep
  ]


  if (length(x) < 3) {

    return(
      tibble(
        rho = NA_real_,
        p_value = NA_real_
      )
    )
  }


  if (
    sd(
      x,
      na.rm = TRUE
    ) == 0 ||

    sd(
      y,
      na.rm = TRUE
    ) == 0
  ) {

    return(
      tibble(
        rho = NA_real_,
        p_value = NA_real_
      )
    )
  }


  test_res <- suppressWarnings(

    cor.test(

      x,

      y,

      method =
        "spearman",

      exact =
        FALSE
    )
  )


  tibble(

    rho =
      unname(
        test_res$estimate
      ),

    p_value =
      test_res$p.value
  )
}


# ======================================================
# 22. Group-specific Spearman correlations
#
# D30, D60, and D90 samples are pooled within
# each treatment group.
# ======================================================

cor_df <- corr_input %>%

  group_by(
    Group,
    Metal,
    ARG_type
  ) %>%

  group_modify(

    ~ {

      cor_res <- safe_spearman(
        .x$MRG_abundance,
        .x$ARG_abundance
      )


      tibble(

        n =
          sum(
            complete.cases(
              .x$MRG_abundance,
              .x$ARG_abundance
            )
          ),

        rho =
          cor_res$rho,

        p_value =
          cor_res$p_value,

        mean_ARG =
          mean(
            .x$ARG_abundance,
            na.rm = TRUE
          ),

        detection_rate_ARG =
          mean(
            .x$ARG_abundance > 0,
            na.rm = TRUE
          )
      )
    }
  ) %>%

  ungroup() %>%

  filter(
    !is.na(rho),
    !is.na(p_value)
  )


# ======================================================
# 23. BH/FDR correction
#
# Separately within each Group × Metal combination.
# ======================================================

cor_df <- cor_df %>%

  group_by(
    Group,
    Metal
  ) %>%

  mutate(

    p_adj =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%

  ungroup() %>%

  mutate(


    signif_label = case_when(

      p_adj < 0.001 ~
        "***",

      p_adj < 0.01 ~
        "**",

      p_adj < 0.05 ~
        "*",

      TRUE ~
        ""
    ),


    raw_signif_label = case_when(

      p_value < 0.001 ~
        "***",

      p_value < 0.01 ~
        "**",

      p_value < 0.05 ~
        "*",

      TRUE ~
        ""
    ),


    neg_log10_FDR =
      -log10(
        pmax(
          p_adj,
          1e-300
        )
      ),


    neg_log10_P =
      -log10(
        pmax(
          p_value,
          1e-300
        )
      )
  )


# ======================================================
# 24. Export complete correlation results
# ======================================================

write.csv(

  cor_df %>%

    arrange(
      Group,
      Metal,
      p_adj
    ),

  file.path(
    output_dir,
    "D30D60D90_stat_CuZnFe_vs_ARGtype_Spearman_byGroup_all.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 25. Export significant results
# ======================================================

write.csv(

  cor_df %>%

    filter(
      p_adj < 0.05
    ) %>%

    arrange(
      Group,
      Metal,
      p_adj
    ),

  file.path(
    output_dir,
    "D30D60D90_stat_CuZnFe_vs_ARGtype_Spearman_byGroup_FDR_significant.csv"
  ),

  row.names = FALSE
)


write.csv(

  cor_df %>%

    filter(
      p_value < 0.05
    ) %>%

    arrange(
      Group,
      Metal,
      p_value
    ),

  file.path(
    output_dir,
    "D30D60D90_stat_CuZnFe_vs_ARGtype_Spearman_byGroup_rawP_significant.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 26. Print significance summaries
# ======================================================

cat(
  "\nFDR < 0.05 correlation counts:\n"
)


print(

  cor_df %>%

    filter(
      p_adj < 0.05
    ) %>%

    count(
      Group,
      Metal,
      name = "n_FDR_lt_0.05"
    )
)


cat(
  "\nRaw P < 0.05 correlation counts:\n"
)


print(

  cor_df %>%

    filter(
      p_value < 0.05
    ) %>%

    count(
      Group,
      Metal,
      name = "n_rawP_lt_0.05"
    )
)


# ======================================================
# 27. Select Top 10 ARG types
# ======================================================

if (rank_by == "FDR") {


  plot_select_df <- cor_df %>%

    group_by(
      Group,
      Metal
    ) %>%

    arrange(
      p_adj,
      desc(
        abs(
          rho
        )
      ),
      .by_group = TRUE
    ) %>%

    slice_head(
      n = top_n
    ) %>%

    ungroup()


} else {


  plot_select_df <- cor_df %>%

    group_by(
      Group,
      Metal
    ) %>%

    arrange(
      p_value,
      desc(
        abs(
          rho
        )
      ),
      .by_group = TRUE
    ) %>%

    slice_head(
      n = top_n
    ) %>%

    ungroup()
}


# Union of selected ARG types across all
# Group × Metal combinations
selected_arg_types <- unique(
  plot_select_df$ARG_type
)


plot_df <- cor_df %>%

  filter(
    ARG_type %in%
      selected_arg_types
  )


# ======================================================
# 28. Order ARG types
# ======================================================

arg_order <- plot_df %>%

  group_by(
    ARG_type
  ) %>%

  summarise(

    max_abs_rho =
      max(
        abs(
          rho
        ),
        na.rm = TRUE
      ),

    min_p_adj =
      min(
        p_adj,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  ) %>%

  arrange(
    desc(
      max_abs_rho
    ),
    min_p_adj
  ) %>%

  pull(
    ARG_type
  )


plot_df <- plot_df %>%

  mutate(

    Group = factor(
      Group,
      levels = c(
        "ITM",
        "OTM1",
        "OTM2"
      )
    ),

    Metal = factor(
      Metal,
      levels = c(
        "Copper",
        "Zinc",
        "Iron"
      )
    ),

    ARG_type = factor(
      ARG_type,
      levels = rev(
        arg_order
      )
    )
  )


# ======================================================
# 29. Export plotting datasets
# ======================================================

write.csv(

  plot_df,

  file.path(
    output_dir,
    "D30D60D90_stat_CuZnFe_vs_ARGtype_plot_data_top10.csv"
  ),

  row.names = FALSE
)


write.csv(

  plot_select_df,

  file.path(
    output_dir,
    "D30D60D90_stat_CuZnFe_vs_ARGtype_top10_selection.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 30. Dot heatmap
# ======================================================

p_dot_D30D90_top10 <- ggplot(

  plot_df,

  aes(
    x = Metal,
    y = ARG_type
  )

) +


  geom_point(

    aes(

      color =
        rho,

      size =
        neg_log10_FDR
    ),

    alpha =
      0.90
  ) +


  geom_text(

    aes(
      label =
        signif_label
    ),

    color =
      "black",

    size =
      3.0,

    fontface =
      "bold"
  ) +


  facet_wrap(

    ~ Group,

    nrow =
      1
  ) +


  scale_color_gradient2(

    low =
      "#2166AC",

    mid =
      "white",

    high =
      "#B2182B",

    midpoint =
      0,

    limits =
      c(
        -1,
        1
      ),

    name =
      "Spearman \u03c1"
  ) +


  scale_size_continuous(

    range =
      c(
        1.8,
        6
      ),

    name =
      "\u2212log10(FDR)"
  ) +


  labs(

    x =
      "MRG category",

    y =
      "ARG type",

    title =
      paste0(
        "D30-D90 group-specific associations between ",
        "Cu/Zn/Fe-related MRGs and ARG types"
      ),

    subtitle =
      paste0(
        "Top 10 ARG types per group and metal; ",
        "stars indicate BH/FDR significance"
      )
  ) +


  theme_bw(
    base_size = 11
  ) +


  theme(

    panel.grid.major =
      element_line(
        linewidth = 0.25,
        color = "grey88"
      ),

    panel.grid.minor =
      element_blank(),


    strip.background =
      element_rect(
        fill = "white",
        color = "black",
        linewidth = 0.6
      ),

    strip.text =
      element_text(
        face = "bold",
        size = 12
      ),


    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5,
        size = 14
      ),

    plot.subtitle =
      element_text(
        hjust = 0.5,
        size = 9.5
      ),


    axis.title =
      element_text(
        face = "bold"
      ),

    axis.text.x =
      element_text(
        face = "bold",
        color = "black",
        size = 10
      ),

    axis.text.y =
      element_text(
        color = "black",
        size = 8.5
      ),


    legend.position =
      "right",

    legend.title =
      element_text(
        face = "bold"
      ),


    plot.margin =
      ggplot2::margin(
        10,
        20,
        10,
        10
      )
  )


print(
  p_dot_D30D90_top10
)


# ======================================================
# 31. Save figure
# ======================================================

figure_height <- max(
  6,
  length(
    selected_arg_types
  ) * 0.22
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_D30D60D90_CuZnFe_vs_ARGtype_Spearman_dotplot.pdf"
  ),

  plot =
    p_dot_D30D90_top10,

  width =
    11,

  height =
    figure_height
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_D30D60D90_CuZnFe_vs_ARGtype_Spearman_dotplot.png"
  ),

  plot =
    p_dot_D30D90_top10,

  width =
    11,

  height =
    figure_height,

  dpi =
    600
)


# ======================================================
# 32. Finish
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "D30-D90 MRG-ARG correlation analysis completed.\n"
)

cat(
  "ARG input:\n"
)

cat(
  arg_type_path,
  "\n"
)

cat(
  "MRG core input:\n"
)

cat(
  mrg_path,
  "\n"
)

cat(
  "Top n = ",
  top_n,
  "\n",
  sep = ""
)

cat(
  "Ranking criterion = ",
  rank_by,
  "\n",
  sep = ""
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