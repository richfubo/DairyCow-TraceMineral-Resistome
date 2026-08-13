########################################################
## Group-specific associations between Cu/Zn/Fe-related
## MRG abundance and MGE type abundance
##
## Figure 3
##
## Input:
##   metadata/metadata.csv
##
##   processed_data/MRG/
##     MRG_merged_subtype_abundance.tsv
##
##   processed_data/MGE/
##     MGE_merged_subtype_abundance.tsv
##
## Data processing:
##
##   MRG merged subtypes are aggregated to MRG
##   type/class according to the prefix before "__".
##
##   Example:
##     Copper__subtype1
##     Copper__subtype2
##          ↓
##        Copper
##
##   MGE merged subtypes are aggregated to MGE
##   type/class according to the prefix before "__".
##
##   Example:
##     integrase__intI1
##     integrase__int2
##          ↓
##       integrase
##
## Analysis:
##   Post-baseline samples only:
##   D30, D60, and D90
##
## Correlation:
##   Spearman correlation between Cu/Zn/Fe-related
##   MRG abundance and individual MGE type abundance.
##
##   Correlations are calculated separately within:
##   ITM, OTM1, and OTM2.
##
## Multiple-testing correction:
##   Benjamini-Hochberg (BH/FDR)
##   within each Group × Metal combination.
##
## Plot:
##   x-axis = Copper / Zinc / Iron
##   y-axis = selected MGE types
##   facets = ITM / OTM1 / OTM2
##   color = Spearman rho
##   size = -log10(FDR)
##   stars = BH/FDR significance
##
## Top features:
##   Top 10 MGE types within each Group × Metal
##   combination ranked by FDR and |rho|.
##
## Output:
##   results/Fig3_MRG_MGE_correlation/
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


mrg_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "MRG",
  "MRG_merged_subtype_abundance.tsv"
)


mge_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "MGE",
  "MGE_merged_subtype_abundance.tsv"
)


output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig3_MRG_MGE_correlation"
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


# Number of MGE types selected within each
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


# Retain D30-D90
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
  "MRG-MGE group-specific correlation analysis\n"
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
# 6. Read MRG and MGE merged-subtype matrices
# ======================================================

mrg_df <- read.delim(
  mrg_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


mge_df <- read.delim(
  mge_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ======================================================
# 7. Function to aggregate merged subtype -> type
#
# Expected input:
#   first column = merged subtype
#   remaining columns = sample IDs
#
# Rule:
#   prefix before "__" is retained as the type/class.
#
# Example:
#   Copper__xxx -> Copper
#   integrase__intI1 -> integrase
# ======================================================

aggregate_merged_subtype <- function(
    feature_df,
    type_name
) {


  colnames(feature_df)[1] <- "feature"


  feature_df$feature <- as.character(
    feature_df$feature
  )


  sample_cols <- colnames(
    feature_df
  )[-1]


  # Convert abundance values to numeric
  feature_df[
    sample_cols
  ] <- lapply(

    feature_df[
      sample_cols
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


  # Extract type/class
  feature_df <- feature_df %>%

    mutate(

      type = if_else(

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


  # Mapping table
  mapping_df <- feature_df %>%

    select(
      feature,
      type
    ) %>%

    distinct()


  # Aggregate subtype abundance to type
  type_df <- feature_df %>%

    select(
      type,
      all_of(
        sample_cols
      )
    ) %>%

    group_by(
      type
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


  colnames(type_df)[1] <- type_name


  return(

    list(

      abundance =
        type_df,

      mapping =
        mapping_df
    )
  )
}


# ======================================================
# 8. Aggregate MRG merged subtypes -> MRG types
# ======================================================

mrg_agg <- aggregate_merged_subtype(
  feature_df = mrg_df,
  type_name = "MRG_type"
)


mrg_type_df <- mrg_agg$abundance


mrg_mapping_df <- mrg_agg$mapping %>%

  rename(
    MRG_type = type
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
  "Merged subtypes:",
  nrow(mrg_df),
  "\n"
)


cat(
  "MRG types/classes:",
  nrow(mrg_type_df),
  "\n"
)


cat(
  "\nMRG types detected:\n"
)


print(
  mrg_type_df$MRG_type
)


# ======================================================
# 9. Aggregate MGE merged subtypes -> MGE types
# ======================================================

mge_agg <- aggregate_merged_subtype(
  feature_df = mge_df,
  type_name = "MGE_type"
)


mge_type_df <- mge_agg$abundance


mge_mapping_df <- mge_agg$mapping %>%

  rename(
    MGE_type = type
  )


cat(
  "\n========================================\n"
)

cat(
  "MGE merged-subtype -> type aggregation\n"
)

cat(
  "========================================\n"
)


cat(
  "Merged subtypes:",
  nrow(mge_df),
  "\n"
)


cat(
  "MGE types/classes:",
  nrow(mge_type_df),
  "\n"
)


cat(
  "\nMGE types detected:\n"
)


print(
  mge_type_df$MGE_type
)


# ======================================================
# 10. Export subtype-to-type mappings
# ======================================================

write.csv(

  mrg_mapping_df,

  file.path(
    output_dir,
    "MRG_merged_subtype_to_type_mapping.csv"
  ),

  row.names = FALSE
)


write.csv(

  mge_mapping_df,

  file.path(
    output_dir,
    "MGE_merged_subtype_to_type_mapping.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 11. Export derived type abundance matrices
# ======================================================

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


write.table(

  mge_type_df,

  file.path(
    output_dir,
    "derived_MGE_type_abundance.tsv"
  ),

  sep = "\t",

  quote = FALSE,

  row.names = FALSE
)


# ======================================================
# 12. Function to construct sample × type matrix
# ======================================================

build_sample_matrix <- function(
    type_df
) {


  feature_names <- as.character(
    type_df[
      [
        1
      ]
    ]
  )


  x <- type_df[
    ,
    -1,
    drop = FALSE
  ]


  mat <- t(
    as.matrix(
      x
    )
  )


  rownames(mat) <- colnames(
    type_df
  )[-1]


  colnames(mat) <- make.unique(
    feature_names
  )


  mat[
    is.na(mat)
  ] <- 0


  mat[
    mat < 0
  ] <- 0


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
# 13. Build MRG and MGE type matrices
# ======================================================

mrg_mat <- build_sample_matrix(
  mrg_type_df
)


mge_mat <- build_sample_matrix(
  mge_type_df
)


cat(
  "\nMRG type matrix: ",
  nrow(mrg_mat),
  " samples × ",
  ncol(mrg_mat),
  " MRG types\n",
  sep = ""
)


cat(
  "MGE type matrix: ",
  nrow(mge_mat),
  " samples × ",
  ncol(mge_mat),
  " MGE types\n",
  sep = ""
)


# ======================================================
# 14. Check Copper / Zinc / Iron categories
# ======================================================

required_metals <- c(
  "Copper",
  "Zinc",
  "Iron"
)


missing_metals <- setdiff(
  required_metals,
  colnames(
    mrg_mat
  )
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
  "\nMRG categories used:\n"
)


print(
  required_metals
)


# ======================================================
# 15. Match MRG, MGE, and metadata samples
# ======================================================

common_samples <- Reduce(

  intersect,

  list(

    rownames(
      mrg_mat
    ),

    rownames(
      mge_mat
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
      "Please check sample IDs in MRG, MGE, and metadata."
    )
  )
}


mrg_mat <- mrg_mat[
  common_samples,
  ,
  drop = FALSE
]


mge_mat <- mge_mat[
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
    rownames(mrg_mat) ==
      rownames(mge_mat)
  )
)


stopifnot(
  all(
    rownames(mrg_mat) ==
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
# 16. Build Cu / Zn / Fe MRG abundance table
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


# ======================================================
# 17. Export metal abundance used in correlations
# ======================================================

write.csv(

  metal_df,

  file.path(
    output_dir,
    "D30D60D90_CuZnFe_MRG_abundance_used.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 18. Convert MGE type matrix to long format
# ======================================================

mge_long <- as.data.frame(
  mge_mat
) %>%

  tibble::rownames_to_column(
    "ID"
  ) %>%

  pivot_longer(

    cols =
      -ID,

    names_to =
      "MGE_type",

    values_to =
      "MGE_abundance"
  )


# ======================================================
# 19. Build correlation input table
# ======================================================

corr_input <- mge_long %>%

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
# 20. Safe Spearman correlation function
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
# 21. Group-specific Spearman correlations
#
# Group × Metal × MGE type
#
# D30, D60, and D90 samples are pooled within each
# treatment group.
# ======================================================

cor_df <- corr_input %>%

  group_by(
    Group,
    Metal,
    MGE_type
  ) %>%

  group_modify(

    ~ {


      cor_res <- safe_spearman(
        .x$MRG_abundance,
        .x$MGE_abundance
      )


      tibble(

        n =
          sum(
            complete.cases(
              .x$MRG_abundance,
              .x$MGE_abundance
            )
          ),


        rho =
          cor_res$rho,


        p_value =
          cor_res$p_value,


        mean_MGE =
          mean(
            .x$MGE_abundance,
            na.rm = TRUE
          ),


        detection_rate_MGE =
          mean(
            .x$MGE_abundance > 0,
            na.rm = TRUE
          ),


        mean_MRG =
          mean(
            .x$MRG_abundance,
            na.rm = TRUE
          ),


        detection_rate_MRG =
          mean(
            .x$MRG_abundance > 0,
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
# 22. BH/FDR correction
#
# Applied separately within each Group × Metal.
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
# 23. Export complete correlation results
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
    "D30D60D90_stat_CuZnFe_vs_MGEtype_Spearman_byGroup_all.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 24. Export significant correlations
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
    "D30D60D90_stat_CuZnFe_vs_MGEtype_Spearman_byGroup_FDR_significant.csv"
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
    "D30D60D90_stat_CuZnFe_vs_MGEtype_Spearman_byGroup_rawP_significant.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 25. Print significance summaries
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
      name =
        "n_FDR_lt_0.05"
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
      name =
        "n_rawP_lt_0.05"
    )
)


cat(
  "\nTop correlations by Group × Metal:\n"
)


print(

  cor_df %>%

    arrange(
      Group,
      Metal,
      p_adj
    ) %>%

    group_by(
      Group,
      Metal
    ) %>%

    slice_head(
      n = 5
    ) %>%

    ungroup()
)


# ======================================================
# 26. Select Top 10 MGE types
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


# ======================================================
# 27. Union of selected MGE types
# ======================================================

selected_mge_types <- unique(
  plot_select_df$MGE_type
)


plot_df <- cor_df %>%

  filter(
    MGE_type %in%
      selected_mge_types
  )


# ======================================================
# 28. Order MGE types
# ======================================================

mge_order <- plot_df %>%

  group_by(
    MGE_type
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
    MGE_type
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

    MGE_type = factor(
      MGE_type,
      levels = rev(
        mge_order
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
    "D30D60D90_stat_CuZnFe_vs_MGEtype_plot_data_top10.csv"
  ),

  row.names = FALSE
)


write.csv(

  plot_select_df,

  file.path(
    output_dir,
    "D30D60D90_stat_CuZnFe_vs_MGEtype_top10_selection.csv"
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
    y = MGE_type
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
      "MGE type",

    title =
      paste0(
        "D30-D90 group-specific associations between ",
        "Cu/Zn/Fe-related MRGs and MGE types"
      ),

    subtitle =
      paste0(
        "Top 10 MGE types per group and metal; ",
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
    selected_mge_types
  ) *
    0.22
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_D30D60D90_CuZnFe_vs_MGEtype_Spearman_dotplot.pdf"
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
    "Figure_D30D60D90_CuZnFe_vs_MGEtype_Spearman_dotplot.png"
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
  "D30-D90 MRG-MGE correlation analysis completed.\n"
)

cat(
  "MRG core input:\n"
)

cat(
  mrg_path,
  "\n"
)

cat(
  "MGE core input:\n"
)

cat(
  mge_path,
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