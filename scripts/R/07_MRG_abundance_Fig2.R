########################################################
## Metal resistance gene (MRG) abundance analysis
##
## Figure 2
##
## Analyses:
##   A-B: Total MRG abundance
##   C-D: Copper-related MRG abundance
##   E-F: Zinc-related MRG abundance
##
## Panel structure:
##   Left:  within-group comparisons over time
##   Right: between-group comparisons at each time point
##
## Statistical analysis:
##   Wilcoxon rank-sum test (unpaired)
##   Benjamini-Hochberg (BH/FDR) correction
##
## Outlier filtering:
##   Upper-tail values > Q3 + 3 × IQR
##   within each Group × Time combination
##
## Input:
##   metadata/metadata.csv
##   processed_data/MRG/
##     MRG_merged_subtype_abundance.tsv
########################################################


# ======================================================
# 1. Load packages
# ======================================================

library(tidyverse)
library(ggpubr)
library(rstatix)
library(colorspace)
library(cowplot)


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


feature_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "MRG",
  "MRG_merged_subtype_abundance.tsv"
)


output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig2_MRG_abundance"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


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


meta_df <- meta_df %>%

  mutate(
    ID = as.character(ID),
    Group = as.character(Group),
    Time = as.character(Time)
  )


# ======================================================
# 4. Check metadata
# ======================================================

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


cat(
  "\n========================================\n"
)

cat(
  "Metadata check\n"
)

cat(
  "========================================\n"
)


cat(
  "Number of metadata rows:",
  nrow(meta_df),
  "\n"
)


cat(
  "Number of unique samples:",
  length(
    unique(meta_df$ID)
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


if (anyDuplicated(meta_df$ID) > 0) {

  stop(
    "Duplicated sequencing sample IDs were found."
  )
}


# ======================================================
# 5. Read cell number-normalized MRG subtype matrix
#
# First column:
# MRG subtype / metal category
#
# Remaining columns:
# sample IDs
# ======================================================

feature_df <- read.delim(
  feature_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ======================================================
# 6. Construct sample × MRG subtype matrix
# ======================================================

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


colnames(mat) <- as.character(
  feature_df[[1]]
)


# ======================================================
# 7. Match samples between metadata and MRG matrix
# ======================================================

common_samples <- intersect(
  meta_df$ID,
  rownames(mat)
)


cat(
  "\nSamples matched between metadata and MRG matrix:",
  length(common_samples),
  "\n"
)


if (length(common_samples) == 0) {

  stop(
    "No matched samples between metadata and MRG matrix."
  )
}


if (length(common_samples) != 120) {

  warning(
    paste0(
      "Expected 120 samples, but ",
      length(common_samples),
      " were matched."
    )
  )
}


missing_samples <- setdiff(
  meta_df$ID,
  rownames(mat)
)


if (length(missing_samples) > 0) {

  cat(
    "\nSamples missing from MRG matrix:\n"
  )

  print(
    missing_samples
  )
}


mat <- mat[
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
    rownames(mat) ==
      meta_df$ID
  )
)


# ======================================================
# 8. Colors
# ======================================================

group_cols <- c(

  ITM =
    "#8ECFC9",

  OTM1 =
    "#FFBE7A",

  OTM2 =
    "#FA7F6F"
)


make_group_gradient <- function(base) {

  c(

    colorspace::lighten(
      base,
      0.65
    ),

    colorspace::lighten(
      base,
      0.40
    ),

    colorspace::lighten(
      base,
      0.20
    ),

    base
  )
}


# ======================================================
# 9. Define comparisons
# ======================================================

time_comparisons <- list(

  c(
    "D0",
    "D30"
  ),

  c(
    "D0",
    "D60"
  ),

  c(
    "D0",
    "D90"
  )
)


group_comparisons <- list(

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
# 10. Plot theme
# ======================================================

theme_main <- theme_classic(
  base_size = 14
) +

  theme(

    legend.position =
      "none",

    axis.title =
      element_text(
        face = "bold"
      ),

    axis.text =
      element_text(
        color = "black"
      ),

    axis.line =
      element_line(
        linewidth = 0.9
      ),

    axis.ticks =
      element_line(
        linewidth = 0.9
      ),

    strip.text =
      element_text(
        face = "bold"
      )
  )


theme_between <- theme(

  axis.text.x =
    element_text(
      angle = 20,
      hjust = 1,
      vjust = 1
    )
)


# ======================================================
# 11. Boxplot and sample-point style
# ======================================================

geom_box_pts <- list(

  geom_boxplot(

    width = 0.65,

    outlier.shape = NA,

    color = "black",

    linewidth = 1
  ),


  geom_point(

    position =
      position_jitter(
        width = 0.12
      ),

    size = 2.2,

    shape = 21,

    stroke = 0.6,

    color = "black",

    alpha = 0.85
  )
)


# ======================================================
# 12. Extract abundance for a specific metal category
#
# Exact matching is attempted first.
# If no exact match is found, case-insensitive pattern
# matching is used and all matching columns are summed.
# ======================================================

extract_metal_sum <- function(
    mat,
    metal_name
) {


  cn <- colnames(
    mat
  )


  # Exact match
  idx_exact <- which(
    tolower(cn) ==
      tolower(metal_name)
  )


  if (length(idx_exact) > 0) {

    v <- rowSums(
      mat[
        ,
        idx_exact,
        drop = FALSE
      ],
      na.rm = TRUE
    )


    names(v) <- rownames(
      mat
    )


    return(
      v
    )
  }


  # Pattern match
  idx_pattern <- grep(
    metal_name,
    cn,
    ignore.case = TRUE
  )


  if (length(idx_pattern) > 0) {

    message(
      "No exact column named '",
      metal_name,
      "' was found; summing ",
      length(idx_pattern),
      " matched subtype(s)."
    )


    v <- rowSums(
      mat[
        ,
        idx_pattern,
        drop = FALSE
      ],
      na.rm = TRUE
    )


    names(v) <- rownames(
      mat
    )


    return(
      v
    )
  }


  stop(
    "No subtype corresponding to '",
    metal_name,
    "' was found in the MRG matrix."
  )
}


# ======================================================
# 13. Build abundance data frame
#
# Original analysis removes only upper-tail outliers:
# Abundance > Q3 + 3 × IQR
# within each Group × Time.
# ======================================================

build_df_from_vector <- function(
    abund_vec,
    meta_df,
    label
) {


  df <- tibble(

    ID = names(
      abund_vec
    ),

    Abundance = as.numeric(
      abund_vec
    )

  ) %>%

    left_join(
      meta_df,
      by = "ID"
    )


  df <- df %>%

    group_by(
      Group,
      Time
    ) %>%

    mutate(

      Q1 = quantile(
        Abundance,
        0.25,
        na.rm = TRUE
      ),

      Q3 = quantile(
        Abundance,
        0.75,
        na.rm = TRUE
      ),

      IQR_value =
        Q3 -
        Q1,

      upper_3IQR =
        Q3 +
        3 *
        IQR_value
    ) %>%

    ungroup() %>%

    filter(
      Abundance <=
        upper_3IQR
    ) %>%

    select(
      -Q1,
      -Q3,
      -IQR_value,
      -upper_3IQR
    )


  df$Metric <- label


  return(
    df
  )
}


# ======================================================
# 14. Add within-group temporal significance
#
# Unpaired Wilcoxon tests are retained to reproduce
# the original analysis.
# ======================================================

add_sig_within_time <- function(
    p,
    data,
    yvar
) {


  stat.test <- ggpubr::compare_means(

    formula =
      as.formula(
        paste0(
          yvar,
          " ~ Time"
        )
      ),

    data = data,

    group.by =
      "Group",

    method =
      "wilcox.test",

    comparisons =
      time_comparisons

  ) %>%

    rstatix::adjust_pvalue(
      method = "BH"
    ) %>%

    filter(
      p.adj < 0.05
    ) %>%

    mutate(

      p.signif =
        case_when(

          p.adj < 0.001 ~
            "***",

          p.adj < 0.01 ~
            "**",

          TRUE ~
            "*"
        )
    )


  if (nrow(stat.test) == 0) {

    return(
      p
    )
  }


  ymax <- data %>%

    group_by(
      Group
    ) %>%

    summarise(

      max_val = max(
        .data[[yvar]],
        na.rm = TRUE
      ),

      .groups = "drop"
    )


  stat.test <- stat.test %>%

    left_join(
      ymax,
      by = "Group"
    ) %>%

    group_by(
      Group
    ) %>%

    mutate(

      y.position =
        max_val +
        (
          0.08 *
            max_val +
            1e-8
        ) *
        row_number()
    ) %>%

    ungroup()


  p +

    ggpubr::stat_pvalue_manual(

      stat.test,

      label =
        "p.signif",

      tip.length =
        0.01,

      size =
        5
    )
}


# ======================================================
# 15. Add between-group significance
# ======================================================

add_sig_between_group <- function(
    p,
    data,
    yvar
) {


  stat.test <- ggpubr::compare_means(

    formula =
      as.formula(
        paste0(
          yvar,
          " ~ Group"
        )
      ),

    data = data,

    group.by =
      "Time",

    method =
      "wilcox.test",

    comparisons =
      group_comparisons

  ) %>%

    rstatix::adjust_pvalue(
      method = "BH"
    ) %>%

    filter(
      p.adj < 0.05
    ) %>%

    mutate(

      p.signif =
        case_when(

          p.adj < 0.001 ~
            "***",

          p.adj < 0.01 ~
            "**",

          TRUE ~
            "*"
        )
    )


  if (nrow(stat.test) == 0) {

    return(
      p
    )
  }


  ymax <- data %>%

    group_by(
      Time
    ) %>%

    summarise(

      max_val = max(
        .data[[yvar]],
        na.rm = TRUE
      ),

      .groups = "drop"
    )


  stat.test <- stat.test %>%

    left_join(
      ymax,
      by = "Time"
    ) %>%

    group_by(
      Time
    ) %>%

    mutate(

      y.position =
        max_val +
        (
          0.08 *
            max_val +
            1e-8
        ) *
        row_number()
    ) %>%

    ungroup()


  p +

    ggpubr::stat_pvalue_manual(

      stat.test,

      label =
        "p.signif",

      tip.length =
        0.01,

      size =
        5
    )
}


# ======================================================
# 16. Generate one pair of abundance plots
#
# Left:
# within-treatment comparisons over time
#
# Right:
# between-treatment comparisons at each time point
# ======================================================

plot_pair <- function(
    df,
    ylab
) {


  grad_cols <- c(

    setNames(

      make_group_gradient(
        group_cols["ITM"]
      ),

      paste0(
        "ITM.",
        levels(df$Time)
      )
    ),


    setNames(

      make_group_gradient(
        group_cols["OTM1"]
      ),

      paste0(
        "OTM1.",
        levels(df$Time)
      )
    ),


    setNames(

      make_group_gradient(
        group_cols["OTM2"]
      ),

      paste0(
        "OTM2.",
        levels(df$Time)
      )
    )
  )


  df$GroupTime <- interaction(
    df$Group,
    df$Time,
    sep = "."
  )


  ymax <- max(
    df$Abundance,
    na.rm = TRUE
  )


  upper_lim <-
    ymax *
    1.35


  # ----------------------------------
  # Left panel
  # ----------------------------------

  p_left <- ggplot(

    df,

    aes(
      Time,
      Abundance,
      fill = GroupTime
    )

  ) +

    geom_box_pts +

    facet_wrap(
      ~ Group,
      nrow = 1
    ) +

    scale_fill_manual(
      values = grad_cols
    ) +

    labs(
      y = ylab,
      x = NULL
    ) +

    coord_cartesian(
      ylim = c(
        0,
        upper_lim
      )
    ) +

    theme_main


  p_left <- add_sig_within_time(
    p_left,
    df,
    "Abundance"
  )


  # ----------------------------------
  # Right panel
  # ----------------------------------

  p_right <- ggplot(

    df,

    aes(
      Group,
      Abundance,
      fill = Group
    )

  ) +

    geom_box_pts +

    facet_wrap(
      ~ Time,
      nrow = 1
    ) +

    scale_fill_manual(
      values = group_cols
    ) +

    labs(
      y = ylab,
      x = NULL
    ) +

    coord_cartesian(
      ylim = c(
        0,
        upper_lim
      )
    ) +

    theme_main +

    theme_between


  p_right <- add_sig_between_group(
    p_right,
    df,
    "Abundance"
  )


  return(
    list(
      left = p_left,
      right = p_right
    )
  )
}


# ======================================================
# 17. Calculate total MRG abundance
# ======================================================

abund_total <- rowSums(
  mat,
  na.rm = TRUE
)


names(
  abund_total
) <- rownames(
  mat
)


df_total <- build_df_from_vector(
  abund_total,
  meta_df,
  "Total"
)


# ======================================================
# 18. Calculate copper-related MRG abundance
# ======================================================

abund_cu <- extract_metal_sum(
  mat,
  "Copper"
)


df_cu <- build_df_from_vector(
  abund_cu,
  meta_df,
  "Copper"
)


# ======================================================
# 19. Calculate zinc-related MRG abundance
# ======================================================

abund_zn <- extract_metal_sum(
  mat,
  "Zinc"
)


df_zn <- build_df_from_vector(
  abund_zn,
  meta_df,
  "Zinc"
)


# ======================================================
# 20. Generate individual plot pairs
# ======================================================

plots_total <- plot_pair(

  df_total,

  paste0(
    "Total metal resistance abundance\n",
    "(cell number-normalized)"
  )
)


plots_cu <- plot_pair(

  df_cu,

  paste0(
    "Copper-related resistance abundance\n",
    "(cell number-normalized)"
  )
)


plots_zn <- plot_pair(

  df_zn,

  paste0(
    "Zinc-related resistance abundance\n",
    "(cell number-normalized)"
  )
)


# ======================================================
# 21. Shared legend
# ======================================================

df_leg <- bind_rows(
  df_total,
  df_cu,
  df_zn
)


p_leg <- ggplot(

  df_leg,

  aes(
    Group,
    Abundance,
    fill = Group
  )

) +

  geom_boxplot() +

  scale_fill_manual(

    values =
      group_cols,

    name =
      "Group"
  ) +

  theme_void() +

  theme(
    legend.position =
      "top"
  )


legend_grob <- cowplot::get_legend(
  p_leg
)


# ======================================================
# 22. Figure A-B: Total MRG abundance
# ======================================================

fig_AB <- cowplot::plot_grid(

  plots_total$left,

  plots_total$right,

  ncol = 2,

  labels = c(
    "A",
    "B"
  ),

  label_size =
    16,

  label_fontface =
    "bold"
)


fig_AB_final <- cowplot::ggdraw(
  fig_AB
) +

  cowplot::draw_grob(

    legend_grob,

    x = 0.55,

    y = 0.95,

    width = 0.40,

    height = 0.08
  )


print(
  fig_AB_final
)


# ======================================================
# 23. Figure C-D: Copper-related MRG abundance
# ======================================================

fig_CD <- cowplot::plot_grid(

  plots_cu$left,

  plots_cu$right,

  ncol = 2,

  labels = c(
    "C",
    "D"
  ),

  label_size =
    16,

  label_fontface =
    "bold"
)


fig_CD_final <- cowplot::ggdraw(
  fig_CD
) +

  cowplot::draw_grob(

    legend_grob,

    x = 0.55,

    y = 0.95,

    width = 0.40,

    height = 0.08
  )


print(
  fig_CD_final
)


# ======================================================
# 24. Figure E-F: Zinc-related MRG abundance
# ======================================================

fig_EF <- cowplot::plot_grid(

  plots_zn$left,

  plots_zn$right,

  ncol = 2,

  labels = c(
    "E",
    "F"
  ),

  label_size =
    16,

  label_fontface =
    "bold"
)


fig_EF_final <- cowplot::ggdraw(
  fig_EF
) +

  cowplot::draw_grob(

    legend_grob,

    x = 0.55,

    y = 0.95,

    width = 0.40,

    height = 0.08
  )


print(
  fig_EF_final
)


# ======================================================
# 25. Combined A-F figure
# ======================================================

combined_core <- cowplot::plot_grid(

  plots_total$left,

  plots_total$right,

  plots_cu$left,

  plots_cu$right,

  plots_zn$left,

  plots_zn$right,

  ncol = 2,

  labels = c(
    "A",
    "B",
    "C",
    "D",
    "E",
    "F"
  ),

  label_size =
    16,

  label_fontface =
    "bold"
)


combined_final <- cowplot::ggdraw(
  combined_core
) +

  cowplot::draw_grob(

    legend_grob,

    x = 0.55,

    y = 0.965,

    width = 0.40,

    height = 0.06
  )


print(
  combined_final
)


# ======================================================
# 26. Save figures
# ======================================================

ggsave(

  filename = file.path(
    output_dir,
    "Figure_AB_Total_MRG_normalized_BH.pdf"
  ),

  plot =
    fig_AB_final,

  width =
    14,

  height =
    5
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_CD_Copper_MRG_normalized_BH.pdf"
  ),

  plot =
    fig_CD_final,

  width =
    14,

  height =
    5
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_EF_Zinc_MRG_normalized_BH.pdf"
  ),

  plot =
    fig_EF_final,

  width =
    14,

  height =
    5
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_ABCDEF_MRG_combined_normalized_BH.pdf"
  ),

  plot =
    combined_final,

  width =
    14,

  height =
    15
)


# ======================================================
# 27. Statistical-result export function
# ======================================================

save_stats <- function(
    df,
    prefix
) {


  # Within-treatment temporal comparisons
  stat_within <- ggpubr::compare_means(

    Abundance ~ Time,

    data =
      df,

    group.by =
      "Group",

    method =
      "wilcox.test",

    comparisons =
      time_comparisons

  ) %>%

    rstatix::adjust_pvalue(
      method = "BH"
    )


  # Between-treatment comparisons
  stat_between <- ggpubr::compare_means(

    Abundance ~ Group,

    data =
      df,

    group.by =
      "Time",

    method =
      "wilcox.test",

    comparisons =
      group_comparisons

  ) %>%

    rstatix::adjust_pvalue(
      method = "BH"
    )


  write.csv(

    stat_within,

    file.path(
      output_dir,
      paste0(
        "stat_",
        prefix,
        "_withinGroup_overTime_BH.csv"
      )
    ),

    row.names =
      FALSE
  )


  write.csv(

    stat_between,

    file.path(
      output_dir,
      paste0(
        "stat_",
        prefix,
        "_betweenGroup_sameTime_BH.csv"
      )
    ),

    row.names =
      FALSE
  )
}


# ======================================================
# 28. Export statistical results
# ======================================================

save_stats(
  df_total,
  "total_MRG_normalized"
)


save_stats(
  df_cu,
  "copper_MRG_normalized"
)


save_stats(
  df_zn,
  "zinc_MRG_normalized"
)


# ======================================================
# 29. Export datasets used in plotting
# ======================================================

write.csv(

  df_total,

  file.path(
    output_dir,
    "data_total_MRG_abundance_used.csv"
  ),

  row.names =
    FALSE
)


write.csv(

  df_cu,

  file.path(
    output_dir,
    "data_copper_MRG_abundance_used.csv"
  ),

  row.names =
    FALSE
)


write.csv(

  df_zn,

  file.path(
    output_dir,
    "data_zinc_MRG_abundance_used.csv"
  ),

  row.names =
    FALSE
)


# ======================================================
# 30. Finish
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "MRG abundance analysis completed.\n"
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