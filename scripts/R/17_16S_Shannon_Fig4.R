########################################################
## 16S alpha diversity: Shannon index
##
## Figure 4
##
## Input:
##   metadata/metadata.csv
##
##   processed_data/16S/
##     Shannon_index.tsv
##
## Shannon_index.tsv format:
##   ID    Shannon
##
## Analysis:
##
##   A. Within-group temporal comparisons:
##      D0 vs D30
##      D0 vs D60
##      D0 vs D90
##
##   B. Between-group comparisons at each time point:
##      ITM vs OTM1
##      ITM vs OTM2
##      OTM1 vs OTM2
##
## Statistical analysis:
##   Wilcoxon rank-sum test (unpaired)
##
## Multiple-testing correction:
##   Benjamini-Hochberg (BH/FDR)
##
## Outlier handling:
##   No samples are removed.
##
## Output:
##   results/Fig4_16S_Shannon/
########################################################


# ======================================================
# 1. Load packages
# ======================================================

library(tidyverse)
library(ggpubr)
library(patchwork)
library(colorspace)
library(cowplot)
library(grid)


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


shannon_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "16S",
  "Shannon_index.tsv"
)


output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig4_16S_Shannon"
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


required_meta_cols <- c(
  "ID",
  "Time",
  "Group"
)


missing_meta_cols <- setdiff(
  required_meta_cols,
  colnames(meta_df)
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


# ======================================================
# 4. Prepare metadata
#
# Expected metadata:
#
# ID | Time | Group | Sample
#
# Group:
# ITM, OTM1, OTM2
#
# Time:
# D0, D30, D60, D90
# ======================================================

meta_df <- meta_df %>%

  mutate(

    ID = as.character(ID),

    Time = as.character(Time),

    Group = as.character(Group)
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


expected_times <- c(
  "D0",
  "D30",
  "D60",
  "D90"
)


unexpected_times <- setdiff(
  unique(meta_df$Time),
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


if (anyDuplicated(meta_df$ID) > 0) {

  stop(
    "Duplicated sequencing sample IDs were found in metadata.csv."
  )
}


# ======================================================
# 5. Read Shannon index
#
# Expected format:
#
# ID    Shannon
# a1    4.12
# a2    4.36
# ...
# ======================================================

shannon_df <- read.delim(
  shannon_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


required_shannon_cols <- c(
  "ID",
  "Shannon"
)


missing_shannon_cols <- setdiff(
  required_shannon_cols,
  colnames(shannon_df)
)


if (length(missing_shannon_cols) > 0) {

  stop(
    paste0(
      "Missing Shannon_index.tsv column(s): ",
      paste(
        missing_shannon_cols,
        collapse = ", "
      ),
      ". Expected columns: ID and Shannon."
    )
  )
}


shannon_df <- shannon_df %>%

  select(
    ID,
    Shannon
  ) %>%

  mutate(

    ID = as.character(ID),

    Shannon = as.numeric(
      Shannon
    )
  )


if (anyDuplicated(shannon_df$ID) > 0) {

  stop(
    "Duplicated sequencing sample IDs were found in Shannon_index.tsv."
  )
}


if (any(is.na(shannon_df$Shannon))) {

  warning(
    "Missing or non-numeric Shannon values were detected."
  )
}


# ======================================================
# 6. Match Shannon data with metadata
# ======================================================

common_samples <- intersect(
  meta_df$ID,
  shannon_df$ID
)


cat(
  "\n========================================\n"
)

cat(
  "16S Shannon alpha-diversity analysis\n"
)

cat(
  "========================================\n"
)


cat(
  "Metadata samples:",
  nrow(meta_df),
  "\n"
)


cat(
  "Shannon samples:",
  nrow(shannon_df),
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


missing_in_shannon <- setdiff(
  meta_df$ID,
  shannon_df$ID
)


if (length(missing_in_shannon) > 0) {

  cat(
    "\nSamples present in metadata but missing from Shannon_index.tsv:\n"
  )

  print(
    missing_in_shannon
  )
}


missing_in_metadata <- setdiff(
  shannon_df$ID,
  meta_df$ID
)


if (length(missing_in_metadata) > 0) {

  cat(
    "\nSamples present in Shannon_index.tsv but missing from metadata:\n"
  )

  print(
    missing_in_metadata
  )
}


if (length(common_samples) < 5) {

  stop(
    paste0(
      "Too few matched samples. ",
      "Please check IDs in metadata.csv and Shannon_index.tsv."
    )
  )
}


# ======================================================
# 7. Merge Shannon index with metadata
# ======================================================

df <- shannon_df %>%

  inner_join(

    meta_df %>%

      select(
        ID,
        Time,
        Group
      ),

    by = "ID"
  ) %>%

  filter(
    !is.na(Shannon),
    !is.na(Time),
    !is.na(Group)
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

    Time = factor(
      Time,
      levels = c(
        "D0",
        "D30",
        "D60",
        "D90"
      )
    ),

    GroupTime = interaction(
      Group,
      Time,
      sep = "."
    )
  )


# ======================================================
# 8. Check final analysis dataset
# ======================================================

cat(
  "\nSamples by Group × Time:\n"
)


print(
  table(
    df$Group,
    df$Time
  )
)


cat(
  "\nFinal number of samples:",
  nrow(df),
  "\n"
)


# Export exact dataset used in analysis
write.csv(

  df %>%

    select(
      ID,
      Time,
      Group,
      Shannon
    ),

  file.path(
    output_dir,
    "Shannon_analysis_input.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 9. Colors
# ======================================================

group_cols <- c(

  ITM =
    "#8ECFC9",

  OTM1 =
    "#FFBE7A",

  OTM2 =
    "#FA7F6F"
)


make_group_gradient <- function(
    base
) {

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


grad_cols <- c(

  setNames(

    make_group_gradient(
      group_cols[
        "ITM"
      ]
    ),

    paste0(
      "ITM.",
      levels(
        df$Time
      )
    )
  ),


  setNames(

    make_group_gradient(
      group_cols[
        "OTM1"
      ]
    ),

    paste0(
      "OTM1.",
      levels(
        df$Time
      )
    )
  ),


  setNames(

    make_group_gradient(
      group_cols[
        "OTM2"
      ]
    ),

    paste0(
      "OTM2.",
      levels(
        df$Time
      )
    )
  )
)


# ======================================================
# 10. Pairwise comparisons
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
# 11. Figure theme
# ======================================================

theme_main <- theme_classic(
  base_size = 9.5
) +

  theme(

    legend.position =
      "none",


    axis.title =
      element_text(
        face = "bold",
        size = 9.5
      ),


    axis.text =
      element_text(
        color = "black",
        size = 8.0
      ),


    axis.line =
      element_line(
        linewidth = 0.38
      ),


    axis.ticks =
      element_line(
        linewidth = 0.38
      ),


    axis.ticks.length =
      unit(
        0.07,
        "cm"
      ),


    strip.text =
      element_text(
        face = "bold",
        size = 8.8,
        margin = margin(
          2,
          2,
          2,
          2
        )
      ),


    strip.background =
      element_rect(
        fill = "white",
        color = "black",
        linewidth = 0.35
      ),


    panel.spacing.x =
      unit(
        0.18,
        "lines"
      ),


    panel.spacing.y =
      unit(
        0.15,
        "lines"
      ),


    plot.margin =
      margin(
        2,
        2,
        2,
        2
      )
  )


theme_group_x <- theme(

  axis.text.x =
    element_text(
      angle = 20,
      hjust = 1,
      vjust = 1,
      size = 7.8
    )
)


# ======================================================
# 12. Boxplots + individual samples
#
# No outlier removal is performed.
#
# outlier.shape = NA only prevents duplicated plotting
# of boxplot outlier symbols because every sample is
# already displayed using jittered points.
# ======================================================

geom_box_pts <- list(

  geom_boxplot(

    width = 0.78,

    outlier.shape = NA,

    color = "black",

    linewidth = 0.38
  ),


  geom_point(

    position =
      position_jitter(
        width = 0.07,
        height = 0
      ),

    size = 1.15,

    shape = 21,

    stroke = 0.28,

    color = "black",

    alpha = 0.82
  )
)


# ======================================================
# 13. Y-axis limits
# ======================================================

get_y_limits <- function(
    data,
    n_sig_max = 0
) {


  y_min <- min(
    data$Shannon,
    na.rm = TRUE
  )


  y_max <- max(
    data$Shannon,
    na.rm = TRUE
  )


  y_range <- y_max - y_min


  if (
    is.na(y_range) ||
    y_range == 0
  ) {

    y_range <- 1
  }


  lower <- ifelse(

    y_min >= 0,

    max(
      0,
      y_min -
        0.025 *
        y_range
    ),

    y_min -
      0.04 *
      y_range
  )


  upper <-
    y_max +
    y_range *
    (
      0.08 +
      0.055 *
      max(
        n_sig_max,
        1
      )
    )


  c(
    lower,
    upper
  )
}


# ======================================================
# 14. Within-group temporal comparisons
#
# D0 vs D30
# D0 vs D60
# D0 vs D90
#
# Wilcoxon rank-sum test (unpaired)
# BH-adjusted P values
# ======================================================

get_sig_within_time <- function(
    data
) {


  stat.test <- ggpubr::compare_means(

    Shannon ~ Time,

    data = data,

    group.by =
      "Group",

    method =
      "wilcox.test",

    paired =
      FALSE,

    comparisons =
      time_comparisons,

    p.adjust.method =
      "BH"
  ) %>%

    filter(
      p.adj < 0.05
    ) %>%

    mutate(

      p.signif = case_when(

        p.adj < 0.001 ~
          "***",

        p.adj < 0.01 ~
          "**",

        TRUE ~
          "*"
      )
    )


  if (
    nrow(
      stat.test
    ) == 0
  ) {

    return(
      stat.test
    )
  }


  y_stats <- data %>%

    group_by(
      Group
    ) %>%

    summarise(

      y_min =
        min(
          Shannon,
          na.rm = TRUE
        ),

      y_max =
        max(
          Shannon,
          na.rm = TRUE
        ),

      y_range =
        y_max -
        y_min,

      .groups =
        "drop"
    ) %>%

    mutate(

      y_range =
        ifelse(
          is.na(
            y_range
          ) |
            y_range == 0,
          1,
          y_range
        )
    )


  stat.test <- stat.test %>%

    left_join(
      y_stats,
      by = "Group"
    ) %>%

    group_by(
      Group
    ) %>%

    arrange(
      p.adj,
      .by_group = TRUE
    ) %>%

    mutate(

      y.position =
        y_max +
        y_range *
        (
          0.045 +
          0.050 *
          (
            row_number() -
            1
          )
        )
    ) %>%

    ungroup()


  return(
    stat.test
  )
}


add_sig_within_time <- function(
    p,
    data
) {


  stat.test <- get_sig_within_time(
    data
  )


  if (
    nrow(
      stat.test
    ) == 0
  ) {

    return(
      p
    )
  }


  p +

    ggpubr::stat_pvalue_manual(

      stat.test,

      label =
        "p.signif",

      tip.length =
        0.004,

      size =
        3.0,

      bracket.size =
        0.28
    )
}


# ======================================================
# 15. Between-group comparisons at each time point
#
# ITM vs OTM1
# ITM vs OTM2
# OTM1 vs OTM2
#
# Wilcoxon rank-sum test (unpaired)
# BH-adjusted P values
# ======================================================

get_sig_between_group <- function(
    data
) {


  stat.test <- ggpubr::compare_means(

    Shannon ~ Group,

    data = data,

    group.by =
      "Time",

    method =
      "wilcox.test",

    paired =
      FALSE,

    comparisons =
      group_comparisons,

    p.adjust.method =
      "BH"
  ) %>%

    filter(
      p.adj < 0.05
    ) %>%

    mutate(

      p.signif = case_when(

        p.adj < 0.001 ~
          "***",

        p.adj < 0.01 ~
          "**",

        TRUE ~
          "*"
      )
    )


  if (
    nrow(
      stat.test
    ) == 0
  ) {

    return(
      stat.test
    )
  }


  y_stats <- data %>%

    group_by(
      Time
    ) %>%

    summarise(

      y_min =
        min(
          Shannon,
          na.rm = TRUE
        ),

      y_max =
        max(
          Shannon,
          na.rm = TRUE
        ),

      y_range =
        y_max -
        y_min,

      .groups =
        "drop"
    ) %>%

    mutate(

      y_range =
        ifelse(
          is.na(
            y_range
          ) |
            y_range == 0,
          1,
          y_range
        )
    )


  stat.test <- stat.test %>%

    left_join(
      y_stats,
      by = "Time"
    ) %>%

    group_by(
      Time
    ) %>%

    arrange(
      p.adj,
      .by_group = TRUE
    ) %>%

    mutate(

      y.position =
        y_max +
        y_range *
        (
          0.045 +
          0.050 *
          (
            row_number() -
            1
          )
        )
    ) %>%

    ungroup()


  return(
    stat.test
  )
}


add_sig_between_group <- function(
    p,
    data
) {


  stat.test <- get_sig_between_group(
    data
  )


  if (
    nrow(
      stat.test
    ) == 0
  ) {

    return(
      p
    )
  }


  p +

    ggpubr::stat_pvalue_manual(

      stat.test,

      label =
        "p.signif",

      tip.length =
        0.004,

      size =
        3.0,

      bracket.size =
        0.28
    )
}


# ======================================================
# 16. Calculate statistics
# ======================================================

sig_A <- get_sig_within_time(
  df
)


sig_B <- get_sig_between_group(
  df
)


n_sig_A <- ifelse(

  nrow(
    sig_A
  ) == 0,

  0,

  max(
    table(
      sig_A$Group
    )
  )
)


n_sig_B <- ifelse(

  nrow(
    sig_B
  ) == 0,

  0,

  max(
    table(
      sig_B$Time
    )
  )
)


ylim_A <- get_y_limits(
  df,
  n_sig_A
)


ylim_B <- get_y_limits(
  df,
  n_sig_B
)


# ======================================================
# 17. Export statistical results
# ======================================================

write.csv(

  sig_A,

  file.path(
    output_dir,
    "stat_Shannon_within_group_time_Wilcoxon_BH_significant.csv"
  ),

  row.names = FALSE
)


write.csv(

  sig_B,

  file.path(
    output_dir,
    "stat_Shannon_between_group_Wilcoxon_BH_significant.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 18. Panel A
#
# Within-group temporal comparisons
# ======================================================

pA <- ggplot(

  df,

  aes(
    x = Time,
    y = Shannon,
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

    x = NULL,

    y = "Shannon index"
  ) +

  coord_cartesian(

    ylim = ylim_A,

    clip = "off"
  ) +

  theme_main


pA <- add_sig_within_time(
  pA,
  df
)


# ======================================================
# 19. Panel B
#
# Between-group comparisons at each time point
# ======================================================

pB <- ggplot(

  df,

  aes(
    x = Group,
    y = Shannon,
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

    x = NULL,

    y = "Shannon index"
  ) +

  coord_cartesian(

    ylim = ylim_B,

    clip = "off"
  ) +

  theme_main +

  theme_group_x


pB <- add_sig_between_group(
  pB,
  df
)


# ======================================================
# 20. Combine panels
# ======================================================

core_plot <- patchwork::wrap_plots(

  pA,

  pB,

  ncol =
    2,

  widths =
    c(
      1,
      1.12
    )
) +

  patchwork::plot_annotation(
    tag_levels = "A"
  ) &

  theme(

    plot.tag =
      element_text(
        face = "bold",
        size = 11
      ),

    plot.margin =
      margin(
        1,
        1,
        1,
        1
      )
  )


# ======================================================
# 21. Group legend
# ======================================================

p_leg <- ggplot(

  df,

  aes(
    x = Group,
    y = Shannon,
    fill = Group
  )

) +

  geom_boxplot(
    linewidth = 0.35
  ) +

  scale_fill_manual(

    values =
      group_cols,

    name =
      "Group"
  ) +

  theme_void(
    base_size = 8.8
  ) +

  theme(

    legend.position =
      "top",

    legend.title =
      element_text(
        face = "bold",
        size = 8.2
      ),

    legend.text =
      element_text(
        color = "black",
        size = 7.8
      ),

    legend.key.size =
      unit(
        0.28,
        "cm"
      ),

    legend.spacing.x =
      unit(
        0.08,
        "cm"
      ),

    legend.margin =
      margin(
        0,
        0,
        0,
        0
      )
  )


legend_grob <- cowplot::get_legend(
  p_leg
)


# ======================================================
# 22. Final figure
# ======================================================

final_plot <- cowplot::ggdraw(
  core_plot
) +

  cowplot::draw_grob(

    legend_grob,

    x = 0.66,

    y = 0.925,

    width = 0.27,

    height = 0.06
  )


print(
  final_plot
)


# ======================================================
# 23. Save figure
# ======================================================

ggsave(

  filename = file.path(
    output_dir,
    "Figure_16S_Shannon_BH.pdf"
  ),

  plot = final_plot,

  width = 10.2,

  height = 4.0
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_16S_Shannon_BH.png"
  ),

  plot = final_plot,

  width = 10.2,

  height = 4.0,

  dpi = 600
)


# ======================================================
# 24. Finish
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "16S Shannon alpha-diversity analysis completed.\n"
)

cat(
  "Panel A: within-group temporal comparisons.\n"
)

cat(
  "Panel B: between-group comparisons at each time point.\n"
)

cat(
  "Wilcoxon rank-sum tests (unpaired); BH adjustment.\n"
)

cat(
  "No outlier removal; all available samples retained.\n"
)

cat(
  "Shannon input:\n"
)

cat(
  shannon_path,
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