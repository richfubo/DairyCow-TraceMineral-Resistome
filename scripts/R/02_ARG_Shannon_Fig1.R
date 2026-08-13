########################################################
## ARG Shannon diversity
##
## Figure 1
##
## Panel A:
## Within-group comparisons over time
## D0 vs D30, D0 vs D60, D0 vs D90
##
## Panel B:
## Between-group comparisons at each time point
##
## Statistical analysis:
## Wilcoxon rank-sum test (unpaired)
## Benjamini-Hochberg (BH/FDR) correction
##
## Outlier filtering:
## 3 × IQR within each Group × Time
##
## Metadata columns:
## ID | Time | Group | Sample
##
## Group:
## ITM, OTM1, OTM2
##
## Sample represents the cow identifier but is not used
## for paired statistical testing in this analysis.
########################################################


# ======================================================
# 1. Load packages
# ======================================================

library(tidyverse)
library(ggpubr)
library(rstatix)
library(patchwork)
library(colorspace)
library(cowplot)
library(vegan)
library(tibble)
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


feature_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "ARG",
  "normalized_cell.subtype_matrix.tsv"
)


output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig1_ARG_Shannon"
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
  "Time",
  "Group",
  "Sample"
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

    Time = as.character(Time),

    Group = as.character(Group),

    Sample = as.character(Sample)
  )


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


# ======================================================
# 4. Metadata checks
# ======================================================

cat(
  "Metadata rows:",
  nrow(meta_df),
  "\n"
)


cat(
  "Unique sequencing samples:",
  length(
    unique(meta_df$ID)
  ),
  "\n"
)


cat(
  "Unique cows:",
  length(
    unique(meta_df$Sample)
  ),
  "\n"
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
# 5. Read normalized ARG subtype abundance matrix
# ======================================================

feature_df <- read.delim(
  feature_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


colnames(feature_df)[1] <- "feature"


# ======================================================
# 6. Construct sample × ARG subtype matrix
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


colnames(mat) <- make.unique(
  as.character(
    feature_df$feature
  )
)


# ======================================================
# 7. Match samples
# ======================================================

common <- intersect(
  meta_df$ID,
  rownames(mat)
)


cat(
  "Samples matched:",
  length(common),
  "\n"
)


if (length(common) != 120) {

  warning(
    paste0(
      "Expected 120 samples, but ",
      length(common),
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
    "Samples missing from ARG table:\n"
  )

  print(
    missing_samples
  )
}


mat <- mat[
  common,
  ,
  drop = FALSE
]


meta_df <- meta_df[
  match(
    common,
    meta_df$ID
  ),
  ,
  drop = FALSE
]


# ======================================================
# 8. Calculate Shannon diversity
# ======================================================

mat[is.na(mat)] <- 0

mat[mat < 0] <- 0


df_raw <- tibble(

  ID = rownames(mat),

  Shannon = vegan::diversity(
    mat,
    index = "shannon"
  )

) %>%

  left_join(
    meta_df,
    by = "ID"
  )


# ======================================================
# 9. Outlier filtering
# 3 × IQR within each Group × Time
# ======================================================

df_marked <- df_raw %>%

  group_by(
    Group,
    Time
  ) %>%

  mutate(

    Q1 = quantile(
      Shannon,
      0.25,
      na.rm = TRUE
    ),

    Q3 = quantile(
      Shannon,
      0.75,
      na.rm = TRUE
    ),

    IQR_value =
      Q3 -
      Q1,

    lower_3IQR =
      Q1 -
      3 * IQR_value,

    upper_3IQR =
      Q3 +
      3 * IQR_value,

    Outlier =
      Shannon < lower_3IQR |
      Shannon > upper_3IQR

  ) %>%

  ungroup()


# ======================================================
# 10. Export outlier information
# ======================================================

outlier_table <- df_marked %>%
  filter(
    Outlier
  )


write.csv(

  outlier_table,

  file.path(
    output_dir,
    "ARG_Shannon_outliers_3IQR.csv"
  ),

  row.names = FALSE
)


cat(
  "Observations removed by 3 × IQR:",
  nrow(outlier_table),
  "\n"
)


# ======================================================
# 11. Dataset used for analysis
# ======================================================

df <- df_marked %>%

  filter(
    !Outlier
  ) %>%

  select(
    ID,
    Sample,
    Group,
    Time,
    Shannon
  )


df$Group <- factor(
  df$Group,
  levels = c(
    "ITM",
    "OTM1",
    "OTM2"
  )
)


df$Time <- factor(
  df$Time,
  levels = c(
    "D0",
    "D30",
    "D60",
    "D90"
  )
)


# ======================================================
# 12. Colors
# ======================================================

group_cols <- c(
  ITM  = "#8ECFC9",
  OTM1 = "#FFBE7A",
  OTM2 = "#FA7F6F"
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


# ======================================================
# 13. Comparisons
# ======================================================

time_comparisons <- list(
  c("D0", "D30"),
  c("D0", "D60"),
  c("D0", "D90")
)


group_comparisons <- list(
  c("ITM", "OTM1"),
  c("ITM", "OTM2"),
  c("OTM1", "OTM2")
)


# ======================================================
# 14. Plot theme
# ======================================================

theme_main <- theme_classic(
  base_size = 13
) +

  theme(

    legend.position = "none",

    axis.title = element_text(
      face = "bold",
      size = 12
    ),

    axis.title.y = element_text(
      face = "bold",
      size = 10.5
    ),

    axis.title.x = element_text(
      face = "bold",
      size = 11
    ),

    axis.text = element_text(
      color = "black",
      size = 9.5
    ),

    axis.text.y = element_text(
      color = "black",
      size = 8.5
    ),

    axis.text.x = element_text(
      color = "black",
      size = 9.5
    ),

    axis.line = element_line(
      linewidth = 0.55
    ),

    axis.ticks = element_line(
      linewidth = 0.55
    ),

    axis.ticks.length = grid::unit(
      0.12,
      "cm"
    ),

    strip.text = element_text(
      face = "bold",
      size = 10.5
    ),

    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.65
    ),

    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.55
    ),

    panel.spacing = grid::unit(
      0.12,
      "lines"
    ),

    plot.margin = ggplot2::margin(
      5,
      5,
      5,
      5
    )
  )


theme_cd <- theme(

  axis.text.x = element_text(
    angle = 20,
    hjust = 1,
    vjust = 1,
    size = 9
  )
)


# ======================================================
# 15. Boxplot and points
# ======================================================

geom_box_pts <- list(

  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.50
  ),

  geom_point(
    position = position_jitter(
      width = 0.10
    ),
    size = 1.65,
    shape = 21,
    stroke = 0.35,
    color = "black",
    alpha = 0.85
  )
)


# ======================================================
# 16. Within-group temporal significance
#
# Unpaired Wilcoxon tests are intentionally used to
# reproduce the original statistical analysis.
# ======================================================

add_sig_within_time <- function(
    p,
    data,
    yvar
) {

  stat.test <- ggpubr::compare_means(

    as.formula(
      paste(
        yvar,
        "~ Time"
      )
    ),

    data = data,

    group.by = "Group",

    method = "wilcox.test",

    comparisons = time_comparisons

  ) %>%

    rstatix::adjust_pvalue(
      method = "BH"
    ) %>%

    filter(
      p.adj < 0.05
    ) %>%

    mutate(

      p.signif = case_when(

        p.adj < 0.001 ~ "***",

        p.adj < 0.01 ~ "**",

        TRUE ~ "*"
      )
    )


  if (nrow(stat.test) == 0) {
    return(p)
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
        max_val *
        0.035 *
        row_number()
    ) %>%

    ungroup()


  p +

    ggpubr::stat_pvalue_manual(

      stat.test,

      label = "p.signif",

      tip.length = 0.01,

      size = 4.2,

      bracket.size = 0.45,

      step.increase = 0.03
    )
}


# ======================================================
# 17. Between-group significance
# ======================================================

add_sig_between_group <- function(
    p,
    data,
    yvar
) {

  stat.test <- ggpubr::compare_means(

    as.formula(
      paste(
        yvar,
        "~ Group"
      )
    ),

    data = data,

    group.by = "Time",

    method = "wilcox.test",

    comparisons = group_comparisons

  ) %>%

    rstatix::adjust_pvalue(
      method = "BH"
    ) %>%

    filter(
      p.adj < 0.05
    ) %>%

    mutate(

      p.signif = case_when(

        p.adj < 0.001 ~ "***",

        p.adj < 0.01 ~ "**",

        TRUE ~ "*"
      )
    )


  if (nrow(stat.test) == 0) {
    return(p)
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
        max_val *
        0.035 *
        row_number()
    ) %>%

    ungroup()


  p +

    ggpubr::stat_pvalue_manual(

      stat.test,

      label = "p.signif",

      tip.length = 0.01,

      size = 4.2,

      bracket.size = 0.45,

      step.increase = 0.03
    )
}


# ======================================================
# 18. Panel A
# ======================================================

pA_shannon <- ggplot(

  df,

  aes(
    Time,
    Shannon,
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
    y = "Shannon diversity"
  ) +

  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.04,
        0.26
      )
    )
  ) +

  theme_main


pA_shannon <- add_sig_within_time(
  pA_shannon,
  df,
  "Shannon"
)


# ======================================================
# 19. Panel B
# ======================================================

pB_shannon <- ggplot(

  df,

  aes(
    Group,
    Shannon,
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
    y = "Shannon diversity"
  ) +

  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.04,
        0.26
      )
    )
  ) +

  theme_main +

  theme_cd


pB_shannon <- add_sig_between_group(
  pB_shannon,
  df,
  "Shannon"
)


# ======================================================
# 20. Combine panels
# ======================================================

core_shannon <- patchwork::wrap_plots(

  pA_shannon,

  pB_shannon,

  ncol = 2

) +

  patchwork::plot_annotation(
    tag_levels = "A"
  )


# ======================================================
# 21. Legend
# ======================================================

p_leg_shannon <- ggplot(

  df,

  aes(
    Group,
    Shannon,
    fill = Group
  )

) +

  geom_boxplot(
    linewidth = 0.50,
    color = "black"
  ) +

  scale_fill_manual(
    values = group_cols,
    name = "Group"
  ) +

  theme_void() +

  theme(

    legend.position = "top",

    legend.title = element_text(
      face = "bold",
      size = 10
    ),

    legend.text = element_text(
      color = "black",
      size = 9
    )
  )


legend_shannon <- cowplot::get_legend(
  p_leg_shannon
)


final_shannon <- cowplot::ggdraw(
  core_shannon
) +

  cowplot::draw_grob(

    legend_shannon,

    x = 0.55,

    y = 0.93,

    width = 0.40,

    height = 0.08
  )


print(
  final_shannon
)


# ======================================================
# 22. Save figure
# ======================================================

ggsave(

  filename = file.path(
    output_dir,
    "Fig1_ARG_Shannon.pdf"
  ),

  plot = final_shannon,

  width = 14,

  height = 5
)


# ======================================================
# 23. Statistical results
# ======================================================

stat_shannon_within <- ggpubr::compare_means(

  Shannon ~ Time,

  data = df,

  group.by = "Group",

  method = "wilcox.test",

  comparisons = time_comparisons

) %>%

  rstatix::adjust_pvalue(
    method = "BH"
  )


stat_shannon_between <- ggpubr::compare_means(

  Shannon ~ Group,

  data = df,

  group.by = "Time",

  method = "wilcox.test",

  comparisons = group_comparisons

) %>%

  rstatix::adjust_pvalue(
    method = "BH"
  )


# ======================================================
# 24. Export results
# ======================================================

write.csv(

  stat_shannon_within,

  file.path(
    output_dir,
    "ARG_Shannon_withinGroup_overTime_3IQR_BH.csv"
  ),

  row.names = FALSE
)


write.csv(

  stat_shannon_between,

  file.path(
    output_dir,
    "ARG_Shannon_betweenGroup_sameTime_3IQR_BH.csv"
  ),

  row.names = FALSE
)


write.csv(

  df,

  file.path(
    output_dir,
    "ARG_Shannon_data_used.csv"
  ),

  row.names = FALSE
)


cat(
  "\nARG Shannon analysis completed.\n"
)
