########################################################
## ARG total abundance analysis
##
## Figure 1
##
## Panel A:
## Within-treatment temporal comparisons
## D0 vs D30, D0 vs D60, D0 vs D90
## Paired Wilcoxon signed-rank test
## Pairing variable: Sample (cow ID)
##
## Panel B:
## Between-treatment comparisons at each time point
## Wilcoxon rank-sum test
##
## Multiple-testing correction:
## Benjamini-Hochberg (BH/FDR)
##
## Outlier filtering:
## 3 × IQR within each Group × Time
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

feature_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "ARG",
  "normalized_cell.subtype_matrix.tsv"
)

output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig1_ARG_total_abundance"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ======================================================
# 3. Read metadata
#
# Expected columns:
# ID | Time | Group | Sample
#
# ID     = sequencing sample ID (a1-a120)
# Time   = D0, D30, D60, D90
# Group  = ITM, OTM1, OTM2
# Sample = cow ID
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
      paste(missing_cols, collapse = ", ")
    )
  )
}

meta_df <- meta_df %>%
  mutate(
    ID = as.character(ID),
    Sample = as.character(Sample),
    Group = as.character(Group),
    Time = as.character(Time)
  )

meta_df$Group <- factor(
  meta_df$Group,
  levels = c("ITM", "OTM1", "OTM2")
)

meta_df$Time <- factor(
  meta_df$Time,
  levels = c("D0", "D30", "D60", "D90")
)


# ======================================================
# 4. Check metadata
# ======================================================

cat("Number of metadata rows:", nrow(meta_df), "\n")
cat("Number of unique sample IDs:", length(unique(meta_df$ID)), "\n")
cat("Number of unique cows:", length(unique(meta_df$Sample)), "\n")

print(
  table(
    meta_df$Group,
    meta_df$Time
  )
)

if (anyDuplicated(meta_df$ID) > 0) {
  stop("Duplicated sample IDs were found in metadata.")
}


# ======================================================
# 5. Read normalized ARG subtype abundance matrix
#
# First column:
# ARG subtype
#
# Remaining columns:
# a1-a120
# ======================================================

feature_df <- read.delim(
  feature_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(feature_df)[1] <- "feature"


# ======================================================
# 6. Convert to sample × ARG matrix
# ======================================================

x <- feature_df[, -1, drop = FALSE]

x[] <- lapply(
  x,
  function(v) as.numeric(as.character(v))
)

mat <- t(as.matrix(x))

rownames(mat) <- colnames(feature_df)[-1]

colnames(mat) <- make.unique(
  as.character(feature_df$feature)
)


# ======================================================
# 7. Match samples with metadata
# ======================================================

common_samples <- intersect(
  meta_df$ID,
  rownames(mat)
)

cat(
  "Samples matched between metadata and ARG table:",
  length(common_samples),
  "\n"
)

if (length(common_samples) != 120) {
  warning(
    paste0(
      "Expected 120 samples, but ",
      length(common_samples),
      " were matched."
    )
  )
}

missing_in_ARG <- setdiff(
  meta_df$ID,
  rownames(mat)
)

if (length(missing_in_ARG) > 0) {
  cat("Samples missing from ARG table:\n")
  print(missing_in_ARG)
}

mat <- mat[
  common_samples,
  ,
  drop = FALSE
]

meta_df <- meta_df[
  match(common_samples, meta_df$ID),
  ,
  drop = FALSE
]


# ======================================================
# 8. Calculate total normalized ARG abundance
# ======================================================

mat[is.na(mat)] <- 0
mat[mat < 0] <- 0

df_raw <- tibble(
  ID = rownames(mat),
  TotalAbundance = rowSums(
    mat,
    na.rm = TRUE
  )
) %>%
  left_join(
    meta_df,
    by = "ID"
  )


# ======================================================
# 9. Outlier filtering: 3 × IQR
#
# Applied independently within each Group × Time
# ======================================================

df_marked <- df_raw %>%
  group_by(
    Group,
    Time
  ) %>%
  mutate(
    Q1 = quantile(
      TotalAbundance,
      0.25,
      na.rm = TRUE
    ),

    Q3 = quantile(
      TotalAbundance,
      0.75,
      na.rm = TRUE
    ),

    IQR_value = Q3 - Q1,

    lower_limit = if_else(
      IQR_value == 0,
      -Inf,
      Q1 - 3 * IQR_value
    ),

    upper_limit = if_else(
      IQR_value == 0,
      Inf,
      Q3 + 3 * IQR_value
    ),

    Outlier =
      TotalAbundance < lower_limit |
      TotalAbundance > upper_limit
  ) %>%
  ungroup()


# Export excluded observations
outlier_table <- df_marked %>%
  filter(Outlier)

write.csv(
  outlier_table,
  file.path(
    output_dir,
    "ARG_total_abundance_outliers_3IQR.csv"
  ),
  row.names = FALSE
)

cat(
  "Number of observations removed by 3 × IQR:",
  nrow(outlier_table),
  "\n"
)


# Dataset used in subsequent analyses
df <- df_marked %>%
  filter(!Outlier) %>%
  select(
    ID,
    Sample,
    Group,
    Time,
    TotalAbundance
  )

df$Group <- factor(
  df$Group,
  levels = c("ITM", "OTM1", "OTM2")
)

df$Time <- factor(
  df$Time,
  levels = c("D0", "D30", "D60", "D90")
)


# ======================================================
# 10. Colors
# ======================================================

group_cols <- c(
  ITM  = "#8ECFC9",
  OTM1 = "#FFBE7A",
  OTM2 = "#FA7F6F"
)

make_group_gradient <- function(base) {
  c(
    colorspace::lighten(base, 0.65),
    colorspace::lighten(base, 0.40),
    colorspace::lighten(base, 0.20),
    base
  )
}

grad_cols <- c(
  setNames(
    make_group_gradient(group_cols["ITM"]),
    paste0("ITM.", levels(df$Time))
  ),

  setNames(
    make_group_gradient(group_cols["OTM1"]),
    paste0("OTM1.", levels(df$Time))
  ),

  setNames(
    make_group_gradient(group_cols["OTM2"]),
    paste0("OTM2.", levels(df$Time))
  )
)

df$GroupTime <- interaction(
  df$Group,
  df$Time,
  sep = "."
)


# ======================================================
# 11. Comparisons
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
# 12. Paired Wilcoxon:
# Within-treatment temporal comparisons
# ======================================================

paired_wilcox_one <- function(
    data,
    group_name,
    time1,
    time2
) {

  d1 <- data %>%
    filter(
      Group == group_name,
      Time == time1
    ) %>%
    select(
      Sample,
      value1 = TotalAbundance
    )

  d2 <- data %>%
    filter(
      Group == group_name,
      Time == time2
    ) %>%
    select(
      Sample,
      value2 = TotalAbundance
    )

  paired_data <- inner_join(
    d1,
    d2,
    by = "Sample"
  )

  if (nrow(paired_data) < 2) {
    p_value <- NA_real_
  } else {
    p_value <- wilcox.test(
      paired_data$value1,
      paired_data$value2,
      paired = TRUE,
      exact = FALSE
    )$p.value
  }

  tibble(
    Group = group_name,
    group1 = time1,
    group2 = time2,
    n_pairs = nrow(paired_data),
    p = p_value
  )
}


stat_total_within <- purrr::map_dfr(
  levels(df$Group),

  function(g) {

    purrr::map_dfr(
      time_comparisons,

      function(comp) {

        paired_wilcox_one(
          data = df,
          group_name = g,
          time1 = comp[1],
          time2 = comp[2]
        )
      }
    )
  }
)


stat_total_within$p.adj <- p.adjust(
  stat_total_within$p,
  method = "BH"
)

stat_total_within <- stat_total_within %>%
  mutate(
    p.signif = case_when(
      is.na(p.adj) ~ NA_character_,
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )


# ======================================================
# 13. Unpaired Wilcoxon:
# Between treatments at each time point
# ======================================================

between_wilcox_one <- function(
    data,
    time_name,
    group1_name,
    group2_name
) {

  x <- data %>%
    filter(
      Time == time_name,
      Group == group1_name
    ) %>%
    pull(TotalAbundance)

  y <- data %>%
    filter(
      Time == time_name,
      Group == group2_name
    ) %>%
    pull(TotalAbundance)

  if (
    length(x) < 2 ||
    length(y) < 2
  ) {
    p_value <- NA_real_
  } else {
    p_value <- wilcox.test(
      x,
      y,
      paired = FALSE,
      exact = FALSE
    )$p.value
  }

  tibble(
    Time = time_name,
    group1 = group1_name,
    group2 = group2_name,
    n1 = length(x),
    n2 = length(y),
    p = p_value
  )
}


stat_total_between <- purrr::map_dfr(
  levels(df$Time),

  function(tt) {

    purrr::map_dfr(
      group_comparisons,

      function(comp) {

        between_wilcox_one(
          data = df,
          time_name = tt,
          group1_name = comp[1],
          group2_name = comp[2]
        )
      }
    )
  }
)


stat_total_between$p.adj <- p.adjust(
  stat_total_between$p,
  method = "BH"
)

stat_total_between <- stat_total_between %>%
  mutate(
    p.signif = case_when(
      is.na(p.adj) ~ NA_character_,
      p.adj < 0.001 ~ "***",
      p.adj < 0.01 ~ "**",
      p.adj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )


# ======================================================
# 14. Significance-bar positions
# ======================================================

within_sig <- stat_total_within %>%
  filter(
    !is.na(p.adj),
    p.adj < 0.05
  )

if (nrow(within_sig) > 0) {

  within_y <- df %>%
    group_by(Group) %>%
    summarise(
      max_val = max(
        TotalAbundance,
        na.rm = TRUE
      ),

      min_val = min(
        TotalAbundance,
        na.rm = TRUE
      ),

      .groups = "drop"
    ) %>%
    mutate(
      y_range = max_val - min_val,

      y_range = if_else(
        y_range == 0,
        pmax(max_val, 1),
        y_range
      )
    )

  within_sig <- within_sig %>%
    left_join(
      within_y,
      by = "Group"
    ) %>%
    group_by(Group) %>%
    mutate(
      y.position =
        max_val +
        y_range *
        0.08 *
        row_number()
    ) %>%
    ungroup()
}


between_sig <- stat_total_between %>%
  filter(
    !is.na(p.adj),
    p.adj < 0.05
  )

if (nrow(between_sig) > 0) {

  between_y <- df %>%
    group_by(Time) %>%
    summarise(
      max_val = max(
        TotalAbundance,
        na.rm = TRUE
      ),

      min_val = min(
        TotalAbundance,
        na.rm = TRUE
      ),

      .groups = "drop"
    ) %>%
    mutate(
      y_range = max_val - min_val,

      y_range = if_else(
        y_range == 0,
        pmax(max_val, 1),
        y_range
      )
    )

  between_sig <- between_sig %>%
    left_join(
      between_y,
      by = "Time"
    ) %>%
    group_by(Time) %>%
    mutate(
      y.position =
        max_val +
        y_range *
        0.08 *
        row_number()
    ) %>%
    ungroup()
}


# ======================================================
# 15. Plot theme
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

    axis.text = element_text(
      color = "black"
    ),

    axis.text.y = element_text(
      size = 8.5
    ),

    axis.text.x = element_text(
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
    )
  )


# ======================================================
# 16. Panel A:
# Within-treatment temporal comparison
# ======================================================

pA_total <- ggplot(
  df,
  aes(
    x = Time,
    y = TotalAbundance,
    fill = GroupTime
  )
) +

  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.50
  ) +

  geom_point(
    position = position_jitter(
      width = 0.10
    ),
    size = 1.65,
    shape = 21,
    stroke = 0.35,
    color = "black",
    alpha = 0.85
  ) +

  facet_wrap(
    ~ Group,
    nrow = 1
  ) +

  scale_fill_manual(
    values = grad_cols
  ) +

  labs(
    x = NULL,
    y = "Total ARG abundance (normalized)"
  ) +

  scale_y_continuous(
    expand = expansion(
      mult = c(0.02, 0.15)
    )
  ) +

  theme_main


if (nrow(within_sig) > 0) {

  pA_total <- pA_total +

    ggpubr::stat_pvalue_manual(
      within_sig,
      label = "p.signif",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      size = 4.2,
      bracket.size = 0.45
    )
}


# ======================================================
# 17. Panel B:
# Between-treatment comparison
# ======================================================

pB_total <- ggplot(
  df,
  aes(
    x = Group,
    y = TotalAbundance,
    fill = Group
  )
) +

  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.50
  ) +

  geom_point(
    position = position_jitter(
      width = 0.10
    ),
    size = 1.65,
    shape = 21,
    stroke = 0.35,
    color = "black",
    alpha = 0.85
  ) +

  facet_wrap(
    ~ Time,
    nrow = 1
  ) +

  scale_fill_manual(
    values = group_cols
  ) +

  labs(
    x = NULL,
    y = "Total ARG abundance (normalized)"
  ) +

  scale_y_continuous(
    expand = expansion(
      mult = c(0.02, 0.15)
    )
  ) +

  theme_main +

  theme(
    axis.text.x = element_text(
      angle = 20,
      hjust = 1,
      vjust = 1,
      size = 9
    )
  )


if (nrow(between_sig) > 0) {

  pB_total <- pB_total +

    ggpubr::stat_pvalue_manual(
      between_sig,
      label = "p.signif",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      size = 4.2,
      bracket.size = 0.45
    )
}


# ======================================================
# 18. Combine panels
# ======================================================

core_total <- patchwork::wrap_plots(
  pA_total,
  pB_total,
  ncol = 2
) +
  patchwork::plot_annotation(
    tag_levels = "A"
  )


# ======================================================
# 19. Group legend
# ======================================================

p_leg_total <- ggplot(
  df,
  aes(
    Group,
    TotalAbundance,
    fill = Group
  )
) +
  geom_boxplot() +
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
      size = 9
    )
  )

legend_total <- cowplot::get_legend(
  p_leg_total
)

final_total <- cowplot::ggdraw(
  core_total
) +
  cowplot::draw_grob(
    legend_total,
    x = 0.55,
    y = 0.93,
    width = 0.40,
    height = 0.08
  )

print(final_total)


# ======================================================
# 20. Save outputs
# ======================================================

ggsave(
  filename = file.path(
    output_dir,
    "Fig1_ARG_total_abundance.pdf"
  ),
  plot = final_total,
  width = 14,
  height = 5
)

write.csv(
  stat_total_within,
  file.path(
    output_dir,
    "ARG_total_within_group_paired_Wilcoxon_BH.csv"
  ),
  row.names = FALSE
)

write.csv(
  stat_total_between,
  file.path(
    output_dir,
    "ARG_total_between_group_Wilcoxon_BH.csv"
  ),
  row.names = FALSE
)

write.csv(
  df,
  file.path(
    output_dir,
    "ARG_total_abundance_data_used.csv"
  ),
  row.names = FALSE
)

cat(
  "\nARG total abundance analysis completed.\n"
)
