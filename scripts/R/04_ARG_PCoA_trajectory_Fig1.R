########################################################
## ARG resistome PCoA succession trajectory
##
## Figure 1
##
## Input:
##   metadata/metadata.csv
##   processed_data/ARG/normalized_cell.subtype_matrix.tsv
##
## Analysis:
##   1. Bray-Curtis dissimilarity
##   2. Principal coordinate analysis (PCoA)
##   3. PERMANOVA: Group + Time
##      999 permutations
##   4. Group × Time centroid calculation
##   5. D0-D30-D60-D90 succession trajectories
##
## Output:
##   results/Fig1_ARG_PCoA/
########################################################


# ======================================================
# 1. Load packages
# ======================================================

library(tidyverse)
library(vegan)
library(ape)
library(ggrepel)
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


out_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig1_ARG_PCoA"
)


dir.create(
  out_dir,
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
  "Group"
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
# 4. Read normalized ARG subtype abundance matrix
# ======================================================

feature_df <- read.delim(
  feature_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


colnames(feature_df)[1] <- "feature"


# ======================================================
# 5. Prepare metadata
#
# Expected metadata:
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

    Group = as.character(Group),

    Time = as.character(Time)
  )


# Check treatment names
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


# Set factor order
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
# 6. Basic metadata checks
# ======================================================

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
  "Number of unique sequencing samples:",
  length(
    unique(meta_df$ID)
  ),
  "\n\n"
)


cat(
  "Samples by Group × Time:\n"
)


print(
  table(
    meta_df$Group,
    meta_df$Time
  )
)


cat("\n")


if (anyDuplicated(meta_df$ID) > 0) {

  stop(
    "Duplicated sequencing sample IDs were found."
  )
}


# ======================================================
# 7. Construct sample × ARG subtype matrix
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


# Missing values to zero
mat[is.na(mat)] <- 0


# Negative abundance values are not allowed
mat[mat < 0] <- 0


# Remove ARG subtypes with zero abundance
# across all samples
mat <- mat[
  ,
  colSums(
    mat,
    na.rm = TRUE
  ) > 0,
  drop = FALSE
]


# ======================================================
# 8. Match samples between metadata and ARG table
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


missing_samples <- setdiff(
  meta_df$ID,
  rownames(mat)
)


if (length(missing_samples) > 0) {

  cat(
    "\nSamples missing from ARG matrix:\n"
  )

  print(
    missing_samples
  )
}


if (length(common_samples) < 5) {

  stop(
    paste0(
      "Too few common samples. ",
      "Please check sample IDs in metadata.csv ",
      "and normalized_cell.subtype_matrix.tsv."
    )
  )
}


# Align matrix and metadata
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


cat(
  "\nNumber of samples used:",
  nrow(mat),
  "\n"
)


cat(
  "Number of ARG subtypes:",
  ncol(mat),
  "\n"
)


cat(
  "\nGroup × Time sample distribution:\n"
)


print(
  table(
    meta_df$Group,
    meta_df$Time
  )
)


# ======================================================
# 9. Bray-Curtis dissimilarity
# ======================================================

bray_dist <- vegan::vegdist(
  mat,
  method = "bray"
)


# ======================================================
# 10. Principal coordinate analysis (PCoA)
# ======================================================

pcoa_res <- ape::pcoa(
  bray_dist
)


pcoa_scores <- as.data.frame(
  pcoa_res$vectors[
    ,
    1:2
  ]
)


colnames(
  pcoa_scores
) <- c(
  "PCoA1",
  "PCoA2"
)


pcoa_scores$ID <- rownames(
  pcoa_scores
)


# Percentage of variation explained
eig <- round(
  100 *
    pcoa_res$values$Relative_eig[
      1:2
    ],
  1
)


# Merge metadata
plot_df <- pcoa_scores %>%

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
# 11. PERMANOVA
#
# Model:
# Bray-Curtis ~ Group + Time
#
# 999 permutations
# ======================================================

set.seed(
  123
)


adonis_res <- vegan::adonis2(

  bray_dist ~ Group + Time,

  data = meta_df,

  permutations = 999,

  by = "terms"
)


print(
  adonis_res
)


# Convert PERMANOVA results to table
adonis_df <- as.data.frame(
  adonis_res
) %>%

  tibble::rownames_to_column(
    "Term"
  ) %>%

  as_tibble()


write.csv(

  adonis_df,

  file.path(
    out_dir,
    "stat_ARG_PCoA_PERMANOVA.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 12. Extract PERMANOVA statistics for figure
# ======================================================

get_adonis_text <- function(
    adonis_table,
    term_name
) {

  row <- adonis_table %>%
    filter(
      Term == term_name
    )


  if (nrow(row) == 0) {

    return(
      NULL
    )
  }


  p_text <- ifelse(

    row$`Pr(>F)` < 0.001,

    "< 0.001",

    sprintf(
      "%.3f",
      row$`Pr(>F)`
    )
  )


  paste0(

    term_name,

    ": R² = ",

    sprintf(
      "%.3f",
      row$R2
    ),

    ", P ",

    ifelse(
      row$`Pr(>F)` < 0.001,
      "=",
      "="
    ),

    " ",

    p_text
  )
}


group_text <- get_adonis_text(
  adonis_df,
  "Group"
)


time_text <- get_adonis_text(
  adonis_df,
  "Time"
)


adonis_label <- paste0(

  "PERMANOVA (Bray-Curtis):\n",

  group_text,

  "\n",

  time_text
)


# ======================================================
# 13. Calculate Group × Time centroids
# ======================================================

centroid_df <- plot_df %>%

  filter(
    !is.na(Group),
    !is.na(Time)
  ) %>%

  group_by(
    Group,
    Time
  ) %>%

  summarise(

    PCoA1 = mean(
      PCoA1,
      na.rm = TRUE
    ),

    PCoA2 = mean(
      PCoA2,
      na.rm = TRUE
    ),

    n = n(),

    .groups = "drop"
  )


# Add temporal order
centroid_df <- centroid_df %>%

  mutate(
    Time_num = as.numeric(
      Time
    )
  ) %>%

  arrange(
    Group,
    Time_num
  )


# ======================================================
# 14. Colors
# ======================================================

group_cols <- c(

  ITM  = "#8ECFC9",

  OTM1 = "#FFBE7A",

  OTM2 = "#FA7F6F"
)


# ======================================================
# 15. Plot theme
# ======================================================

theme_pcoa <- theme_classic(
  base_size = 13
) +

  theme(

    legend.position = "right",

    axis.title = element_text(
      face = "bold",
      size = 12
    ),

    axis.text = element_text(
      color = "black",
      size = 10
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

    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 15
    ),

    legend.title = element_text(
      face = "bold",
      size = 10.5
    ),

    legend.text = element_text(
      color = "black",
      size = 9.5
    ),

    plot.margin = ggplot2::margin(
      5,
      8,
      5,
      5
    )
  )


# ======================================================
# 16. Determine PERMANOVA annotation position
# ======================================================

x_range <- range(
  plot_df$PCoA1,
  na.rm = TRUE
)


y_range <- range(
  plot_df$PCoA2,
  na.rm = TRUE
)


label_x <-
  x_range[1] +
  0.03 *
  diff(x_range)


label_y <-
  y_range[1] +
  0.08 *
  diff(y_range)


# ======================================================
# 17. Plot PCoA succession trajectory
# ======================================================

p_pcoa_traj <- ggplot(

  plot_df,

  aes(
    x = PCoA1,
    y = PCoA2
  )

) +


  # Individual samples
  geom_point(

    aes(
      color = Group
    ),

    size = 1.8,

    alpha = 0.35
  ) +


  # Group × Time centroid trajectories
  geom_path(

    data = centroid_df,

    aes(
      x = PCoA1,
      y = PCoA2,
      color = Group,
      group = Group
    ),

    linewidth = 0.9,

    alpha = 0.95,

    arrow = arrow(

      length = grid::unit(
        0.12,
        "cm"
      ),

      type = "closed"
    )
  ) +


  # Centroid points
  geom_point(

    data = centroid_df,

    aes(
      x = PCoA1,
      y = PCoA2,
      fill = Group
    ),

    shape = 21,

    size = 3.2,

    stroke = 0.55,

    color = "black"
  ) +


  # Sampling-time labels
  ggrepel::geom_text_repel(

    data = centroid_df,

    aes(
      x = PCoA1,
      y = PCoA2,
      label = Time,
      color = Group
    ),

    size = 3.5,

    fontface = "bold",

    show.legend = FALSE,

    max.overlaps = Inf,

    box.padding = 0.25,

    point.padding = 0.25,

    segment.size = 0.25,

    min.segment.length = 0
  ) +


  # PERMANOVA annotation
  annotate(

    "text",

    x = label_x,

    y = label_y,

    label = adonis_label,

    hjust = 0,

    vjust = 0,

    size = 4,

    fontface = "italic"
  ) +


  scale_color_manual(

    values = group_cols,

    name = "Group"
  ) +


  scale_fill_manual(

    values = group_cols,

    name = "Group"
  ) +


  labs(

    title =
      "Succession trajectory of ARG resistome profile",

    x = paste0(
      "PCoA1 (",
      eig[1],
      "%)"
    ),

    y = paste0(
      "PCoA2 (",
      eig[2],
      "%)"
    )
  ) +


  theme_pcoa


print(
  p_pcoa_traj
)


# ======================================================
# 18. Save figure
# ======================================================

ggsave(

  filename = file.path(
    out_dir,
    "Figure_ARG_PCoA_Trajectory.pdf"
  ),

  plot = p_pcoa_traj,

  width = 8.5,

  height = 5.2
)


ggsave(

  filename = file.path(
    out_dir,
    "Figure_ARG_PCoA_Trajectory.png"
  ),

  plot = p_pcoa_traj,

  width = 8.5,

  height = 5.2,

  dpi = 600
)


# ======================================================
# 19. Export PCoA coordinates
# ======================================================

write.csv(

  plot_df,

  file.path(
    out_dir,
    "stat_ARG_PCoA_coordinates.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 20. Export trajectory centroids
# ======================================================

write.csv(

  centroid_df,

  file.path(
    out_dir,
    "stat_ARG_PCoA_centroids.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 21. Finish
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "ARG resistome PCoA analysis completed.\n"
)

cat(
  "Output directory:\n"
)

cat(
  out_dir,
  "\n"
)

cat(
  "========================================\n"
)