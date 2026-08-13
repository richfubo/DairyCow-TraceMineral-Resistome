########################################################
## MRG beta diversity analysis
##
## Figure 2
##
## Analysis level:
##   MRG type/class
##
## Input:
##   metadata/metadata.csv
##
##   processed_data/MRG/
##     MRG_merged_subtype_abundance.tsv
##
## Data processing:
##   MRG merged subtypes are aggregated to the
##   corresponding MRG type/class before analysis.
##
##   Example:
##     Copper__subtype1
##     Copper__subtype2
##          ↓
##        Copper
##
## Analysis:
##   1. Aggregate merged subtypes to MRG types/classes
##   2. Bray-Curtis dissimilarity
##   3. Principal coordinate analysis (PCoA)
##   4. PERMANOVA:
##        Bray-Curtis ~ Group + Time
##      999 permutations
##   5. Group × Time centroid calculation
##   6. Temporal succession trajectories
##
## Output:
##   results/Fig2_MRG_PCoA/
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
  "MRG",
  "MRG_merged_subtype_abundance.tsv"
)


output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig2_MRG_PCoA"
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


# ======================================================
# 4. Prepare metadata
#
# Expected:
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
  "Metadata check\n"
)

cat(
  "========================================\n"
)


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
  "\nSamples by Group × Time:\n"
)


print(
  table(
    meta_df$Group,
    meta_df$Time
  )
)


# ======================================================
# 6. Read MRG merged-subtype abundance matrix
#
# Expected format:
#
# First column:
#   merged MRG subtype
#
# Remaining columns:
#   sequencing sample IDs
#
# Example feature names:
#
# Copper__xxx
# Copper__yyy
# Zinc__xxx
# Iron__xxx
#
# MRG type/class is extracted from the part before "__".
# ======================================================

feature_df <- read.delim(
  feature_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


colnames(feature_df)[1] <- "feature"


feature_df$feature <- as.character(
  feature_df$feature
)


sample_cols <- colnames(
  feature_df
)[-1]


# ======================================================
# 7. Convert abundance columns to numeric
# ======================================================

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


    # Missing values -> zero
    v[
      is.na(v)
    ] <- 0


    # Negative abundance values are not allowed
    v[
      v < 0
    ] <- 0


    return(
      v
    )
  }
)


# ======================================================
# 8. Extract MRG type/class from merged-subtype name
#
# Rule:
#
# Copper__xxx  -> Copper
# Zinc__xxx    -> Zinc
# Iron__xxx    -> Iron
#
# If "__" is absent, the original feature name is kept.
# ======================================================

feature_df <- feature_df %>%

  mutate(

    MRG_type = if_else(

      stringr::str_detect(
        feature,
        fixed("__")
      ),

      stringr::str_replace(
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
  "MRG merged-subtype → type aggregation\n"
)

cat(
  "========================================\n"
)


cat(
  "Number of merged subtypes:",
  nrow(feature_df),
  "\n"
)


cat(
  "Number of MRG types/classes:",
  length(
    unique(feature_df$MRG_type)
  ),
  "\n"
)


cat(
  "\nMRG types/classes detected:\n"
)


print(
  sort(
    unique(
      feature_df$MRG_type
    )
  )
)


# ======================================================
# 9. Export merged-subtype → type mapping
#
# This makes the aggregation step fully reproducible.
# ======================================================

mapping_df <- feature_df %>%

  select(
    feature,
    MRG_type
  ) %>%

  distinct()


write.csv(

  mapping_df,

  file.path(
    output_dir,
    "MRG_merged_subtype_to_type_mapping.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 10. Aggregate merged subtypes to MRG type/class
#
# Abundances belonging to the same MRG type are summed
# within each sequencing sample.
# ======================================================

mrg_type_df <- feature_df %>%

  select(
    MRG_type,
    all_of(sample_cols)
  ) %>%

  group_by(
    MRG_type
  ) %>%

  summarise(

    across(
      all_of(sample_cols),
      ~ sum(
        .x,
        na.rm = TRUE
      )
    ),

    .groups = "drop"
  )


cat(
  "\nAggregated MRG type abundance matrix:\n"
)


cat(
  nrow(mrg_type_df),
  "MRG types ×",
  length(sample_cols),
  "samples\n"
)


# ======================================================
# 11. Export derived MRG type abundance table
#
# This is a derived output rather than an additional
# required input file.
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


# ======================================================
# 12. Construct sample × MRG type matrix
# ======================================================

x <- mrg_type_df[
  ,
  -1,
  drop = FALSE
]


mat <- t(
  as.matrix(x)
)


rownames(mat) <- colnames(
  mrg_type_df
)[-1]


colnames(mat) <- make.unique(
  as.character(
    mrg_type_df$MRG_type
  )
)


# Replace missing values with zero
mat[
  is.na(mat)
] <- 0


# Negative abundance values are not allowed
mat[
  mat < 0
] <- 0


# ======================================================
# 13. Match metadata and abundance samples
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


if (length(common_samples) < 5) {

  stop(
    paste0(
      "Too few common samples. ",
      "Please check sample IDs in metadata.csv and ",
      "MRG_merged_subtype_abundance.tsv."
    )
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
    meta_df$ID ==
      rownames(mat)
  )
)


# ======================================================
# 14. Remove MRG types absent from all matched samples
# ======================================================

mat <- mat[
  ,
  colSums(
    mat,
    na.rm = TRUE
  ) > 0,
  drop = FALSE
]


cat(
  "\nNumber of samples:",
  nrow(mat),
  "\n"
)


cat(
  "Number of MRG types/classes:",
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
# 15. Check zero-abundance samples
# ======================================================

zero_samples <- rownames(
  mat
)[
  rowSums(
    mat,
    na.rm = TRUE
  ) == 0
]


if (length(zero_samples) > 0) {

  stop(
    paste0(
      "Zero-abundance MRG profiles detected for sample(s): ",
      paste(
        zero_samples,
        collapse = ", "
      )
    )
  )
}


# ======================================================
# 16. Bray-Curtis dissimilarity
# ======================================================

bray_dist <- vegan::vegdist(
  mat,
  method = "bray"
)


# ======================================================
# 17. Principal coordinate analysis (PCoA)
# ======================================================

pcoa_res <- ape::pcoa(
  bray_dist
)


pcoa_df <- data.frame(

  ID = rownames(
    mat
  ),

  PCoA1 =
    pcoa_res$vectors[
      ,
      1
    ],

  PCoA2 =
    pcoa_res$vectors[
      ,
      2
    ]

) %>%

  left_join(
    meta_df,
    by = "ID"
  )


# Variation explained by first two axes
eig_percent <- round(

  100 *

    pcoa_res$values$Relative_eig[
      1:2
    ],

  1
)


# ======================================================
# 18. Calculate Group × Time centroids
# ======================================================

centroid_df <- pcoa_df %>%

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
  ) %>%

  arrange(
    Group,
    Time
  )


# ======================================================
# 19. PERMANOVA
#
# Model:
# Bray-Curtis ~ Group + Time
#
# 999 permutations
#
# The same statistical model as the original analysis
# is retained.
# ======================================================

set.seed(
  123
)


permanova_res <- vegan::adonis2(

  bray_dist ~ Group + Time,

  data = meta_df,

  permutations = 999,

  by = "terms"
)


print(
  permanova_res
)


permanova_df <- as.data.frame(
  permanova_res
) %>%

  tibble::rownames_to_column(
    "Factor"
  )


print(
  permanova_df
)


# ======================================================
# 20. Export PERMANOVA results
# ======================================================

write.csv(

  permanova_df,

  file.path(
    output_dir,
    "stat_MRG_type_PERMANOVA_BrayCurtis_Group_Time.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 21. Extract Group and Time R2 and P values
# ======================================================

group_R2 <- permanova_df %>%

  filter(
    Factor == "Group"
  ) %>%

  pull(
    R2
  )


group_P <- permanova_df %>%

  filter(
    Factor == "Group"
  ) %>%

  pull(
    `Pr(>F)`
  )


time_R2 <- permanova_df %>%

  filter(
    Factor == "Time"
  ) %>%

  pull(
    R2
  )


time_P <- permanova_df %>%

  filter(
    Factor == "Time"
  ) %>%

  pull(
    `Pr(>F)`
  )


format_p <- function(p) {

  ifelse(

    is.na(p),

    "NA",

    ifelse(

      p < 0.001,

      "< 0.001",

      sprintf(
        "%.3f",
        p
      )
    )
  )
}


permanova_label <- paste0(

  "PERMANOVA (Bray\u2212Curtis):\n",

  "Group: R\u00b2 = ",
  sprintf(
    "%.3f",
    group_R2
  ),

  ", P = ",
  format_p(
    group_P
  ),

  "\n",

  "Time: R\u00b2 = ",
  sprintf(
    "%.3f",
    time_R2
  ),

  ", P = ",
  format_p(
    time_P
  )
)


cat(
  "\nPERMANOVA annotation:\n"
)


cat(
  permanova_label,
  "\n"
)


# ======================================================
# 22. Colors
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
# 23. Determine PERMANOVA annotation position
# ======================================================

x_range <- range(
  pcoa_df$PCoA1,
  na.rm = TRUE
)


y_range <- range(
  pcoa_df$PCoA2,
  na.rm = TRUE
)


label_x <-
  x_range[1] +
  0.03 *
  diff(
    x_range
  )


label_y <-
  y_range[1] +
  0.08 *
  diff(
    y_range
  )


# ======================================================
# 24. Plot PCoA succession trajectory
# ======================================================

p_pcoa <- ggplot(

  pcoa_df,

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

    size = 2.2,

    alpha = 0.35
  ) +


  # Group centroid trajectories over time
  geom_path(

    data = centroid_df,

    aes(
      x = PCoA1,
      y = PCoA2,
      color = Group,
      group = Group
    ),

    linewidth = 1.2,

    arrow = arrow(

      length = grid::unit(
        0.18,
        "cm"
      ),

      type = "closed"
    )
  ) +


  # Group × Time centroids
  geom_point(

    data = centroid_df,

    aes(
      x = PCoA1,
      y = PCoA2,
      color = Group
    ),

    size = 5,

    alpha = 0.95
  ) +


  # Sampling-time labels
  ggrepel::geom_text_repel(

    data = centroid_df,

    aes(
      x = PCoA1,
      y = PCoA2,
      label = Time
    ),

    color = "black",

    size = 4,

    fontface = "bold",

    show.legend = FALSE,

    max.overlaps = Inf,

    box.padding = 0.25,

    point.padding = 0.25,

    segment.color = NA
  ) +


  # PERMANOVA annotation
  annotate(

    "text",

    x = label_x,

    y = label_y,

    label = permanova_label,

    hjust = 0,

    vjust = 0,

    size = 4.2,

    fontface = "italic"
  ) +


  scale_color_manual(
    values = group_cols
  ) +


  labs(

    x = paste0(
      "PCoA1 (",
      eig_percent[1],
      "%)"
    ),

    y = paste0(
      "PCoA2 (",
      eig_percent[2],
      "%)"
    ),

    color = NULL,

    title =
      "Succession trajectory of metal resistome profile"
  ) +


  theme_bw(
    base_size = 14
  ) +


  theme(

    panel.grid.major =
      element_blank(),

    panel.grid.minor =
      element_blank(),


    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5,
        size = 18
      ),


    axis.title =
      element_text(
        face = "bold",
        size = 14
      ),


    axis.text =
      element_text(
        color = "black",
        size = 12
      ),


    axis.line =
      element_line(
        linewidth = 0.8,
        color = "black"
      ),


    axis.ticks =
      element_line(
        linewidth = 0.8,
        color = "black"
      ),


    legend.position =
      "right",

    legend.text =
      element_text(
        size = 12
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
  p_pcoa
)


# ======================================================
# 25. Save PCoA figure
# ======================================================

ggsave(

  filename = file.path(
    output_dir,
    "Figure_MRG_type_PCoA_BrayCurtis_Group_Time.pdf"
  ),

  plot = p_pcoa,

  width = 8,

  height = 6
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_MRG_type_PCoA_BrayCurtis_Group_Time.png"
  ),

  plot = p_pcoa,

  width = 8,

  height = 6,

  dpi = 600
)


# ======================================================
# 26. Export PCoA coordinates
# ======================================================

write.csv(

  pcoa_df,

  file.path(
    output_dir,
    "stat_MRG_type_PCoA_coordinates.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 27. Export Group × Time centroids
# ======================================================

write.csv(

  centroid_df,

  file.path(
    output_dir,
    "stat_MRG_type_PCoA_centroids.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 28. Finish
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "MRG type-level PCoA analysis completed.\n"
)

cat(
  "Core input:\n"
)

cat(
  feature_path,
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