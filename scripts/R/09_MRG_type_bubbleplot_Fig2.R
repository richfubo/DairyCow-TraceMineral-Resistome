########################################################
## MRG type/class bubble plot
##
## Figure 2
##
## Input:
##   metadata/metadata.csv
##
##   processed_data/MRG/
##     MRG_merged_subtype_abundance.tsv
##
## Data processing:
##   MRG merged subtypes are aggregated to the
##   corresponding MRG type/class before plotting.
##
##   Example:
##     Copper__subtype1
##     Copper__subtype2
##          ↓
##        Copper
##
## Plot:
##   Facets: D0 / D30 / D60 / D90
##   x-axis: ITM / OTM1 / OTM2
##   y-axis: MRG type/class
##
## Bubble size:
##   Median cell number-normalized abundance
##
## Bubble color:
##   log10(median abundance + 1e-6)
##
## Output:
##   results/Fig2_MRG_type_bubbleplot/
########################################################


# ======================================================
# 1. Load packages
# ======================================================

library(tidyverse)


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
  "Fig2_MRG_type_bubbleplot"
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
# First column:
#   merged MRG subtype
#
# Remaining columns:
#   sequencing sample IDs
#
# Example:
#   Copper__xxx
#   Copper__yyy
#   Zinc__xxx
#   Iron__xxx
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
# 8. Extract MRG type/class
#
# Rule:
# Copper__xxx -> Copper
# Zinc__xxx   -> Zinc
# Iron__xxx   -> Iron
#
# If "__" is absent, retain the original feature name.
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
# All merged subtypes belonging to the same type/class
# are summed within each sequencing sample.
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
  "\nAggregated MRG type matrix:\n"
)


cat(
  nrow(mrg_type_df),
  "MRG types ×",
  length(sample_cols),
  "samples\n"
)


# ======================================================
# 11. Export derived MRG type abundance table
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
# 12. Convert aggregated MRG type matrix to long format
# ======================================================

mrg_long <- mrg_type_df %>%

  pivot_longer(

    cols = -MRG_type,

    names_to = "ID",

    values_to = "Abundance"
  ) %>%

  mutate(

    ID = as.character(ID),

    Abundance = as.numeric(
      Abundance
    ),

    Abundance = ifelse(
      is.na(Abundance),
      0,
      Abundance
    ),

    Abundance = ifelse(
      Abundance < 0,
      0,
      Abundance
    )
  )


# ======================================================
# 13. Check sample matching
# ======================================================

matrix_samples <- unique(
  mrg_long$ID
)


common_samples <- intersect(
  meta_df$ID,
  matrix_samples
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
  matrix_samples
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
      "Too few matched samples. ",
      "Please check sample IDs in metadata.csv and ",
      "MRG_merged_subtype_abundance.tsv."
    )
  )
}


# ======================================================
# 14. Merge metadata
# ======================================================

plot_df_raw <- mrg_long %>%

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
    !is.na(Time),
    !is.na(Group)
  )


# ======================================================
# 15. Check sample numbers after merging
# ======================================================

cat(
  "\nSamples by Time × Group after merging:\n"
)


print(

  plot_df_raw %>%

    distinct(
      ID,
      Time,
      Group
    ) %>%

    count(
      Time,
      Group
    )
)


# ======================================================
# 16. Summarize MRG type abundance
#
# For each:
# Time × Group × MRG type
#
# Calculate:
#   median abundance
#   mean abundance
#   detection rate
# ======================================================

plot_df <- plot_df_raw %>%

  group_by(
    Time,
    Group,
    MRG_type
  ) %>%

  summarise(

    median_abundance = median(
      Abundance,
      na.rm = TRUE
    ),

    mean_abundance = mean(
      Abundance,
      na.rm = TRUE
    ),

    detection_rate = mean(
      Abundance > 0,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%

  mutate(

    log10_abundance = log10(
      median_abundance +
        1e-6
    )
  )


# ======================================================
# 17. Export summarized abundance data
# ======================================================

write.csv(

  plot_df,

  file.path(
    output_dir,
    "stat_MRG_type_median_abundance_TimeGroup.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 18. Order MRG types/classes
#
# Types with higher overall median abundance are shown
# toward the top of the y-axis.
# ======================================================

mrg_order <- plot_df %>%

  group_by(
    MRG_type
  ) %>%

  summarise(

    total_median = sum(
      median_abundance,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%

  arrange(
    total_median
  ) %>%

  pull(
    MRG_type
  )


plot_df <- plot_df %>%

  mutate(

    Time = factor(
      Time,
      levels = c(
        "D0",
        "D30",
        "D60",
        "D90"
      )
    ),

    Group = factor(
      Group,
      levels = c(
        "ITM",
        "OTM1",
        "OTM2"
      )
    ),

    MRG_type = factor(
      MRG_type,
      levels = mrg_order
    )
  )


# ======================================================
# 19. Bubble plot
# ======================================================

p_bubble_facet <- ggplot(

  plot_df,

  aes(
    x = Group,
    y = MRG_type
  )

) +


  geom_point(

    aes(
      size = median_abundance,
      color = log10_abundance
    ),

    alpha = 0.9
  ) +


  facet_grid(

    . ~ Time,

    scales = "free_x",

    space = "free_x"
  ) +


  scale_size_area(

    max_size = 10,

    name = "Median abundance"
  ) +


  scale_color_gradient(

    low = "#D9ECF2",

    high = "#2166AC",

    name = "log10(abundance)"
  ) +


  scale_x_discrete(
    drop = FALSE
  ) +


  labs(

    x = "Group",

    y = "MRG type",

    title =
      "MRG type composition across treatment groups and time points"
  ) +


  theme_bw(
    base_size = 13
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
        size = 16
      ),


    axis.title =
      element_text(
        face = "bold"
      ),


    axis.text.x =
      element_text(
        angle = 0,
        hjust = 0.5,
        color = "black",
        size = 10
      ),


    axis.text.y =
      element_text(
        color = "black",
        size = 10,
        face = "bold"
      ),


    legend.position =
      "right"
  )


print(
  p_bubble_facet
)


# ======================================================
# 20. Save figure
# ======================================================

figure_height <- max(

  5,

  length(
    unique(
      plot_df$MRG_type
    )
  ) *
    0.35
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_MRG_type_bubbleplot_facetTime.pdf"
  ),

  plot = p_bubble_facet,

  width = 9,

  height = figure_height
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_MRG_type_bubbleplot_facetTime.png"
  ),

  plot = p_bubble_facet,

  width = 9,

  height = figure_height,

  dpi = 600
)


# ======================================================
# 21. Finish
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "MRG type/class bubble plot completed.\n"
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