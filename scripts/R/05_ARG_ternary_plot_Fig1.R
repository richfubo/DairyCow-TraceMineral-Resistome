########################################################
## ARG ternary plot
##
## Figure 1
##
## Input:
##   metadata/metadata.csv
##   processed_data/ARG/normalized_cell.subtype_matrix.tsv
##
## Analysis:
##   1. log10(x + 1) transformation of ARG abundance
##   2. Bray-Curtis dissimilarity
##   3. PCoA
##   4. Group centroid distances at each time point
##   5. Ternary distribution of each ARG subtype among
##      ITM, OTM1, and OTM2
##
## Each point in the ternary plot represents one ARG
## subtype at one sampling time.
##
## Point size:
##   Total mean normalized abundance across three groups
##
## Point color:
##   Relative contribution of the ITM group (%)
########################################################


# ======================================================
# 1. Load packages
# ======================================================

library(tidyverse)
library(ggtern)
library(vegan)


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
  "Fig1_ARG_ternary"
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


if (anyDuplicated(meta_df$ID) > 0) {

  stop(
    "Duplicated sequencing sample IDs were found."
  )
}


# ======================================================
# 6. Read normalized ARG subtype abundance matrix
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


# Convert abundance columns to numeric
feature_df[
  ,
  -1
] <- lapply(

  feature_df[
    ,
    -1,
    drop = FALSE
  ],

  function(x) {
    as.numeric(
      as.character(x)
    )
  }
)


# Replace missing and negative values
feature_df[
  ,
  -1
][
  is.na(
    feature_df[
      ,
      -1
    ]
  )
] <- 0


for (i in 2:ncol(feature_df)) {

  feature_df[[i]][
    feature_df[[i]] < 0
  ] <- 0
}


# ======================================================
# 7. Construct sample × ARG subtype matrix
# ======================================================

x <- feature_df[
  ,
  -1,
  drop = FALSE
]


mat_raw <- t(
  as.matrix(x)
)


rownames(mat_raw) <- colnames(
  feature_df
)[-1]


colnames(mat_raw) <- make.unique(
  feature_df$feature
)


# ======================================================
# 8. Match samples
# ======================================================

common_samples <- intersect(
  meta_df$ID,
  rownames(mat_raw)
)


cat(
  "\nSamples matched between metadata and ARG table:",
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
  rownames(mat_raw)
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
    "Too few matched samples."
  )
}


mat_raw <- mat_raw[
  common_samples,
  ,
  drop = FALSE
]


meta_df_pcoa <- meta_df[
  match(
    common_samples,
    meta_df$ID
  ),
  ,
  drop = FALSE
]


stopifnot(
  all(
    rownames(mat_raw) ==
      meta_df_pcoa$ID
  )
)


# Remove ARG subtypes absent from all samples
mat_raw <- mat_raw[
  ,
  colSums(
    mat_raw,
    na.rm = TRUE
  ) > 0,
  drop = FALSE
]


cat(
  "ARG subtypes used:",
  ncol(mat_raw),
  "\n"
)


# ======================================================
# 9. Calculate PCoA for centroid-distance estimation
#
# Original analysis:
# log10(x + 1) -> Bray-Curtis -> PCoA
# ======================================================

mat_log <- log10(
  mat_raw + 1
)


dist_matrix <- vegan::vegdist(
  mat_log,
  method = "bray"
)


pcoa <- cmdscale(
  dist_matrix,
  k = 2
)


plot_data <- as.data.frame(
  pcoa
)


colnames(
  plot_data
) <- c(
  "PCoA1",
  "PCoA2"
)


plot_data$ID <- rownames(
  plot_data
)


plot_data <- plot_data %>%

  left_join(
    meta_df_pcoa,
    by = "ID"
  )


# ======================================================
# 10. Calculate group centroids and centroid distances
#     at each sampling time
# ======================================================

all_times <- levels(
  meta_df$Time
)


facet_labels <- tibble()


centroid_distance_df <- tibble()


for (t in all_times) {


  pts <- plot_data %>%

    filter(
      Time == t
    ) %>%

    group_by(
      Group
    ) %>%

    summarise(

      X = mean(
        PCoA1,
        na.rm = TRUE
      ),

      Y = mean(
        PCoA2,
        na.rm = TRUE
      ),

      .groups = "drop"
    )


  required_groups <- c(
    "ITM",
    "OTM1",
    "OTM2"
  )


  if (
    all(
      required_groups %in%
        pts$Group
    )
  ) {


    c_i <- as.numeric(

      pts[
        pts$Group == "ITM",
        c(
          "X",
          "Y"
        )
      ]
    )


    c_o1 <- as.numeric(

      pts[
        pts$Group == "OTM1",
        c(
          "X",
          "Y"
        )
      ]
    )


    c_o2 <- as.numeric(

      pts[
        pts$Group == "OTM2",
        c(
          "X",
          "Y"
        )
      ]
    )


    d_i_o1 <- sqrt(
      sum(
        (
          c_i -
            c_o1
        )^2
      )
    )


    d_i_o2 <- sqrt(
      sum(
        (
          c_i -
            c_o2
        )^2
      )
    )


    d_o1_o2 <- sqrt(
      sum(
        (
          c_o1 -
            c_o2
        )^2
      )
    )


    # Save full-precision values
    centroid_distance_df <- bind_rows(

      centroid_distance_df,

      tibble(

        Time = t,

        ITM_OTM1_distance =
          d_i_o1,

        ITM_OTM2_distance =
          d_i_o2,

        OTM1_OTM2_distance =
          d_o1_o2
      )
    )


    # Rounded values for facet labels
    new_label <- paste0(

      t,

      "\n",

      "I-O1=",
      round(
        d_i_o1,
        3
      ),

      " | I-O2=",
      round(
        d_i_o2,
        3
      ),

      " | O1-O2=",
      round(
        d_o1_o2,
        3
      )
    )


    facet_labels <- bind_rows(

      facet_labels,

      tibble(

        Time = t,

        FullLabel =
          new_label
      )
    )
  }
}


# Set time order
facet_labels$Time <- factor(
  facet_labels$Time,
  levels = all_times
)


centroid_distance_df$Time <- factor(
  centroid_distance_df$Time,
  levels = all_times
)


# ======================================================
# 11. Prepare ternary plot data
#
# Each point:
# one ARG subtype at one sampling time
#
# Mean normalized abundance is calculated separately
# for ITM, OTM1, and OTM2.
# ======================================================

tern_data <- feature_df %>%

  pivot_longer(

    cols = -feature,

    names_to = "ID",

    values_to = "Abundance"
  ) %>%

  mutate(

    ID = as.character(
      ID
    ),

    Abundance = as.numeric(
      Abundance
    )
  ) %>%

  inner_join(

    meta_df %>%
      select(
        ID,
        Time,
        Group
      ),

    by = "ID"
  ) %>%

  group_by(
    feature,
    Time,
    Group
  ) %>%

  summarise(

    MeanAbundance = mean(
      Abundance,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%

  pivot_wider(

    names_from = Group,

    values_from = MeanAbundance,

    values_fill = 0
  )


# ======================================================
# 12. Convert abundance to ternary proportions
# ======================================================

tern_final_data <- tern_data %>%

  mutate(

    Total =
      ITM +
      OTM1 +
      OTM2
  ) %>%

  filter(
    Total > 0
  ) %>%

  mutate(

    p_ITM =
      (
        ITM /
          Total
      ) *
      100,

    p_OTM1 =
      (
        OTM1 /
          Total
      ) *
      100,

    p_OTM2 =
      (
        OTM2 /
          Total
      ) *
      100
  ) %>%

  left_join(
    facet_labels,
    by = "Time"
  )


# ======================================================
# 13. Plot ternary distribution
# ======================================================

p_final <- ggtern(

  data = tern_final_data,

  aes(
    x = p_OTM1,
    y = p_ITM,
    z = p_OTM2
  )

) +


  geom_point(

    aes(
      size = Total,
      color = p_ITM
    ),

    alpha = 0.5
  ) +


  scale_color_gradientn(

    colors = c(
      "#FA7F6F",
      "#FFBE7A",
      "#8ECFC9"
    ),

    name = "ITM preference (%)"
  ) +


  scale_size(

    range = c(
      1,
      6
    ),

    name = "Total abundance"
  ) +


  facet_wrap(

    ~ FullLabel,

    ncol = 2
  ) +


  theme_rgbw() +

  theme_showarrows() +


  Tlab(
    "ITM"
  ) +

  Llab(
    "OTM1"
  ) +

  Rlab(
    "OTM2"
  ) +


  theme(

    tern.axis.title.T =
      element_text(
        color = "#8ECFC9",
        face = "bold",
        size = 12
      ),

    tern.axis.title.L =
      element_text(
        color = "#FFBE7A",
        face = "bold",
        size = 12
      ),

    tern.axis.title.R =
      element_text(
        color = "#FA7F6F",
        face = "bold",
        size = 12
      ),

    strip.background =
      element_rect(
        fill = "white",
        color = NA
      ),

    strip.text =
      element_text(
        size = 11,
        face = "bold",
        color = "black"
      ),

    legend.position =
      "right",

    legend.box =
      "vertical",

    plot.margin =
      ggplot2::margin(
        1,
        1,
        1,
        1,
        "cm"
      )
  ) +


  labs(
    title = NULL
  )


print(
  p_final
)


# ======================================================
# 14. Save figure
# ======================================================

ggsave(

  filename = file.path(
    output_dir,
    "Figure_ARG_Ternary_FacetTitleDist.pdf"
  ),

  plot = p_final,

  width = 11,

  height = 8
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_ARG_Ternary_FacetTitleDist.png"
  ),

  plot = p_final,

  width = 11,

  height = 8,

  dpi = 600
)


# ======================================================
# 15. Export data used for ternary plot
# ======================================================

write.csv(

  tern_final_data,

  file.path(
    output_dir,
    "stat_ARG_ternary_data.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 16. Export centroid distances
# ======================================================

write.csv(

  centroid_distance_df,

  file.path(
    output_dir,
    "stat_ARG_PCoA_centroid_distances.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 17. Finish
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "ARG ternary plot analysis completed.\n"
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