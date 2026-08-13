########################################################
## 16S beta diversity analysis
##
## Figure 4
##
## Bray-Curtis PCoA + pairwise PERMANOVA
##
## Input:
##   metadata/metadata.csv
##
##   processed_data/16S/
##     BrayCurtis_distance_matrix.tsv
##
## Analysis:
##   PCoA is performed separately at:
##   D0, D30, D60, and D90
##
## Pairwise PERMANOVA within each time point:
##   ITM vs OTM1
##   ITM vs OTM2
##   OTM1 vs OTM2
##
## PERMANOVA:
##   999 permutations
##
## Output:
##   results/Fig4_16S_PCoA/
########################################################


# ======================================================
# 1. Load packages
# ======================================================

library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)


# ======================================================
# 2. Project paths
#
# Run this script from the repository root:
# DairyCow-TraceMineral-Resistome/
# ======================================================

PROJECT_DIR <- "."


distance_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "16S",
  "BrayCurtis_distance_matrix.tsv"
)


meta_path <- file.path(
  PROJECT_DIR,
  "metadata",
  "metadata.csv"
)


output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig4_16S_PCoA"
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
      "Missing metadata column(s): ",
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

    Group = as.character(Group)
  )


# ======================================================
# 4. Check metadata
# ======================================================

expected_groups <- c(
  "ITM",
  "OTM1",
  "OTM2"
)


expected_times <- c(
  "D0",
  "D30",
  "D60",
  "D90"
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


meta_df$Group <- factor(
  meta_df$Group,
  levels = expected_groups
)


meta_df$Time <- factor(
  meta_df$Time,
  levels = expected_times
)


cat(
  "\n========================================\n"
)

cat(
  "16S beta-diversity analysis\n"
)

cat(
  "========================================\n"
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
# 5. Read Bray-Curtis distance matrix
#
# Expected format:
#
#       a1     a2     a3 ...
# a1    0      ...
# a2    ...    0
# a3    ...           0
#
# Row names and column names must be sample IDs.
# ======================================================

dist_matrix <- read.table(

  distance_path,

  header = TRUE,

  row.names = 1,

  sep = "\t",

  check.names = FALSE,

  stringsAsFactors = FALSE
)


dist_matrix <- as.matrix(
  dist_matrix
)


storage.mode(
  dist_matrix
) <- "numeric"


# ======================================================
# 6. Check distance matrix
# ======================================================

if (
  nrow(dist_matrix) !=
    ncol(dist_matrix)
) {

  stop(
    "The Bray-Curtis distance matrix is not square."
  )
}


if (
  !identical(
    rownames(dist_matrix),
    colnames(dist_matrix)
  )
) {

  stop(
    paste0(
      "Row names and column names of the distance ",
      "matrix are not identical."
    )
  )
}


if (
  any(
    is.na(
      dist_matrix
    )
  )
) {

  stop(
    "NA values were detected in the distance matrix."
  )
}


if (
  any(
    dist_matrix < 0
  )
) {

  stop(
    "Negative values were detected in the distance matrix."
  )
}


# ======================================================
# 7. Match metadata and distance matrix samples
# ======================================================

common_samples <- intersect(
  meta_df$ID,
  rownames(
    dist_matrix
  )
)


cat(
  "\nMetadata samples:",
  nrow(meta_df),
  "\n"
)


cat(
  "Distance-matrix samples:",
  nrow(dist_matrix),
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


missing_in_distance <- setdiff(
  meta_df$ID,
  rownames(
    dist_matrix
  )
)


if (length(missing_in_distance) > 0) {

  cat(
    "\nSamples in metadata but missing from distance matrix:\n"
  )

  print(
    missing_in_distance
  )
}


missing_in_metadata <- setdiff(
  rownames(
    dist_matrix
  ),
  meta_df$ID
)


if (length(missing_in_metadata) > 0) {

  cat(
    "\nSamples in distance matrix but missing from metadata:\n"
  )

  print(
    missing_in_metadata
  )
}


if (length(common_samples) < 5) {

  stop(
    paste0(
      "Too few matched samples. ",
      "Please check sample IDs."
    )
  )
}


# ======================================================
# 8. Restrict to matched samples
# ======================================================

dist_matrix <- dist_matrix[
  common_samples,
  common_samples,
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
      rownames(
        dist_matrix
      )
  )
)


# ======================================================
# 9. Analysis settings
# ======================================================

times <- c(
  "D0",
  "D30",
  "D60",
  "D90"
)


group_levels <- c(
  "ITM",
  "OTM1",
  "OTM2"
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


group_cols <- c(

  ITM =
    "#8ECFC9",

  OTM1 =
    "#FFBE7A",

  OTM2 =
    "#FA7F6F"
)


# ======================================================
# 10. Initialize result tables
# ======================================================

all_results <- tibble()


all_coordinates <- tibble()


# ======================================================
# 11. P-value formatting
# ======================================================

format_p <- function(p) {

  if (is.na(p)) {

    return(
      "NA"
    )
  }


  if (p < 0.001) {

    return(
      "< 0.001"
    )
  }


  sprintf(
    "%.3f",
    p
  )
}


# ======================================================
# 12. Main loop: analyze each time point separately
# ======================================================

for (time in times) {


  cat(
    "\n========================================\n"
  )


  cat(
    "Processing:",
    time,
    "\n"
  )


  cat(
    "========================================\n"
  )


  # ----------------------------------------------------
  # 12.1 Metadata for current time point
  # ----------------------------------------------------

  meta_sub <- meta_df %>%

    filter(
      Time == time
    ) %>%

    droplevels()


  meta_sub$Group <- factor(
    meta_sub$Group,
    levels = group_levels
  )


  samples <- meta_sub$ID


  cat(
    "Samples:",
    length(samples),
    "\n"
  )


  print(
    table(
      meta_sub$Group
    )
  )


  if (length(samples) < 5) {

    warning(
      paste0(
        "Too few samples at ",
        time,
        ". This time point was skipped."
      )
    )

    next
  }


  # ----------------------------------------------------
  # 12.2 Subset distance matrix
  # ----------------------------------------------------

  dist_sub_matrix <- dist_matrix[
    samples,
    samples,
    drop = FALSE
  ]


  dist_sub <- as.dist(
    dist_sub_matrix
  )


  # ----------------------------------------------------
  # 12.3 PCoA
  #
  # Original analysis uses cmdscale().
  # ----------------------------------------------------

  pcoa <- cmdscale(

    dist_sub,

    eig = TRUE,

    k = 2
  )


  points <- as.data.frame(
    pcoa$points
  )


  colnames(points) <- c(
    "PC1",
    "PC2"
  )


  points <- points %>%

    rownames_to_column(
      "ID"
    ) %>%

    left_join(

      meta_sub %>%

        select(
          ID,
          Group,
          Time
        ),

      by = "ID"
    )


  # ----------------------------------------------------
  # 12.4 Variation explained
  #
  # Retain the calculation used in the original script.
  # ----------------------------------------------------

  eig <- pcoa$eig


  pc1_var <- round(

    eig[1] /
      sum(eig) *
      100,

    2
  )


  pc2_var <- round(

    eig[2] /
      sum(eig) *
      100,

    2
  )


  cat(
    "PC1 explained variation:",
    pc1_var,
    "%\n"
  )


  cat(
    "PC2 explained variation:",
    pc2_var,
    "%\n"
  )


  # Save coordinates
  all_coordinates <- bind_rows(

    all_coordinates,

    points %>%

      mutate(

        PC1_variance_percent =
          pc1_var,

        PC2_variance_percent =
          pc2_var
      )
  )


  # ----------------------------------------------------
  # 12.5 Pairwise PERMANOVA
  # ----------------------------------------------------

  label_texts <- c()


  for (
    pair in group_comparisons
  ) {


    pair_meta <- meta_sub %>%

      filter(
        Group %in% pair
      ) %>%

      droplevels()


    pair_samples <- pair_meta$ID


    # Require samples from both groups
    if (
      length(
        unique(
          pair_meta$Group
        )
      ) < 2
    ) {

      next
    }


    pair_dist_matrix <- dist_sub_matrix[
      pair_samples,
      pair_samples,
      drop = FALSE
    ]


    pair_dist <- as.dist(
      pair_dist_matrix
    )


    set.seed(
      123
    )


    adonis_res <- vegan::adonis2(

      pair_dist ~ Group,

      data = pair_meta,

      permutations = 999
    )


    pval <- adonis_res$`Pr(>F)`[
      1
    ]


    fval <- adonis_res$F[
      1
    ]


    r2val <- adonis_res$R2[
      1
    ]


    all_results <- bind_rows(

      all_results,

      tibble(

        Time =
          time,

        Group1 =
          pair[1],

        Group2 =
          pair[2],

        Comparison =
          paste0(
            pair[1],
            " vs ",
            pair[2]
          ),

        F =
          fval,

        R2 =
          r2val,

        P_value =
          pval,

        Permutations =
          999,

        N_group1 =
          sum(
            pair_meta$Group ==
              pair[1]
          ),

        N_group2 =
          sum(
            pair_meta$Group ==
              pair[2]
          )
      )
    )


    label_texts <- c(

      label_texts,

      paste0(
        pair[1],
        " vs ",
        pair[2],
        ": P ",
        ifelse(
          pval < 0.001,
          "< 0.001",
          paste0(
            "= ",
            format_p(
              pval
            )
          )
        )
      )
    )
  }


  # ----------------------------------------------------
  # 12.6 Statistical annotation
  # ----------------------------------------------------

  label_line <- paste(
    label_texts,
    collapse = " | "
  )


  # ----------------------------------------------------
  # 12.7 PCoA plot
  # ----------------------------------------------------

  p <- ggplot(

    points,

    aes(
      x = PC1,
      y = PC2,