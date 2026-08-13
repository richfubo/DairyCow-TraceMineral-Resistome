########################################################
## ARG-MRG overall structure association analysis
##
## Figure 3
##
## Procrustes analysis + Mantel test
##
## Input:
##   metadata/metadata.csv
##
##   processed_data/ARG/
##     normalized_cell.subtype_matrix.tsv
##
##   processed_data/MRG/
##     MRG_merged_subtype_abundance.tsv
##
## Analysis level:
##   ARG subtype profile
##   MRG merged-subtype profile
##
## Analysis:
##   Post-baseline samples only:
##   D30, D60, and D90
##
##   1. Bray-Curtis dissimilarity for ARG profiles
##   2. Bray-Curtis dissimilarity for MRG profiles
##   3. ARG and MRG PCoA
##   4. Symmetric Procrustes analysis using PCoA1/2
##   5. Procrustes permutation test (999 permutations)
##   6. Spearman Mantel test (999 permutations)
##
## Purpose:
##   Evaluate whether the overall ARG and MRG resistome
##   structures covary during the post-treatment period.
##
## Output:
##   results/Fig3_ARG_MRG_association/
########################################################


# ======================================================
# 1. Load packages
# ======================================================

library(tidyverse)
library(vegan)
library(ape)
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


arg_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "ARG",
  "normalized_cell.subtype_matrix.tsv"
)


mrg_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "MRG",
  "MRG_merged_subtype_abundance.tsv"
)


output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig3_ARG_MRG_association"
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


# Retain D30-D90 only
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
  "ARG-MRG structure association analysis\n"
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
# 6. Read ARG and MRG abundance matrices
# ======================================================

arg_df <- read.delim(
  arg_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


mrg_df <- read.delim(
  mrg_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ======================================================
# 7. Function to construct sample × feature matrix
#
# Input:
#   feature × sample abundance matrix
#
# Output:
#   sample × feature abundance matrix
# ======================================================

build_sample_matrix <- function(
    feature_df,
    matrix_name
) {


  colnames(feature_df)[1] <- "feature"


  feature_df$feature <- as.character(
    feature_df$feature
  )


  # Check duplicated feature names
  if (anyDuplicated(feature_df$feature) > 0) {

    warning(
      paste0(
        "Duplicated feature names detected in ",
        matrix_name,
        ". make.unique() will be applied."
      )
    )
  }


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
    feature_df$feature
  )


  # Missing values -> zero
  mat[
    is.na(mat)
  ] <- 0


  # Negative abundance values are not allowed
  mat[
    mat < 0
  ] <- 0


  # Remove features absent from all samples
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
# 8. Build ARG subtype and MRG merged-subtype matrices
# ======================================================

arg_mat <- build_sample_matrix(
  feature_df = arg_df,
  matrix_name = "ARG subtype matrix"
)


mrg_mat <- build_sample_matrix(
  feature_df = mrg_df,
  matrix_name = "MRG merged-subtype matrix"
)


cat(
  "\nARG subtype matrix:",
  nrow(arg_mat),
  "samples ×",
  ncol(arg_mat),
  "ARG subtypes\n"
)


cat(
  "MRG merged-subtype matrix:",
  nrow(mrg_mat),
  "samples ×",
  ncol(mrg_mat),
  "MRG merged subtypes\n"
)


# ======================================================
# 9. Match ARG, MRG, and metadata samples
# ======================================================

common_samples <- Reduce(

  intersect,

  list(

    rownames(
      arg_mat
    ),

    rownames(
      mrg_mat
    ),

    meta_df$ID
  )
)


cat(
  "\nMatched D30-D90 samples:",
  length(common_samples),
  "\n"
)


# Expected:
# 30 cows × 3 post-baseline time points = 90 samples
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
      "Please check sample IDs in ARG, MRG, ",
      "and metadata."
    )
  )
}


# Align ARG
arg_mat <- arg_mat[
  common_samples,
  ,
  drop = FALSE
]


# Align MRG
mrg_mat <- mrg_mat[
  common_samples,
  ,
  drop = FALSE
]


# Align metadata
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
    rownames(arg_mat) ==
      rownames(mrg_mat)
  )
)


stopifnot(
  all(
    rownames(arg_mat) ==
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
# 10. Remove features absent from matched samples
# ======================================================

arg_mat <- arg_mat[
  ,
  colSums(
    arg_mat,
    na.rm = TRUE
  ) > 0,
  drop = FALSE
]


mrg_mat <- mrg_mat[
  ,
  colSums(
    mrg_mat,
    na.rm = TRUE
  ) > 0,
  drop = FALSE
]


cat(
  "\nFinal ARG subtypes:",
  ncol(arg_mat),
  "\n"
)


cat(
  "Final MRG merged subtypes:",
  ncol(mrg_mat),
  "\n"
)


# ======================================================
# 11. Check zero-abundance samples
#
# Bray-Curtis is not informative when a sample has
# zero abundance across all features.
# ======================================================

sample_keep <-

  rowSums(
    arg_mat,
    na.rm = TRUE
  ) > 0 &

  rowSums(
    mrg_mat,
    na.rm = TRUE
  ) > 0


removed_samples <- rownames(
  arg_mat
)[
  !sample_keep
]


if (length(removed_samples) > 0) {

  cat(
    "\nSamples removed because ARG or MRG profile ",
    "contained zero total abundance:\n",
    sep = ""
  )


  print(
    removed_samples
  )
}


arg_mat <- arg_mat[
  sample_keep,
  ,
  drop = FALSE
]


mrg_mat <- mrg_mat[
  sample_keep,
  ,
  drop = FALSE
]


meta_df <- meta_df[
  sample_keep,
  ,
  drop = FALSE
]


cat(
  "\nFinal number of samples:",
  nrow(arg_mat),
  "\n"
)


# ======================================================
# 12. Export excluded samples if applicable
# ======================================================

if (length(removed_samples) > 0) {

  write.csv(

    tibble(
      ID = removed_samples,
      Reason = "Zero total abundance in ARG and/or MRG profile"
    ),

    file.path(
      output_dir,
      "stat_D30D60D90_ARG_MRG_excluded_zero_samples.csv"
    ),

    row.names = FALSE
  )
}


# ======================================================
# 13. Bray-Curtis dissimilarity
# ======================================================

arg_bray <- vegan::vegdist(
  arg_mat,
  method = "bray"
)


mrg_bray <- vegan::vegdist(
  mrg_mat,
  method = "bray"
)


# ======================================================
# 14. Principal coordinate analysis
# ======================================================

arg_pcoa <- ape::pcoa(
  arg_bray
)


mrg_pcoa <- ape::pcoa(
  mrg_bray
)


arg_scores <- arg_pcoa$vectors[
  ,
  1:2,
  drop = FALSE
]


mrg_scores <- mrg_pcoa$vectors[
  ,
  1:2,
  drop = FALSE
]


rownames(arg_scores) <- rownames(
  arg_mat
)


rownames(mrg_scores) <- rownames(
  mrg_mat
)


# Variation explained by first two axes
arg_eig <- round(

  100 *
    arg_pcoa$values$Relative_eig[
      1:2
    ],

  1
)


mrg_eig <- round(

  100 *
    mrg_pcoa$values$Relative_eig[
      1:2
    ],

  1
)


cat(
  "\nARG PCoA1/2 explained variation:",
  arg_eig[1],
  "%,",
  arg_eig[2],
  "%\n"
)


cat(
  "MRG PCoA1/2 explained variation:",
  mrg_eig[1],
  "%,",
  mrg_eig[2],
  "%\n"
)


# ======================================================
# 15. Symmetric Procrustes analysis
#
# ARG ordination = X
# MRG ordination = Y
# ======================================================

set.seed(
  123
)


proc_res <- vegan::procrustes(

  X = arg_scores,

  Y = mrg_scores,

  symmetric = TRUE
)


print(
  proc_res
)


# ======================================================
# 16. Procrustes permutation test
#
# 999 permutations
# ======================================================

set.seed(
  123
)


protest_res <- vegan::protest(

  X = arg_scores,

  Y = mrg_scores,

  permutations = 999,

  symmetric = TRUE
)


print(
  protest_res
)


proc_M2 <- proc_res$ss


proc_r <- protest_res$t0


proc_p <- protest_res$signif


cat(
  "\nProcrustes results:\n"
)


cat(
  "M2 =",
  proc_M2,
  "\n"
)


cat(
  "Correlation r =",
  proc_r,
  "\n"
)


cat(
  "P =",
  proc_p,
  "\n"
)


# ======================================================
# 17. Mantel test
#
# Spearman correlation between:
#
# ARG Bray-Curtis distance matrix
# and
# MRG Bray-Curtis distance matrix
#
# 999 permutations
# ======================================================

set.seed(
  123
)


mantel_res <- vegan::mantel(

  xdis = arg_bray,

  ydis = mrg_bray,

  method = "spearman",

  permutations = 999
)


print(
  mantel_res
)


mantel_r <- as.numeric(
  mantel_res$statistic
)


mantel_p <- mantel_res$signif


cat(
  "\nMantel results:\n"
)


cat(
  "Mantel r =",
  mantel_r,
  "\n"
)


cat(
  "P =",
  mantel_p,
  "\n"
)


# ======================================================
# 18. Prepare Procrustes coordinates
# ======================================================

proc_plot_df <- data.frame(

  ID =
    rownames(
      proc_res$X
    ),

  ARG1 =
    proc_res$X[
      ,
      1
    ],

  ARG2 =
    proc_res$X[
      ,
      2
    ],

  MRG1 =
    proc_res$Yrot[
      ,
      1
    ],

  MRG2 =
    proc_res$Yrot[
      ,
      2
    ]

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
# 19. Calculate Procrustes residual lengths
# ======================================================

proc_plot_df <- proc_plot_df %>%

  mutate(

    residual_length = sqrt(

      (
        ARG1 -
          MRG1
      )^2 +

        (
          ARG2 -
            MRG2
        )^2
    )
  )


# ======================================================
# 20. Convert coordinates to long format
# ======================================================

point_df <- bind_rows(


  proc_plot_df %>%

    transmute(

      ID,

      Group,

      Time,

      Axis1 =
        ARG1,

      Axis2 =
        ARG2,

      Profile =
        "ARG profile"
    ),


  proc_plot_df %>%

    transmute(

      ID,

      Group,

      Time,

      Axis1 =
        MRG1,

      Axis2 =
        MRG2,

      Profile =
        "Fitted MRG profile"
    )
)


point_df$Profile <- factor(

  point_df$Profile,

  levels = c(
    "ARG profile",
    "Fitted MRG profile"
  )
)


# ======================================================
# 21. Compile statistical summary
# ======================================================

stat_df <- tibble(

  Analysis = c(
    "Procrustes",
    "Protest",
    "Mantel"
  ),

  Statistic = c(
    "M2",
    "r",
    "Mantel r"
  ),

  Value = c(
    proc_M2,
    proc_r,
    mantel_r
  ),

  P_value = c(
    NA,
    proc_p,
    mantel_p
  ),

  Method = c(

    "Symmetric Procrustes sum of squares",

    "Symmetric Procrustes permutation test (999 permutations)",

    "Spearman Mantel test (999 permutations)"
  ),

  Samples_used =
    paste(
      use_times,
      collapse = " + "
    ),

  N_samples =
    nrow(
      arg_mat
    ),

  ARG_level =
    "Subtype",

  MRG_level =
    "Merged subtype"
)


# ======================================================
# 22. Export statistical results
# ======================================================

write.csv(

  stat_df,

  file.path(
    output_dir,
    "stat_D30D60D90_ARG_MRG_Procrustes_Mantel.csv"
  ),

  row.names = FALSE
)


capture.output(

  list(

    Procrustes =
      proc_res,

    Protest =
      protest_res,

    Mantel =
      mantel_res
  ),

  file = file.path(
    output_dir,
    "stat_D30D60D90_ARG_MRG_Protest_summary.txt"
  )
)


write.csv(

  proc_plot_df,

  file.path(
    output_dir,
    "stat_D30D60D90_ARG_MRG_Procrustes_plot_coordinates.csv"
  ),

  row.names = FALSE
)


write.csv(

  point_df,

  file.path(
    output_dir,
    "stat_D30D60D90_ARG_MRG_Procrustes_point_coordinates_long.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 23. Plot colors
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
# 24. Format P values
# ======================================================

format_p <- function(
    p
) {


  if (is.na(p)) {

    return(
      "NA"
    )
  }


  if (p < 0.001) {

    return(
      "< 0.001"
    )

  } else {

    return(
      sprintf(
        "%.3f",
        p
      )
    )
  }
}


# ======================================================
# 25. Statistical annotation
# ======================================================

proc_label <- paste0(

  "Procrustes:\n",

  "M\u00b2 = ",
  sprintf(
    "%.3f",
    proc_M2
  ),

  "\n",

  "r = ",
  sprintf(
    "%.3f",
    proc_r
  ),

  "\n",

  "P = ",
  format_p(
    proc_p
  ),

  "\n\n",

  "Mantel:\n",

  "r = ",
  sprintf(
    "%.3f",
    mantel_r
  ),

  "\n",

  "P = ",
  format_p(
    mantel_p
  )
)


# ======================================================
# 26. Annotation position
# ======================================================

x_range <- range(

  c(
    proc_plot_df$ARG1,
    proc_plot_df$MRG1
  ),

  na.rm = TRUE
)


y_range <- range(

  c(
    proc_plot_df$ARG2,
    proc_plot_df$MRG2
  ),

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
# 27. Procrustes plot
# ======================================================

p_proc <- ggplot() +


  # Residual connections
  geom_segment(

    data =
      proc_plot_df,

    aes(

      x =
        ARG1,

      y =
        ARG2,

      xend =
        MRG1,

      yend =
        MRG2,

      color =
        Group
    ),

    linewidth =
      0.45,

    alpha =
      0.55
  ) +


  # ARG and fitted MRG points
  geom_point(

    data =
      point_df,

    aes(

      x =
        Axis1,

      y =
        Axis2,

      color =
        Group,

      shape =
        Profile
    ),

    size =
      2.8,

    alpha =
      0.88
  ) +


  annotate(

    "text",

    x =
      label_x,

    y =
      label_y,

    label =
      proc_label,

    hjust =
      0,

    vjust =
      0,

    size =
      4,

    fontface =
      "italic"
  ) +


  scale_color_manual(

    values =
      group_cols,

    name =
      "Group"
  ) +


  scale_shape_manual(

    values = c(

      "ARG profile" =
        16,

      "Fitted MRG profile" =
        17
    ),

    name =
      "Profile"
  ) +


  labs(

    x = paste0(
      "ARG PCoA1 (",
      arg_eig[1],
      "%)"
    ),

    y = paste0(
      "ARG PCoA2 (",
      arg_eig[2],
      "%)"
    ),

    title =
      "Procrustes analysis between ARG and MRG profiles",

    subtitle =
      "Post-baseline samples only (D30-D90)"
  ) +


  guides(

    color =
      guide_legend(
        order = 1
      ),

    shape =
      guide_legend(
        order = 2
      )
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
        size = 16
      ),


    plot.subtitle =
      element_text(
        hjust = 0.5,
        size = 11
      ),


    axis.title =
      element_text(
        face = "bold"
      ),


    axis.text =
      element_text(
        color = "black"
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
  p_proc
)


# ======================================================
# 28. Save figure
# ======================================================

ggsave(

  filename = file.path(
    output_dir,
    "Figure_D30D60D90_ARG_MRG_Procrustes.pdf"
  ),

  plot =
    p_proc,

  width =
    8,

  height =
    6
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_D30D60D90_ARG_MRG_Procrustes.png"
  ),

  plot =
    p_proc,

  width =
    8,

  height =
    6,

  dpi =
    600
)


# ======================================================
# 29. Finish
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "D30-D90 ARG-MRG Procrustes + Mantel analysis completed.\n"
)

cat(
  "ARG input:\n"
)

cat(
  arg_path,
  "\n"
)

cat(
  "MRG input:\n"
)

cat(
  mrg_path,
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