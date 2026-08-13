########################################################
## MGE-MRG overall structure association analysis
##
## Figure 3
##
## Procrustes analysis + Mantel test
##
## Input:
##   metadata/metadata.csv
##
##   processed_data/MGE/
##     MGE_merged_subtype_abundance.tsv
##
##   processed_data/MRG/
##     MRG_merged_subtype_abundance.tsv
##
## Analysis level:
##   MGE merged-subtype profile
##   MRG merged-subtype profile
##
## Analysis:
##   Post-baseline samples only:
##   D30, D60, and D90
##
##   1. Bray-Curtis dissimilarity for MGE profiles
##   2. Bray-Curtis dissimilarity for MRG profiles
##   3. MGE and MRG PCoA
##   4. Symmetric Procrustes analysis using PCoA1/2
##   5. Procrustes permutation test (999 permutations)
##   6. Spearman Mantel test (999 permutations)
##
## Samples with zero total abundance in either the MGE
## or MRG profile are excluded before distance analysis.
##
## Plot symbols:
##   MGE profile        = filled square
##   Fitted MRG profile = filled circle
##
## Output:
##   results/Fig3_MGE_MRG_association/
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


mge_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "MGE",
  "MGE_merged_subtype_abundance.tsv"
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
  "Fig3_MGE_MRG_association"
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


# Retain D30-D90 samples
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
  "MGE-MRG structure association analysis\n"
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
# 6. Read MGE and MRG abundance matrices
#
# Both files are expected to contain:
#
# First column:
#   merged subtype / feature
#
# Remaining columns:
#   sequencing sample IDs
# ======================================================

mge_df <- read.delim(
  mge_path,
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
# 8. Build MGE and MRG merged-subtype matrices
# ======================================================

mge_mat <- build_sample_matrix(
  feature_df = mge_df,
  matrix_name = "MGE merged-subtype matrix"
)


mrg_mat <- build_sample_matrix(
  feature_df = mrg_df,
  matrix_name = "MRG merged-subtype matrix"
)


cat(
  "\nMGE merged-subtype matrix:",
  nrow(mge_mat),
  "samples ×",
  ncol(mge_mat),
  "MGE merged subtypes\n"
)


cat(
  "MRG merged-subtype matrix:",
  nrow(mrg_mat),
  "samples ×",
  ncol(mrg_mat),
  "MRG merged subtypes\n"
)


# ======================================================
# 9. Match MGE, MRG, and metadata samples
# ======================================================

common_samples <- Reduce(

  intersect,

  list(

    rownames(
      mge_mat
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
      "Please check sample IDs in the MGE matrix, ",
      "MRG matrix, and metadata."
    )
  )
}


# Align MGE matrix
mge_mat <- mge_mat[
  common_samples,
  ,
  drop = FALSE
]


# Align MRG matrix
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
    rownames(mge_mat) ==
      rownames(mrg_mat)
  )
)


stopifnot(
  all(
    rownames(mge_mat) ==
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

mge_mat <- mge_mat[
  ,
  colSums(
    mge_mat,
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
  "\nFinal MGE merged subtypes:",
  ncol(mge_mat),
  "\n"
)


cat(
  "Final MRG merged subtypes:",
  ncol(mrg_mat),
  "\n"
)


# ======================================================
# 11. Remove samples with zero total abundance
#
# A sample is retained only when both MGE and MRG
# profiles contain at least one non-zero feature.
# ======================================================

sample_keep <-

  rowSums(
    mge_mat,
    na.rm = TRUE
  ) > 0 &

  rowSums(
    mrg_mat,
    na.rm = TRUE
  ) > 0


removed_zero_samples <- rownames(
  mge_mat
)[
  !sample_keep
]


if (length(removed_zero_samples) > 0) {

  cat(
    "\nSamples removed because the MGE or MRG profile ",
    "contained zero total abundance:\n",
    sep = ""
  )


  print(
    removed_zero_samples
  )
}


mge_mat <- mge_mat[
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
  "\nFinal number of samples used:",
  nrow(mge_mat),
  "\n"
)


cat(
  "\nFinal Group × Time distribution:\n"
)


print(
  table(
    meta_df$Group,
    meta_df$Time
  )
)


# ======================================================
# 12. Export excluded samples if present
# ======================================================

if (length(removed_zero_samples) > 0) {

  write.csv(

    tibble(

      ID =
        removed_zero_samples,

      Reason =
        "Zero total abundance in MGE and/or MRG profile"
    ),

    file.path(
      output_dir,
      "stat_D30D60D90_MGE_MRG_excluded_zero_samples.csv"
    ),

    row.names = FALSE
  )
}


# ======================================================
# 13. Bray-Curtis dissimilarity
# ======================================================

mge_bray <- vegan::vegdist(
  mge_mat,
  method = "bray"
)


mrg_bray <- vegan::vegdist(
  mrg_mat,
  method = "bray"
)


# ======================================================
# 14. Principal coordinate analysis
# ======================================================

mge_pcoa <- ape::pcoa(
  mge_bray
)


mrg_pcoa <- ape::pcoa(
  mrg_bray
)


# First two PCoA axes
mge_scores <- mge_pcoa$vectors[
  ,
  1:2,
  drop = FALSE
]


mrg_scores <- mrg_pcoa$vectors[
  ,
  1:2,
  drop = FALSE
]


rownames(mge_scores) <- rownames(
  mge_mat
)


rownames(mrg_scores) <- rownames(
  mrg_mat
)


# Variation explained
mge_eig <- round(

  100 *
    mge_pcoa$values$Relative_eig[
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
  "\nMGE PCoA1/2 explained variation:",
  mge_eig[1],
  "%,",
  mge_eig[2],
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
# MGE ordination is used as X.
# MRG ordination is rotated/fitted to MGE.
# ======================================================

set.seed(
  123
)


proc_res <- vegan::procrustes(

  X = mge_scores,

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

  X = mge_scores,

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
# MGE Bray-Curtis distance matrix
# and
# MRG Bray-Curtis distance matrix
#
# 999 permutations
# ======================================================

set.seed(
  123
)


mantel_res <- vegan::mantel(

  xdis = mge_bray,

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
# 18. Prepare Procrustes plotting coordinates
#
# proc_res$X:
# MGE ordination coordinates
#
# proc_res$Yrot:
# fitted/rotated MRG coordinates
# ======================================================

proc_plot_df <- data.frame(

  ID =
    rownames(
      proc_res$X
    ),

  MGE1 =
    proc_res$X[
      ,
      1
    ],

  MGE2 =
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
        MGE1 -
          MRG1
      )^2 +

        (
          MGE2 -
            MRG2
        )^2
    )
  )


# ======================================================
# 20. Prepare MGE and fitted MRG point coordinates
# ======================================================

mge_point_df <- proc_plot_df %>%

  transmute(

    ID,

    Group,

    Time,

    Axis1 =
      MGE1,

    Axis2 =
      MGE2,

    Profile =
      "MGE profile"
  )


mrg_point_df <- proc_plot_df %>%

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


point_df <- bind_rows(

  mge_point_df,

  mrg_point_df
)


point_df$Profile <- factor(

  point_df$Profile,

  levels = c(
    "MGE profile",
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
      mge_mat
    ),

  MGE_level =
    "Merged subtype",

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
    "stat_D30D60D90_MGE_MRG_Procrustes_Mantel.csv"
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
    "stat_D30D60D90_MGE_MRG_Protest_summary.txt"
  )
)


write.csv(

  proc_plot_df,

  file.path(
    output_dir,
    "stat_D30D60D90_MGE_MRG_Procrustes_plot_coordinates.csv"
  ),

  row.names = FALSE
)


write.csv(

  point_df,

  file.path(
    output_dir,
    "stat_D30D60D90_MGE_MRG_Procrustes_point_coordinates_long.csv"
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
# 24. P-value formatting
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
# 26. Determine annotation position
# ======================================================

x_range <- range(

  c(
    proc_plot_df$MGE1,
    proc_plot_df$MRG1
  ),

  na.rm = TRUE
)


y_range <- range(

  c(
    proc_plot_df$MGE2,
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
#
# MGE profile:
# filled square (shape 15)
#
# Fitted MRG profile:
# filled circle (shape 16)
# ======================================================

p_proc <- ggplot() +


  # Residual lines
  geom_segment(

    data =
      proc_plot_df,

    aes(

      x =
        MGE1,

      y =
        MGE2,

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


  # MGE profile
  geom_point(

    data =
      mge_point_df,

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
      2.9,

    alpha =
      0.90
  ) +


  # Fitted MRG profile
  geom_point(

    data =
      mrg_point_df,

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
      2.9,

    alpha =
      0.90
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

      "MGE profile" =
        15,

      "Fitted MRG profile" =
        16
    ),

    breaks = c(
      "MGE profile",
      "Fitted MRG profile"
    ),

    name =
      "Profile"
  ) +


  labs(

    x = paste0(
      "MGE PCoA1 (",
      mge_eig[1],
      "%)"
    ),

    y = paste0(
      "MGE PCoA2 (",
      mge_eig[2],
      "%)"
    ),

    title =
      "Procrustes analysis between MGE and MRG profiles",

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
    "Figure_D30D60D90_MGE_MRG_Procrustes_squareMGE.pdf"
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
    "Figure_D30D60D90_MGE_MRG_Procrustes_squareMGE.png"
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
  "D30-D90 MGE-MRG Procrustes + Mantel analysis completed.\n"
)

cat(
  "MGE input:\n"
)

cat(
  mge_path,
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