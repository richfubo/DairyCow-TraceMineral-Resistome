########################################################
## Differential ARG type abundance during D30-D90
##
## Figure 1
##
## Input:
##   metadata/metadata.csv
##   processed_data/ARG/normalized_cell.type_matrix.tsv
##
## Analysis:
##   D30, D60, and D90 samples are pooled.
##
## Pairwise comparisons:
##   OTM1 vs ITM
##   OTM2 vs ITM
##   OTM2 vs OTM1
##
## Statistical test:
##   Wilcoxon rank-sum test (unpaired)
##
## Multiple-testing correction:
##   Benjamini-Hochberg (BH/FDR)
##   applied separately within each pairwise comparison
##
## Volcano plot:
##   x-axis = log2 fold change
##   y-axis = -log10(FDR)
##   labels = ARG types with FDR < 0.05
########################################################


# ======================================================
# 1. Load packages
# ======================================================

library(tidyverse)
library(ggrepel)


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
  "normalized_cell.type_matrix.tsv"
)


output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig1_ARG_type_volcano"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Sampling times included in this analysis
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


# Set factor levels
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


# Retain D30, D60, and D90
meta_df <- meta_df %>%

  filter(
    Time %in% use_times
  ) %>%

  droplevels()


# ======================================================
# 5. Metadata checks
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "ARG type differential analysis\n"
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
  "\nTotal samples per group:\n"
)


print(
  table(
    meta_df$Group
  )
)


if (anyDuplicated(meta_df$ID) > 0) {

  stop(
    "Duplicated sequencing sample IDs were found."
  )
}


# ======================================================
# 6. Read normalized ARG type abundance matrix
#
# First column:
# ARG type
#
# Remaining columns:
# sample IDs (a1-a120)
# ======================================================

feature_df <- read.delim(
  feature_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


colnames(feature_df)[1] <- "ARG_type"


feature_df$ARG_type <- as.character(
  feature_df$ARG_type
)


# ======================================================
# 7. Construct sample × ARG type matrix
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
    feature_df$ARG_type
  )
)


# Replace missing values with zero
mat[is.na(mat)] <- 0


# Negative abundance values are not allowed
mat[mat < 0] <- 0


# Remove ARG types absent from all samples
mat <- mat[
  ,
  colSums(
    mat,
    na.rm = TRUE
  ) > 0,
  drop = FALSE
]


cat(
  "\nFull ARG matrix:",
  nrow(mat),
  "samples ×",
  ncol(mat),
  "ARG types\n"
)


# ======================================================
# 8. Match samples
# ======================================================

common_samples <- intersect(
  meta_df$ID,
  rownames(mat)
)


cat(
  "\nMatched D30-D90 samples:",
  length(common_samples),
  "\n"
)


if (length(common_samples) < 5) {

  stop(
    paste0(
      "Too few matched samples. ",
      "Please check sample IDs in metadata.csv ",
      "and normalized_cell.type_matrix.tsv."
    )
  )
}


missing_samples <- setdiff(
  meta_df$ID,
  rownames(mat)
)


if (length(missing_samples) > 0) {

  cat(
    "\nSamples missing from ARG type matrix:\n"
  )

  print(
    missing_samples
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
    rownames(mat) ==
      meta_df$ID
  )
)


cat(
  "\nSamples used after matching:\n"
)


print(
  table(
    meta_df$Group,
    meta_df$Time
  )
)


# ======================================================
# 9. Convert matrix to long format
# ======================================================

df_long <- as.data.frame(
  mat
) %>%

  tibble::rownames_to_column(
    "ID"
  ) %>%

  pivot_longer(

    cols = -ID,

    names_to = "ARG_type",

    values_to = "Abundance"
  ) %>%

  left_join(

    meta_df %>%
      select(
        ID,
        Group,
        Time
      ),

    by = "ID"
  ) %>%

  filter(
    !is.na(Group),
    !is.na(Time)
  )


# ======================================================
# 10. Pairwise Wilcoxon function
#
# g1 = denominator/reference group
# g2 = numerator/comparison group
#
# log2FC = log2(mean_g2 / mean_g1)
# ======================================================

pairwise_wilcox_one <- function(
    data,
    g1,
    g2,
    pseudo = 1e-6
) {


  data_sub <- data %>%

    filter(
      Group %in% c(
        g1,
        g2
      )
    )


  res <- data_sub %>%

    group_by(
      ARG_type
    ) %>%

    summarise(

      n_g1 = sum(
        Group == g1
      ),

      n_g2 = sum(
        Group == g2
      ),


      mean_g1 = mean(
        Abundance[
          Group == g1
        ],
        na.rm = TRUE
      ),

      mean_g2 = mean(
        Abundance[
          Group == g2
        ],
        na.rm = TRUE
      ),


      median_g1 = median(
        Abundance[
          Group == g1
        ],
        na.rm = TRUE
      ),

      median_g2 = median(
        Abundance[
          Group == g2
        ],
        na.rm = TRUE
      ),


      p_value = tryCatch(

        wilcox.test(

          Abundance[
            Group == g2
          ],

          Abundance[
            Group == g1
          ],

          paired = FALSE,

          exact = FALSE

        )$p.value,

        error = function(e) {
          NA_real_
        }
      ),

      .groups = "drop"
    ) %>%

    mutate(

      comparison = paste0(
        g2,
        " vs ",
        g1
      ),


      log2FC = log2(

        (
          mean_g2 +
            pseudo
        ) /

          (
            mean_g1 +
              pseudo
          )
      ),


      higher_group = case_when(

        log2FC > 0 ~ g2,

        log2FC < 0 ~ g1,

        TRUE ~ "Equal"
      )
    )


  return(
    res
  )
}


# ======================================================
# 11. Run three pairwise comparisons
#
# log2FC direction:
#
# OTM1 vs ITM:
# log2(OTM1 / ITM)
#
# OTM2 vs ITM:
# log2(OTM2 / ITM)
#
# OTM2 vs OTM1:
# log2(OTM2 / OTM1)
# ======================================================

pairwise_all <- bind_rows(


  pairwise_wilcox_one(
    df_long,
    "ITM",
    "OTM1"
  ),


  pairwise_wilcox_one(
    df_long,
    "ITM",
    "OTM2"
  ),


  pairwise_wilcox_one(
    df_long,
    "OTM1",
    "OTM2"
  )

) %>%

  filter(
    !is.na(p_value)
  ) %>%

  group_by(
    comparison
  ) %>%

  mutate(

    p_adj = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%

  ungroup()


# ======================================================
# 12. Prepare volcano plot data
# ======================================================

pairwise_plot_fdr <- pairwise_all %>%

  filter(

    !is.na(p_adj),

    !is.na(log2FC),

    is.finite(log2FC)
  ) %>%

  mutate(

    neg_log10_FDR =
      -log10(
        pmax(
          p_adj,
          1e-300
        )
      ),


    direction_FDR = case_when(


      p_adj < 0.05 &
        log2FC > 0 ~

        paste0(
          "Higher in ",
          str_replace(
            comparison,
            " vs .*",
            ""
          )
        ),


      p_adj < 0.05 &
        log2FC < 0 ~

        paste0(
          "Higher in ",
          str_replace(
            comparison,
            ".* vs ",
            ""
          )
        ),


      TRUE ~
        "Not significant"
    ),


    direction_FDR = factor(

      direction_FDR,

      levels = c(

        "Higher in ITM",

        "Higher in OTM1",

        "Higher in OTM2",

        "Not significant"
      )
    )
  )


# ======================================================
# 13. Select labels
#
# Label FDR < 0.05 ARG types.
# If >20 significant ARG types occur in one comparison,
# retain the 20 ARG types with the lowest FDR values.
# ======================================================

label_df_fdr <- pairwise_plot_fdr %>%

  filter(
    p_adj < 0.05
  ) %>%

  group_by(
    comparison
  ) %>%

  arrange(
    p_adj,
    .by_group = TRUE
  ) %>%

  slice_head(
    n = 20
  ) %>%

  ungroup()


# ======================================================
# 14. Summarize significant ARG types
# ======================================================

fdr_sig_summary <- pairwise_plot_fdr %>%

  filter(
    p_adj < 0.05
  ) %>%

  count(
    comparison,
    name = "n_FDR_lt_0.05"
  )


cat(
  "\nNumber of ARG types with FDR < 0.05:\n"
)


print(
  fdr_sig_summary
)


cat(
  "\nTop results:\n"
)


print(

  pairwise_plot_fdr %>%

    arrange(
      comparison,
      p_adj
    ) %>%

    group_by(
      comparison
    ) %>%

    slice_head(
      n = 10
    ) %>%

    ungroup()
)


# ======================================================
# 15. Export complete statistical results
# ======================================================

write.csv(

  pairwise_plot_fdr %>%

    arrange(
      comparison,
      p_adj
    ),

  file.path(
    output_dir,
    "D30D60D90_ARG_type_pairwise_Wilcoxon_FDR_all.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 16. Export significant results
# ======================================================

write.csv(

  pairwise_plot_fdr %>%

    filter(
      p_adj < 0.05
    ) %>%

    arrange(
      comparison,
      p_adj
    ),

  file.path(
    output_dir,
    "D30D60D90_ARG_type_pairwise_Wilcoxon_FDR_significant.csv"
  ),

  row.names = FALSE
)


# Export significant-count summary
write.csv(

  fdr_sig_summary,

  file.path(
    output_dir,
    "D30D60D90_ARG_type_FDR_significant_summary.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 17. Colors
# ======================================================

group_cols <- c(

  "Higher in ITM" =
    "#8ECFC9",

  "Higher in OTM1" =
    "#FFBE7A",

  "Higher in OTM2" =
    "#FA7F6F",

  "Not significant" =
    "grey75"
)


# ======================================================
# 18. Volcano plot
# ======================================================

p_volcano_fdr <- ggplot(

  pairwise_plot_fdr,

  aes(
    x = log2FC,
    y = neg_log10_FDR
  )

) +


  geom_point(

    aes(
      color = direction_FDR
    ),

    size = 2.3,

    alpha = 0.85
  ) +


  # No-change reference
  geom_vline(

    xintercept = 0,

    linetype = "dashed",

    linewidth = 0.35
  ) +


  # FDR = 0.05 threshold
  geom_hline(

    yintercept =
      -log10(
        0.05
      ),

    linetype = "dashed",

    linewidth = 0.35
  ) +


  # Significant ARG type labels
  ggrepel::geom_text_repel(

    data = label_df_fdr,

    aes(
      label = ARG_type
    ),

    size = 3,

    max.overlaps = Inf,

    box.padding = 0.35,

    point.padding = 0.25,

    segment.size = 0.25,

    min.segment.length = 0
  ) +


  facet_wrap(

    ~ comparison,

    nrow = 1,

    scales = "free_x"
  ) +


  scale_color_manual(

    values = group_cols,

    drop = FALSE
  ) +


  labs(

    title =
      "Differential ARG type abundance during D30-D90",

    subtitle =
      "Pairwise Wilcoxon test with BH/FDR correction; labels indicate FDR < 0.05",

    x =
      "log2 fold change",

    y =
      "\u2212log10(FDR)",

    color =
      ""
  ) +


  theme_bw(
    base_size = 12
  ) +


  theme(

    panel.grid.minor =
      element_blank(),

    panel.grid.major =
      element_line(
        linewidth = 0.2,
        color = "grey88"
      ),


    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5,
        size = 15
      ),

    plot.subtitle =
      element_text(
        hjust = 0.5,
        size = 10
      ),


    strip.background =
      element_rect(
        fill = "grey90",
        color = "black",
        linewidth = 0.5
      ),

    strip.text =
      element_text(
        face = "bold",
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
      "bottom",

    legend.text =
      element_text(
        size = 10
      )
  )


print(
  p_volcano_fdr
)


# ======================================================
# 19. Save volcano plot
# ======================================================

ggsave(

  filename = file.path(
    output_dir,
    "D30D60D90_ARG_type_pairwise_volcano_FDR_labeled.pdf"
  ),

  plot = p_volcano_fdr,

  width = 15,

  height = 5.5
)


# ======================================================
# 20. Finish
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "ARG type differential analysis completed.\n"
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