########################################################
## D30-D90 MRG type differential abundance analysis
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
##   corresponding MRG type/class before analysis.
##
##   Example:
##     Copper__subtype1
##     Copper__subtype2
##          ↓
##        Copper
##
## Analysis:
##   Only post-treatment samples are included:
##   D30, D60, and D90
##
##   D0 baseline samples are excluded.
##
## Pairwise comparisons:
##   OTM1 vs ITM
##   OTM2 vs ITM
##   OTM2 vs OTM1
##
## Statistical analysis:
##   Wilcoxon rank-sum test (unpaired)
##
## Multiple-testing correction:
##   Benjamini-Hochberg (BH/FDR)
##   applied separately within each pairwise comparison
##
## Fold change:
##   log2FC calculated from group mean abundances
##
## Volcano plot:
##   x-axis = log2 fold change
##   y-axis = -log10(FDR)
##   labels = up to 20 MRG types with FDR < 0.05
##
## Output:
##   results/Fig2_MRG_type_volcano/
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
  "MRG",
  "MRG_merged_subtype_abundance.tsv"
)


output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig2_MRG_type_volcano"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Post-treatment sampling times
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


# Keep D30-D90 only
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
  "MRG type differential abundance analysis\n"
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
# 8. Extract MRG type from merged-subtype name
#
# Rule:
#
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
  "MRG merged-subtype -> type aggregation\n"
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
    unique(
      feature_df$MRG_type
    )
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
# 9. Export merged-subtype -> type mapping
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
# 10. Aggregate merged subtypes to MRG type
#
# Abundance of all merged subtypes belonging to the
# same MRG type is summed for each sequencing sample.
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
# 11. Export derived MRG type abundance matrix
#
# This is generated from the core merged-subtype input
# and is not required as an independent repository input.
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
  "\nMatched D30-D90 samples:",
  length(common_samples),
  "\n"
)


# Expected:
# 30 cows × 3 post-treatment time points = 90 samples
if (length(common_samples) != 90) {

  warning(
    paste0(
      "Expected 90 D30-D90 samples, but ",
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
# 14. Merge metadata and MRG type abundance
# ======================================================

analysis_df <- mrg_long %>%

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

    Time %in%
      use_times,

    !is.na(Group)
  ) %>%

  mutate(

    Group = factor(
      Group,
      levels = c(
        "ITM",
        "OTM1",
        "OTM2"
      )
    ),

    Time = factor(
      Time,
      levels = use_times
    )
  )


# ======================================================
# 15. Check analysis dataset
# ======================================================

cat(
  "\nSamples used in D30-D90 analysis:\n"
)


print(

  analysis_df %>%

    distinct(
      ID,
      Group,
      Time
    ) %>%

    count(
      Group,
      Time
    )
)


cat(
  "\nNumber of MRG types/classes:",
  length(
    unique(
      analysis_df$MRG_type
    )
  ),
  "\n"
)


# ======================================================
# 16. Pairwise Wilcoxon test function
#
# group2 vs group1:
#
# log2FC > 0:
#   higher abundance in group2
#
# log2FC < 0:
#   higher abundance in group1
#
# D30, D60, and D90 samples are pooled.
#
# Unpaired Wilcoxon tests are intentionally retained
# to reproduce the original analysis.
# ======================================================

pairwise_mrgtype_test <- function(
    data,
    group1,
    group2,
    pseudo = 1e-6
) {


  res <- data %>%

    filter(
      Group %in%
        c(
          group1,
          group2
        )
    ) %>%

    group_by(
      MRG_type
    ) %>%

    summarise(


      n_group1 = sum(
        Group == group1 &
          !is.na(Abundance)
      ),


      n_group2 = sum(
        Group == group2 &
          !is.na(Abundance)
      ),


      mean_group1 = mean(
        Abundance[
          Group == group1
        ],
        na.rm = TRUE
      ),


      mean_group2 = mean(
        Abundance[
          Group == group2
        ],
        na.rm = TRUE
      ),


      median_group1 = median(
        Abundance[
          Group == group1
        ],
        na.rm = TRUE
      ),


      median_group2 = median(
        Abundance[
          Group == group2
        ],
        na.rm = TRUE
      ),


      p_value = tryCatch(

        wilcox.test(

          Abundance[
            Group == group2
          ],

          Abundance[
            Group == group1
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
        group2,
        " vs ",
        group1
      ),


      # Fold change based on mean abundance
      log2FC = log2(

        (
          mean_group2 +
            pseudo
        ) /

          (
            mean_group1 +
              pseudo
          )
      ),


      # BH/FDR correction is applied separately
      # within each pairwise comparison.
      p_adj = p.adjust(
        p_value,
        method = "BH"
      ),


      neg_log10_FDR = -log10(

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
            group2
          ),


        p_adj < 0.05 &
          log2FC < 0 ~

          paste0(
            "Higher in ",
            group1
          ),


        TRUE ~
          "Not significant"
      )
    ) %>%

    arrange(
      p_adj
    )


  return(
    res
  )
}


# ======================================================
# 17. Run pairwise comparisons
# ======================================================

res_OTM1_vs_ITM <- pairwise_mrgtype_test(

  data = analysis_df,

  group1 = "ITM",

  group2 = "OTM1"
)


res_OTM2_vs_ITM <- pairwise_mrgtype_test(

  data = analysis_df,

  group1 = "ITM",

  group2 = "OTM2"
)


res_OTM2_vs_OTM1 <- pairwise_mrgtype_test(

  data = analysis_df,

  group1 = "OTM1",

  group2 = "OTM2"
)


# ======================================================
# 18. Combine statistical results
# ======================================================

pairwise_mrg_D30D90 <- bind_rows(

  res_OTM1_vs_ITM,

  res_OTM2_vs_ITM,

  res_OTM2_vs_OTM1

) %>%

  filter(

    !is.na(p_value),

    !is.na(p_adj),

    !is.na(log2FC),

    is.finite(log2FC)
  ) %>%

  mutate(

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
# 19. Summarize significant MRG types
# ======================================================

fdr_sig_summary <- pairwise_mrg_D30D90 %>%

  filter(
    p_adj < 0.05
  ) %>%

  count(
    comparison,
    name = "n_FDR_lt_0.05"
  )


cat(
  "\nFDR < 0.05 summary:\n"
)


print(
  fdr_sig_summary
)


cat(
  "\nTop FDR results:\n"
)


print(

  pairwise_mrg_D30D90 %>%

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
# 20. Export complete statistical results
# ======================================================

write.csv(

  pairwise_mrg_D30D90 %>%

    arrange(
      comparison,
      p_adj
    ),

  file.path(
    output_dir,
    "D30D60D90_MRGtype_pairwise_Wilcoxon_FDR_all.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 21. Export significant results
# ======================================================

write.csv(

  pairwise_mrg_D30D90 %>%

    filter(
      p_adj < 0.05
    ) %>%

    arrange(
      comparison,
      p_adj
    ),

  file.path(
    output_dir,
    "D30D60D90_MRGtype_pairwise_Wilcoxon_FDR_lt_0.05.csv"
  ),

  row.names = FALSE
)


write.csv(

  fdr_sig_summary,

  file.path(
    output_dir,
    "D30D60D90_MRGtype_FDR_significant_summary.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 22. Select labels for volcano plot
#
# FDR < 0.05 features are eligible.
#
# If more than 20 significant MRG types occur within
# one comparison, retain the 20 with the lowest FDR.
# ======================================================

label_df_FDR <- pairwise_mrg_D30D90 %>%

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
# 23. Colors
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
# 24. Volcano plot
# ======================================================

p_volcano_FDR <- ggplot(

  pairwise_mrg_D30D90,

  aes(
    x = log2FC,
    y = neg_log10_FDR
  )

) +


  geom_point(

    aes(
      color = direction_FDR
    ),

    size = 2.8,

    alpha = 0.85
  ) +


  geom_vline(

    xintercept = 0,

    linetype = "dashed",

    linewidth = 0.35
  ) +


  geom_hline(

    yintercept =
      -log10(
        0.05
      ),

    linetype = "dashed",

    linewidth = 0.35
  ) +


  ggrepel::geom_text_repel(

    data =
      label_df_FDR,

    aes(
      label =
        MRG_type
    ),

    size = 3.2,

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

    values =
      group_cols,

    drop =
      FALSE
  ) +


  labs(

    title =
      "FDR-adjusted differential MRG types during D30-D90",

    subtitle =
      paste0(
        "Pairwise Wilcoxon test with ",
        "Benjamini-Hochberg correction; ",
        "labels indicate FDR < 0.05"
      ),

    x =
      "log2 fold change",

    y =
      "\u2212log10(FDR)",

    color =
      NULL
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
  p_volcano_FDR
)


# ======================================================
# 25. Save volcano plot
# ======================================================

ggsave(

  filename = file.path(
    output_dir,
    "Figure_D30D60D90_MRGtype_volcano_FDR_labeled.pdf"
  ),

  plot =
    p_volcano_FDR,

  width =
    15,

  height =
    5.5
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_D30D60D90_MRGtype_volcano_FDR_labeled.png"
  ),

  plot =
    p_volcano_FDR,

  width =
    15,

  height =
    5.5,

  dpi =
    600
)


# ======================================================
# 26. Finish
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "D30-D90 MRG type differential analysis completed.\n"
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