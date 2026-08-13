########################################################
## Associations between functional MetaCyc pathway
## categories and metal-related MRG abundance
##
## Figure 5
##
## Input:
##   metadata/metadata.csv
##
##   processed_data/HUMAnN/
##     pathabundance_relab_unstratified.tsv
##
##   processed_data/MRG/
##     MRG_type_abundance.tsv
##
##   results/Fig5_HUMAnN_pathway_volcano/
##     stat_MetaCyc_pathways_all_comparisons_D30D60D90.csv
##
## Samples:
##   D30 + D60 + D90
##
## Pathway selection:
##   MetaCyc pathways significant at FDR < 0.05
##   in at least one of the three pairwise treatment
##   comparisons are retained.
##
## Functional categories:
##   1. Carbohydrate metabolism / fermentation
##   2. Amino acid / polyamine metabolism
##   3. Nucleotide metabolism
##   4. Cofactor / redox metabolism
##   5. Cell envelope / membrane
##   6. Energy metabolism
##
## For each sample:
##   pathway-category abundance =
##   sum of retained pathway abundances assigned
##   to the corresponding functional category.
##
## MRG targets:
##   Copper
##   Iron
##   Zinc
##
## Analysis:
##   Group-specific linear regression:
##
##     MRG abundance ~ pathway-category abundance
##
##   performed separately for:
##     ITM
##     OTM1
##     OTM2
##
## Multiple-testing correction:
##   BH/FDR correction across the three treatment
##   groups within each Functional category × Metal
##   combination.
##
## Transformation:
##   log10(abundance + 1e-8)
##
## Output:
##   results/Fig5_pathway_category_MRG_association/
########################################################


# ======================================================
# 1. Load packages
# ======================================================

library(tidyverse)
library(patchwork)


# ======================================================
# 2. Project paths
#
# Run this script from the repository root:
# DairyCow-TraceMineral-Resistome/
# ======================================================

PROJECT_DIR <- "."


pathway_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "HUMAnN",
  "pathabundance_relab_unstratified.tsv"
)


meta_path <- file.path(
  PROJECT_DIR,
  "metadata",
  "metadata.csv"
)


diff_path <- file.path(
  PROJECT_DIR,
  "results",
  "Fig5_HUMAnN_pathway_volcano",
  "stat_MetaCyc_pathways_all_comparisons_D30D60D90.csv"
)


mrg_type_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "MRG",
  "MRG_type_abundance.tsv"
)


output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig5_pathway_category_MRG_association"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ======================================================
# 3. Analysis settings
# ======================================================

use_times <- c(
  "D30",
  "D60",
  "D90"
)


group_levels <- c(
  "ITM",
  "OTM1",
  "OTM2"
)


metal_targets <- c(
  "Copper",
  "Iron",
  "Zinc"
)


# Log10 transformation
log_transform <- TRUE


pseudo_count <- 1e-8


# ======================================================
# 4. Read metadata
# ======================================================

meta <- read.csv(
  meta_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


required_meta_cols <- c(
  "ID",
  "Time",
  "Group"
)


missing_meta_cols <- setdiff(
  required_meta_cols,
  colnames(meta)
)


if (length(missing_meta_cols) > 0) {

  stop(
    paste0(
      "Missing metadata column(s): ",
      paste(
        missing_meta_cols,
        collapse = ", "
      )
    )
  )
}


meta <- meta %>%

  mutate(

    ID = as.character(ID),

    Time = as.character(Time),

    Group = as.character(Group)
  )


# ======================================================
# 5. Check metadata
# ======================================================

unexpected_groups <- setdiff(
  unique(meta$Group),
  group_levels
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


expected_times <- c(
  "D0",
  "D30",
  "D60",
  "D90"
)


unexpected_times <- setdiff(
  unique(meta$Time),
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


if (anyDuplicated(meta$ID) > 0) {

  stop(
    "Duplicated sequencing sample IDs were found in metadata.csv."
  )
}


meta$Group <- factor(
  meta$Group,
  levels = group_levels
)


meta$Time <- factor(
  meta$Time,
  levels = expected_times
)


# Retain D30-D90
meta_post <- meta %>%

  filter(
    Time %in% use_times
  ) %>%

  droplevels()


rownames(meta_post) <- meta_post$ID


cat(
  "\n========================================\n"
)

cat(
  "Pathway-category × MRG association analysis\n"
)

cat(
  "========================================\n"
)


cat(
  "\nSamples retained:\n"
)


print(
  table(
    meta_post$Group,
    meta_post$Time
  )
)


cat(
  "\nTotal post-baseline samples:",
  nrow(meta_post),
  "\n"
)


if (nrow(meta_post) != 90) {

  warning(
    paste0(
      "Expected 90 D30-D90 samples, but ",
      nrow(meta_post),
      " were retained."
    )
  )
}


# ======================================================
# 6. Read HUMAnN MetaCyc pathway abundance
# ======================================================

pathway <- read.delim(

  pathway_path,

  sep = "\t",

  header = TRUE,

  check.names = FALSE,

  comment.char = "",

  quote = "",

  stringsAsFactors = FALSE
)


colnames(pathway)[1] <- "Pathway"


pathway$Pathway <- as.character(
  pathway$Pathway
)


# Remove HUMAnN special rows
pathway <- pathway %>%

  filter(
    !Pathway %in%
      c(
        "UNMAPPED",
        "UNINTEGRATED",
        "UNGROUPED"
      )
  )


# ======================================================
# 7. Convert pathway abundance columns to numeric
# ======================================================

pathway_sample_cols <- setdiff(
  colnames(pathway),
  "Pathway"
)


pathway[
  pathway_sample_cols
] <- lapply(

  pathway[
    pathway_sample_cols
  ],

  function(x) {

    x <- as.numeric(
      as.character(x)
    )


    x[
      is.na(x)
    ] <- 0


    x[
      x < 0
    ] <- 0


    x
  }
)


# ======================================================
# 8. Merge duplicated pathways
# ======================================================

pathway <- pathway %>%

  group_by(
    Pathway
  ) %>%

  summarise(

    across(
      all_of(
        pathway_sample_cols
      ),
      ~ sum(
        .x,
        na.rm = TRUE
      )
    ),

    .groups = "drop"
  )


# Pathway × sample
pathway_mat <- pathway %>%

  column_to_rownames(
    "Pathway"
  ) %>%

  as.matrix()


storage.mode(
  pathway_mat
) <- "numeric"


# sample × pathway
pathway_mat <- t(
  pathway_mat
)


# ======================================================
# 9. Clean HUMAnN sample names
# ======================================================

original_pathway_sample_names <- rownames(
  pathway_mat
)


clean_pathway_sample_names <- original_pathway_sample_names


clean_pathway_sample_names <- gsub(
  "_Abundance$",
  "",
  clean_pathway_sample_names
)


clean_pathway_sample_names <- gsub(
  "_merged.*$",
  "",
  clean_pathway_sample_names
)


clean_pathway_sample_names <- gsub(
  "\\.tsv$",
  "",
  clean_pathway_sample_names
)


rownames(pathway_mat) <- clean_pathway_sample_names


if (
  anyDuplicated(
    rownames(pathway_mat)
  ) > 0
) {

  stop(
    paste0(
      "Duplicated sample IDs were generated after ",
      "cleaning HUMAnN sample names."
    )
  )
}


cat(
  "\nHUMAnN pathway matrix:\n"
)


cat(
  "Samples:",
  nrow(pathway_mat),
  "\n"
)


cat(
  "Pathways:",
  ncol(pathway_mat),
  "\n"
)


# ======================================================
# 10. Read differential pathway results from script 22
# ======================================================

diff_all <- read.csv(
  diff_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


required_diff_cols <- c(
  "Pathway",
  "Description",
  "FDR"
)


missing_diff_cols <- setdiff(
  required_diff_cols,
  colnames(diff_all)
)


if (length(missing_diff_cols) > 0) {

  stop(
    paste0(
      "Missing differential-result column(s): ",
      paste(
        missing_diff_cols,
        collapse = ", "
      )
    )
  )
}


# ======================================================
# 11. Functional-category assignment
#
# Retained from the original analysis.
#
# First matching category is retained.
# ======================================================

assign_function_category <- function(
    description
) {


  x <- tolower(
    as.character(
      description
    )
  )


  case_when(


    str_detect(

      x,

      paste0(

        "glycolysis|glucose|gluconeogenesis|pentose|",

        "starch|glycogen|sucrose|fructose|mannose|",

        "galactose|xylose|arabinose|cellulose|",

        "hemicellulose|carbohydrate|sugar|",

        "fermentation|lactate|lactic|butanol|",

        "ethanol|acetate|propionate|butyrate|",

        "pyruvate|inositol|neuraminate|glucosamine|",

        "mannosamine|entner|doudoroff|",

        "bifidobacterium shunt"
      )

    ) ~

      "Carbohydrate metabolism / fermentation",


    str_detect(

      x,

      paste0(

        "amino acid|alanine|arginine|aspartate|",

        "asparagine|cysteine|glutamate|glutamine|",

        "glycine|histidine|isoleucine|leucine|",

        "lysine|methionine|phenylalanine|proline|",

        "serine|threonine|tryptophan|tyrosine|",

        "valine|ornithine|urea|polyamine|",

        "putrescine|citrulline|sulfur amino"
      )

    ) ~

      "Amino acid / polyamine metabolism",


    str_detect(

      x,

      paste0(

        "purine|pyrimidine|nucleotide|nucleoside|",

        "adenosine|guanosine|cytidine|uridine|",

        "thymidine|inosine|ump|dna|rna|trna|",

        "queuosine|ribonucleotide|deoxyribonucleotide"
      )

    ) ~

      "Nucleotide metabolism",


    str_detect(

      x,

      paste0(

        "folate|cobalamin|vitamin|biotin|riboflavin|",

        "thiamin|thiamine|pyridoxal|pantothenate|",

        "coa|coenzyme|nad|nadh|nadp|fad|heme|",

        "menaquinone|menaquinol|quinone|tetrapyrrole"
      )

    ) ~

      "Cofactor / redox metabolism",


    str_detect(

      x,

      paste0(

        "cell wall|cell envelope|membrane|",

        "peptidoglycan|lipopolysaccharide|lps|",

        "o-antigen|capsule|exopolysaccharide|",

        "outer membrane|murein|udp-n-acetyl"
      )

    ) ~

      "Cell envelope / membrane",


    str_detect(

      x,

      paste0(

        "tca|tricarboxylic|citric acid cycle|",

        "respiration|electron transfer|",

        "electron transport|oxidative phosphorylation|",

        "atp|glyoxylate"
      )

    ) ~

      "Energy metabolism",


    TRUE ~

      "Other / unclear"
  )
}


diff_all <- diff_all %>%

  mutate(

    Functional_category =
      assign_function_category(
        Description
      )
  )


# ======================================================
# 12. Six functional categories used in the figure
# ======================================================

focus_categories <- c(

  "Carbohydrate metabolism / fermentation",

  "Amino acid / polyamine metabolism",

  "Nucleotide metabolism",

  "Cofactor / redox metabolism",

  "Cell envelope / membrane",

  "Energy metabolism"
)


category_labels <- setNames(
  focus_categories,
  focus_categories
)


# ======================================================
# 13. Select significant pathways
#
# A pathway is retained when:
#
#   FDR < 0.05 in at least one pairwise comparison
#
# and it belongs to one of the six selected categories.
#
# The union across all pairwise comparisons is used.
# ======================================================

selected_info <- diff_all %>%

  filter(

    !is.na(FDR),

    FDR < 0.05,

    Functional_category %in%
      focus_categories
  ) %>%

  select(
    Pathway,
    Functional_category
  ) %>%

  distinct()


cat(
  "\nSelected FDR-significant pathways by category:\n"
)


print(
  table(
    selected_info$Functional_category
  )
)


# Export exact pathway list used
write.csv(

  selected_info,

  file.path(
    output_dir,
    "stat_selected_MetaCyc_pathways_by_functional_category.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 14. Match selected pathways to HUMAnN matrix
# ======================================================

selected_pathways <- intersect(

  selected_info$Pathway,

  colnames(
    pathway_mat
  )
)


cat(
  "\nSelected pathways matched to HUMAnN matrix:",
  length(selected_pathways),
  "\n"
)


if (length(selected_pathways) == 0) {

  stop(
    paste0(
      "No selected differential pathways matched ",
      "the HUMAnN pathway abundance matrix."
    )
  )
}


missing_selected_pathways <- setdiff(
  selected_info$Pathway,
  colnames(pathway_mat)
)


if (length(missing_selected_pathways) > 0) {

  cat(
    "\nSelected pathways missing from HUMAnN matrix:\n"
  )

  print(
    missing_selected_pathways
  )
}


# ======================================================
# 15. Match HUMAnN samples to post-baseline metadata
# ======================================================

common_pathway_samples <- intersect(

  rownames(
    pathway_mat
  ),

  meta_post$ID
)


pathway_post <- pathway_mat[
  common_pathway_samples,
  selected_pathways,
  drop = FALSE
]


meta_pathway <- meta_post[
  match(
    common_pathway_samples,
    meta_post$ID
  ),
  ,
  drop = FALSE
]


# ======================================================
# 16. Calculate category-level pathway abundance
#
# Category abundance =
# sum of all selected significant pathways assigned
# to the category for each sample.
# ======================================================

selected_info_matched <- selected_info %>%

  filter(
    Pathway %in%
      selected_pathways
  )


pathway_long <- as.data.frame(
  pathway_post
) %>%

  rownames_to_column(
    "ID"
  ) %>%

  pivot_longer(

    cols =
      -ID,

    names_to =
      "Pathway",

    values_to =
      "Pathway_abundance"
  ) %>%

  left_join(

    selected_info_matched,

    by =
      "Pathway"
  )


category_abundance_long <- pathway_long %>%

  group_by(
    ID,
    Functional_category
  ) %>%

  summarise(

    Category_abundance =
      sum(
        Pathway_abundance,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  )


# Complete missing category × sample combinations with zero
category_abundance_long <- tidyr::complete(

  category_abundance_long,

  ID =
    common_pathway_samples,

  Functional_category =
    focus_categories,

  fill =
    list(
      Category_abundance = 0
    )
)


category_mat <- category_abundance_long %>%

  pivot_wider(

    names_from =
      Functional_category,

    values_from =
      Category_abundance,

    values_fill =
      0
  ) %>%

  column_to_rownames(
    "ID"
  ) %>%

  as.matrix()


storage.mode(
  category_mat
) <- "numeric"


cat(
  "\nCategory-level pathway matrix:\n"
)


cat(
  "Samples:",
  nrow(category_mat),
  "\n"
)


cat(
  "Functional categories:",
  ncol(category_mat),
  "\n"
)


print(
  colnames(
    category_mat
  )
)


# ======================================================
# 17. Read MRG type abundance table
#
# Expected first column:
#   MRG type/category
#
# Expected categories include:
#   Copper
#   Iron
#   Zinc
# ======================================================

mrg <- read.delim(

  mrg_type_path,

  sep = "\t",

  header = TRUE,

  check.names = FALSE,

  comment.char = "",

  quote = "",

  stringsAsFactors = FALSE
)


colnames(mrg)[1] <- "MRG_type"


mrg$MRG_type <- as.character(
  mrg$MRG_type
)


mrg_sample_cols <- setdiff(
  colnames(mrg),
  "MRG_type"
)


mrg[
  mrg_sample_cols
] <- lapply(

  mrg[
    mrg_sample_cols
  ],

  function(x) {

    x <- as.numeric(
      as.character(x)
    )


    x[
      is.na(x)
    ] <- 0


    x[
      x < 0
    ] <- 0


    x
  }
)


# ======================================================
# 18. Aggregate duplicated MRG types
# ======================================================

mrg_type_df <- mrg %>%

  group_by(
    MRG_type
  ) %>%

  summarise(

    across(
      all_of(
        mrg_sample_cols
      ),
      ~ sum(
        .x,
        na.rm = TRUE
      )
    ),

    .groups =
      "drop"
  )


# type × sample
mrg_mat <- mrg_type_df %>%

  column_to_rownames(
    "MRG_type"
  ) %>%

  as.matrix()


storage.mode(
  mrg_mat
) <- "numeric"


# sample × MRG type
mrg_mat <- t(
  mrg_mat
)


# ======================================================
# 19. Check Copper / Iron / Zinc
# ======================================================

missing_metals <- setdiff(

  metal_targets,

  colnames(
    mrg_mat
  )
)


if (length(missing_metals) > 0) {

  stop(
    paste0(
      "Required MRG type(s) missing: ",
      paste(
        missing_metals,
        collapse = ", "
      )
    )
  )
}


mrg_targets <- mrg_mat[
  ,
  metal_targets,
  drop = FALSE
]


cat(
  "\nMRG targets used:\n"
)


print(
  colnames(
    mrg_targets
  )
)


# ======================================================
# 20. Match pathway, MRG, and metadata samples
# ======================================================

common_samples <- Reduce(

  intersect,

  list(

    rownames(
      category_mat
    ),

    rownames(
      mrg_targets
    ),

    meta_post$ID
  )
)


cat(
  "\nMatched post-baseline samples:",
  length(common_samples),
  "\n"
)


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
    "Too few common samples for pathway-MRG association analysis."
  )
}


category_mat <- category_mat[
  common_samples,
  ,
  drop = FALSE
]


mrg_targets <- mrg_targets[
  common_samples,
  ,
  drop = FALSE
]


meta_matched <- meta_post[
  match(
    common_samples,
    meta_post$ID
  ),
  ,
  drop = FALSE
]


stopifnot(
  all(
    rownames(category_mat) ==
      rownames(mrg_targets)
  )
)


stopifnot(
  all(
    rownames(category_mat) ==
      meta_matched$ID
  )
)


cat(
  "\nFinal matched sample distribution:\n"
)


print(
  table(
    meta_matched$Group,
    meta_matched$Time
  )
)


# ======================================================
# 21. Prepare long-format pathway data
# ======================================================

category_long <- as.data.frame(
  category_mat
) %>%

  rownames_to_column(
    "ID"
  ) %>%

  pivot_longer(

    cols =
      -ID,

    names_to =
      "Functional_category",

    values_to =
      "Category_abundance"
  )


# ======================================================
# 22. Prepare long-format MRG data
# ======================================================

mrg_long <- as.data.frame(
  mrg_targets
) %>%

  rownames_to_column(
    "ID"
  ) %>%

  pivot_longer(

    cols =
      -ID,

    names_to =
      "Metal",

    values_to =
      "MRG_abundance"
  )


# ======================================================
# 23. Build final analysis dataset
# ======================================================

plot_all <- category_long %>%

  inner_join(
    mrg_long,
    by = "ID"
  ) %>%

  left_join(

    meta_matched %>%

      select(
        ID,
        Group,
        Time
      ),

    by =
      "ID"
  ) %>%

  filter(

    Functional_category %in%
      focus_categories,

    Metal %in%
      metal_targets
  ) %>%

  mutate(

    Functional_category = factor(
      Functional_category,
      levels =
        focus_categories
    ),

    Metal = factor(
      Metal,
      levels =
        metal_targets
    ),

    Group = factor(
      Group,
      levels =
        group_levels
    )
  )


# ======================================================
# 24. Apply abundance transformation
# ======================================================

if (log_transform) {


  plot_all <- plot_all %>%

    mutate(

      X_value =
        log10(
          Category_abundance +
            pseudo_count
        ),

      Y_value =
        log10(
          MRG_abundance +
            pseudo_count
        )
    )


  x_lab_text <-
    "Pathway-category abundance (log10)"


  y_lab_text <-
    "MRG abundance (log10)"


} else {


  plot_all <- plot_all %>%

    mutate(

      X_value =
        Category_abundance,

      Y_value =
        MRG_abundance
    )


  x_lab_text <-
    "Pathway-category abundance"


  y_lab_text <-
    "MRG abundance"
}


# Export exact analysis dataset
write.csv(

  plot_all,

  file.path(
    output_dir,
    "plot_pathway_category_MRG_association_data.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 25. Group-specific linear regression
#
# Model:
#
#   Y_value ~ X_value
#
# separately within:
#
# Functional category × Metal × Group
# ======================================================

lm_stats <- plot_all %>%

  group_by(
    Functional_category,
    Metal,
    Group
  ) %>%

  group_modify(

    ~ {


      data_use <- .x %>%

        filter(
          is.finite(X_value),
          is.finite(Y_value)
        )


      if (

        nrow(data_use) < 6 ||

        sd(
          data_use$X_value,
          na.rm = TRUE
        ) == 0 ||

        sd(
          data_use$Y_value,
          na.rm = TRUE
        ) == 0

      ) {


        return(

          tibble(

            N =
              nrow(data_use),

            R2 =
              NA_real_,

            P_value =
              NA_real_
          )
        )
      }


      fit <- lm(

        Y_value ~
          X_value,

        data =
          data_use
      )


      model_summary <- summary(
        fit
      )


      tibble(

        N =
          nrow(data_use),

        R2 =
          model_summary$r.squared,

        P_value =
          model_summary$coefficients[
            "X_value",
            "Pr(>|t|)"
          ]
      )
    }
  ) %>%

  ungroup()


# ======================================================
# 26. BH/FDR correction
#
# Correction is performed across ITM, OTM1, and OTM2
# within each Functional category × Metal combination.
# ======================================================

lm_stats <- lm_stats %>%

  group_by(
    Functional_category,
    Metal
  ) %>%

  mutate(

    FDR =
      p.adjust(
        P_value,
        method = "BH"
      )
  ) %>%

  ungroup() %>%

  mutate(

    Label = case_when(

      is.na(R2) ~

        paste0(
          as.character(Group),
          ": NA"
        ),

      TRUE ~

        paste0(

          as.character(Group),

          ": R\u00b2=",

          sprintf(
            "%.2f",
            R2
          ),

          ", FDR=",

          ifelse(

            FDR < 0.001,

            "<0.001",

            sprintf(
              "%.3f",
              FDR
            )
          )
        )
    )
  )


cat(
  "\nLinear model statistics:\n"
)


print(

  lm_stats %>%

    arrange(
      Functional_category,
      Metal,
      Group
    )
)


# ======================================================
# 27. Export linear-model statistics
# ======================================================

write.csv(

  lm_stats,

  file.path(
    output_dir,
    "stat_pathway_category_MRG_group_specific_linear_models.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 28. Prepare annotation positions
# ======================================================

annotation_ranges <- plot_all %>%

  group_by(
    Functional_category,
    Metal
  ) %>%

  summarise(

    x_min =
      min(
        X_value,
        na.rm = TRUE
      ),

    x_max =
      max(
        X_value,
        na.rm = TRUE
      ),

    y_min =
      min(
        Y_value,
        na.rm = TRUE
      ),

    y_max =
      max(
        Y_value,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  )


anno_df <- annotation_ranges %>%

  left_join(

    lm_stats,

    by = c(
      "Functional_category",
      "Metal"
    )
  ) %>%

  mutate(

    x_range =
      x_max -
      x_min,

    y_range =
      y_max -
      y_min,


    x_range =
      ifelse(
        !is.finite(x_range) |
          x_range == 0,
        1,
        x_range
      ),


    y_range =
      ifelse(
        !is.finite(y_range) |
          y_range == 0,
        1,
        y_range
      ),


    x_pos =
      x_min +
      0.03 *
      x_range,


    y_pos =
      y_max -

      case_when(

        Group ==
          "ITM" ~

          0.06 *
          y_range,

        Group ==
          "OTM1" ~

          0.18 *
          y_range,

        Group ==
          "OTM2" ~

          0.30 *
          y_range,

        TRUE ~

          0.06 *
          y_range
      )
  )


# ======================================================
# 29. Colors
# ======================================================

group_colors <- c(

  ITM =
    "#8ECFC9",

  OTM1 =
    "#FFBE7A",

  OTM2 =
    "#FA7F6F"
)


# ======================================================
# 30. Plot function
#
# One main panel = one functional category
#
# Each panel contains three facets:
#   Copper
#   Iron
#   Zinc
# ======================================================

plot_one_category <- function(
    category_name,
    tag_label
) {


  data_plot <- plot_all %>%

    filter(
      Functional_category ==
        category_name
    )


  data_anno <- anno_df %>%

    filter(
      Functional_category ==
        category_name
    )


  # ----------------------------------------------------
  # Handle a functional category containing no selected
  # pathways
  # ----------------------------------------------------

  if (
    nrow(
      data_plot
    ) == 0
  ) {


    return(

      ggplot() +

        annotate(

          "text",

          x = 0,

          y = 0,

          label =
            "No selected pathways",

          fontface =
            "italic",

          color =
            "grey40"
        ) +

        xlim(
          -1,
          1
        ) +

        ylim(
          -1,
          1
        ) +

        labs(

          title =
            category_labels[
              category_name
            ],

          tag =
            tag_label
        ) +

        theme_void() +

        theme(

          plot.title =
            element_text(
              face = "bold",
              hjust = 0.5,
              size = 10
            ),

          plot.tag =
            element_text(
              face = "bold",
              size = 12
            )
        )
    )
  }


  p <- ggplot(

    data_plot,

    aes(

      x =
        X_value,

      y =
        Y_value,

      color =
        Group
    )

  ) +


    geom_point(

      size =
        1.8,

      alpha =
        0.70
    ) +


    geom_smooth(

      method =
        "lm",

      se =
        FALSE,

      linewidth =
        0.7,

      na.rm =
        TRUE
    ) +


    geom_text(

      data =
        data_anno,

      aes(

        x =
          x_pos,

        y =
          y_pos,

        label =
          Label,

        color =
          Group
      ),

      inherit.aes =
        FALSE,

      hjust =
        0,

      size =
        2.3,

      show.legend =
        FALSE
    ) +


    facet_wrap(

      ~ Metal,

      nrow =
        1,

      scales =
        "free_y"
    ) +


    scale_color_manual(

      values =
        group_colors,

      limits =
        group_levels,

      drop =
        FALSE
    ) +


    labs(

      x =
        x_lab_text,

      y =
        y_lab_text,

      title =
        category_labels[
          category_name
        ],

      tag =
        tag_label,

      color =
        NULL
    ) +


    theme_bw(
      base_size = 9
    ) +


    theme(

      panel.grid.major =
        element_line(
          linewidth = 0.20,
          color = "grey90"
        ),

      panel.grid.minor =
        element_blank(),


      axis.text =
        element_text(
          color = "black",
          size = 7
        ),


      axis.title =
        element_text(
          face = "bold",
          size = 8
        ),


      strip.text =
        element_text(
          face = "bold",
          size = 8
        ),


      legend.title =
        element_blank(),


      legend.position =
        "bottom",


      plot.title =
        element_text(
          face = "bold",
          hjust = 0.5,
          size = 10
        ),


      plot.tag =
        element_text(
          face = "bold",
          size = 12
        )
    )


  return(
    p
  )
}


# ======================================================
# 31. Generate six functional-category panels
# ======================================================

p1 <- plot_one_category(

  "Carbohydrate metabolism / fermentation",

  "A"
)


p2 <- plot_one_category(

  "Amino acid / polyamine metabolism",

  "B"
)


p3 <- plot_one_category(

  "Nucleotide metabolism",

  "C"
)


p4 <- plot_one_category(

  "Cofactor / redox metabolism",

  "D"
)


p5 <- plot_one_category(

  "Cell envelope / membrane",

  "E"
)


p6 <- plot_one_category(

  "Energy metabolism",

  "F"
)


# ======================================================
# 32. Combine as 2 × 3
# ======================================================

combined_scatter <- (

  p1 |
    p2 |
    p3

) / (

  p4 |
    p5 |
    p6

) +


  patchwork::plot_layout(

    guides =
      "collect",

    heights =
      c(
        1,
        1
      )
  ) &


  theme(

    legend.position =
      "bottom"
  )


combined_scatter <- combined_scatter +

  patchwork::plot_annotation(

    title =
      paste0(
        "Associations between pathway categories ",
        "and metal-related MRGs"
      ),

    subtitle =
      "D30 + D60 + D90; group-specific linear fits",

    theme =
      theme(

        plot.title =
          element_text(
            face = "bold",
            hjust = 0.5,
            size = 14
          ),

        plot.subtitle =
          element_text(
            hjust = 0.5,
            size = 11
          )
      )
  )


print(
  combined_scatter
)


# ======================================================
# 33. Save combined figure
# ======================================================

ggsave(

  filename = file.path(
    output_dir,
    "Figure_pathway_category_MRG_association_2x3.pdf"
  ),

  plot =
    combined_scatter,

  width =
    16,

  height =
    10
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_pathway_category_MRG_association_2x3.png"
  ),

  plot =
    combined_scatter,

  width =
    16,

  height =
    10,

  dpi =
    600
)


# ======================================================
# 34. Finish
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "Pathway-category × MRG association analysis completed.\n"
)

cat(
  "Samples used: D30 + D60 + D90.\n"
)

cat(
  "MRG targets: Copper, Iron, Zinc.\n"
)

cat(
  "Group-specific linear regression was used.\n"
)

cat(
  paste0(
    "BH/FDR correction was applied across the three ",
    "groups within each pathway category × metal combination.\n"
  )
)

cat(
  "Log10 transformation:",
  log_transform,
  "\n"
)

cat(
  "Pathway input:\n"
)

cat(
  pathway_path,
  "\n"
)

cat(
  "Differential pathway input:\n"
)

cat(
  diff_path,
  "\n"
)

cat(
  "MRG type input:\n"
)

cat(
  mrg_type_path,
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