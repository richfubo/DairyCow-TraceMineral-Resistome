########################################################
## Functional categories of differential MetaCyc pathways
##
## Figure 5
##
## Input:
##   results/Fig5_HUMAnN_pathway_volcano/
##     stat_MetaCyc_pathways_all_comparisons_D30D60D90.csv
##
## Source analysis:
##   D30 + D60 + D90 combined
##
## Pairwise comparisons:
##   OTM1 vs ITM
##   OTM2 vs ITM
##   OTM2 vs OTM1
##
## Differential-pathway criterion:
##   FDR < 0.05
##
## Functional classification:
##   Significant MetaCyc pathways are assigned to broad
##   functional categories according to predefined
##   keyword patterns in pathway descriptions.
##
## Plot:
##   x-axis = number of differential pathways
##   y-axis = functional category
##   facets = pairwise comparisons
##   fill = group in which pathway abundance is higher
##
## Output:
##   results/Fig5_HUMAnN_pathway_functional_category/
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


diff_path <- file.path(
  PROJECT_DIR,
  "results",
  "Fig5_HUMAnN_pathway_volcano",
  "stat_MetaCyc_pathways_all_comparisons_D30D60D90.csv"
)


output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig5_HUMAnN_pathway_functional_category"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ======================================================
# 3. Plot setting
#
# "fixed":
#   same x-axis among comparisons
#
# "free":
#   independent x-axis for each comparison
# ======================================================

axis_mode <- "fixed"

# Alternative:
# axis_mode <- "free"


# ======================================================
# 4. Read differential pathway results
# ======================================================

diff_all <- read.csv(
  diff_path,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


required_cols <- c(
  "Comparison",
  "Group_ref",
  "Group_test",
  "Description",
  "Log2FC",
  "P_value",
  "FDR"
)


missing_cols <- setdiff(
  required_cols,
  colnames(diff_all)
)


if (length(missing_cols) > 0) {

  stop(
    paste0(
      "Missing column(s) in differential pathway file: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  )
}


cat(
  "\n========================================\n"
)

cat(
  "MetaCyc pathway functional-category analysis\n"
)

cat(
  "========================================\n"
)


cat(
  "\nTotal rows in differential result:",
  nrow(diff_all),
  "\n"
)


cat(
  "\nComparisons in input file:\n"
)


print(
  table(
    diff_all$Comparison
  )
)


# ======================================================
# 5. Comparison order
# ======================================================

comparison_levels <- c(
  "OTM1 vs ITM",
  "OTM2 vs ITM",
  "OTM2 vs OTM1"
)


unexpected_comparisons <- setdiff(
  unique(diff_all$Comparison),
  comparison_levels
)


if (length(unexpected_comparisons) > 0) {

  warning(
    paste0(
      "Unexpected comparison(s) found: ",
      paste(
        unexpected_comparisons,
        collapse = ", "
      )
    )
  )
}


diff_all$Comparison <- factor(
  diff_all$Comparison,
  levels = comparison_levels
)


# ======================================================
# 6. Functional-category assignment
#
# Classification is based on keyword matching against
# the MetaCyc pathway Description field.
#
# The first matching category is retained.
# Therefore, category order below is important.
# ======================================================

assign_function_category <- function(
    description
) {


  x <- tolower(
    as.character(
      description
    )
  )


  dplyr::case_when(


    # --------------------------------------------------
    # Carbohydrate metabolism / fermentation
    # --------------------------------------------------

    stringr::str_detect(

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


    # --------------------------------------------------
    # Amino acid / polyamine metabolism
    # --------------------------------------------------

    stringr::str_detect(

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


    # --------------------------------------------------
    # Nucleotide metabolism
    # --------------------------------------------------

    stringr::str_detect(

      x,

      paste0(

        "purine|pyrimidine|nucleotide|nucleoside|",

        "adenosine|guanosine|cytidine|uridine|",

        "thymidine|inosine|ump|dna|rna|trna|",

        "queuosine|ribonucleotide|deoxyribonucleotide"
      )

    ) ~

      "Nucleotide metabolism",


    # --------------------------------------------------
    # Cofactor / redox metabolism
    # --------------------------------------------------

    stringr::str_detect(

      x,

      paste0(

        "folate|cobalamin|vitamin|biotin|",

        "riboflavin|thiamin|thiamine|pyridoxal|",

        "pantothenate|coa|coenzyme|nad|nadh|nadp|",

        "fad|heme|menaquinone|menaquinol|quinone|",

        "tetrapyrrole"
      )

    ) ~

      "Cofactor / redox metabolism",


    # --------------------------------------------------
    # Lipid metabolism
    # --------------------------------------------------

    stringr::str_detect(

      x,

      paste0(

        "fatty acid|lipid|phospholipid|",

        "glycerolipid|beta-oxidation|β-oxidation"
      )

    ) ~

      "Lipid metabolism",


    # --------------------------------------------------
    # Energy metabolism
    # --------------------------------------------------

    stringr::str_detect(

      x,

      paste0(

        "tca|tricarboxylic|citric acid cycle|",

        "respiration|electron transfer|",

        "electron transport|oxidative phosphorylation|",

        "atp|glyoxylate"
      )

    ) ~

      "Energy metabolism",


    # --------------------------------------------------
    # Cell envelope / membrane
    # --------------------------------------------------

    stringr::str_detect(

      x,

      paste0(

        "cell wall|cell envelope|membrane|",

        "peptidoglycan|lipopolysaccharide|lps|",

        "o-antigen|capsule|exopolysaccharide|",

        "outer membrane|murein|udp-n-acetyl"
      )

    ) ~

      "Cell envelope / membrane",


    # --------------------------------------------------
    # Secondary metabolism / degradation
    # --------------------------------------------------

    stringr::str_detect(

      x,

      paste0(

        "degradation|catabolism|aromatic|benzoate|",

        "toluene|xylene|xenobiotic|",

        "secondary metabolite|acetylene|",

        "formaldehyde|chitin"
      )

    ) ~

      "Secondary metabolism / degradation",


    # --------------------------------------------------
    # Other biosynthesis
    # --------------------------------------------------

    stringr::str_detect(

      x,

      "biosynthesis|superpathway of .* biosynthesis|formation"

    ) ~

      "Biosynthesis, other",


    # --------------------------------------------------
    # Unassigned
    # --------------------------------------------------

    TRUE ~

      "Other / unclear"
  )
}


# ======================================================
# 7. Assign functional categories
# ======================================================

diff_all <- diff_all %>%

  mutate(

    Functional_category =
      assign_function_category(
        Description
      )
  )


# ======================================================
# 8. Determine enriched group
#
# Log2FC was defined in script 22 as:
#
#   group_test / group_ref
#
# Therefore:
#
#   Log2FC > 0 -> enriched in Group_test
#   Log2FC < 0 -> enriched in Group_ref
# ======================================================

diff_all <- diff_all %>%

  mutate(

    Enriched_group = case_when(

      FDR < 0.05 &
        Log2FC < 0 ~

        Group_ref,

      FDR < 0.05 &
        Log2FC > 0 ~

        Group_test,

      TRUE ~

        "Not significant"
    )
  )


# ======================================================
# 9. Retain FDR-significant pathways
# ======================================================

sig_df <- diff_all %>%

  filter(
    !is.na(FDR),
    FDR < 0.05
  )


cat(
  "\nNumber of FDR-significant pathways by comparison:\n"
)


print(
  table(
    sig_df$Comparison
  )
)


cat(
  "\nDirection of FDR-significant pathways:\n"
)


print(
  table(
    sig_df$Comparison,
    sig_df$Enriched_group
  )
)


# ======================================================
# 10. Export pathway-level functional assignments
#
# This table shows exactly how every significant
# pathway was assigned to a broad category.
# ======================================================

write.csv(

  sig_df %>%

    select(

      Comparison,

      Pathway,

      Pathway_ID,

      Description,

      Group_ref,

      Group_test,

      Log2FC,

      P_value,

      FDR,

      Enriched_group,

      Functional_category
    ) %>%

    arrange(
      Comparison,
      Functional_category,
      FDR
    ),

  file.path(
    output_dir,
    "stat_MetaCyc_significant_pathway_functional_assignment.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 11. Summarize by comparison × category × direction
# ======================================================

category_summary_raw <- sig_df %>%

  group_by(

    Comparison,

    Functional_category,

    Enriched_group
  ) %>%

  summarise(

    N_pathways =
      n(),

    Mean_abs_log2FC =
      mean(
        abs(
          Log2FC
        ),
        na.rm = TRUE
      ),

    .groups =
      "drop"
  )


# ======================================================
# 12. Functional-category ordering
#
# Categories with more differential pathways overall
# are placed toward the top of the horizontal figure.
# ======================================================

category_order <- category_summary_raw %>%

  group_by(
    Functional_category
  ) %>%

  summarise(

    Total =
      sum(
        N_pathways
      ),

    .groups =
      "drop"
  ) %>%

  arrange(
    Total
  ) %>%

  pull(
    Functional_category
  )


# Handle the unlikely case that no significant
# pathways are found in any comparison.
if (
  length(
    category_order
  ) == 0
) {

  category_order <- "No significant pathways"
}


# ======================================================
# 13. Create complete plotting grid
#
# Explicitly construct:
#
# Comparison × Functional category × Enriched group
#
# Missing combinations are filled with zero.
# ======================================================

plot_grid <- tidyr::expand_grid(

  Comparison =
    comparison_levels,

  Functional_category =
    category_order,

  Enriched_group =
    c(
      "ITM",
      "OTM1",
      "OTM2"
    )
)


category_summary <- plot_grid %>%

  left_join(

    category_summary_raw,

    by = c(
      "Comparison",
      "Functional_category",
      "Enriched_group"
    )
  ) %>%

  mutate(

    N_pathways =
      replace_na(
        N_pathways,
        0
      ),

    Comparison = factor(
      Comparison,
      levels =
        comparison_levels
    ),

    Functional_category = factor(
      Functional_category,
      levels =
        category_order
    ),

    Enriched_group = factor(
      Enriched_group,
      levels = c(
        "ITM",
        "OTM1",
        "OTM2"
      )
    )
  )


# ======================================================
# 14. Identify comparisons with no significant pathways
# ======================================================

comparison_totals <- category_summary %>%

  group_by(
    Comparison
  ) %>%

  summarise(

    Total_significant_pathways =
      sum(
        N_pathways
      ),

    .groups =
      "drop"
  )


empty_comparisons <- comparison_totals %>%

  filter(
    Total_significant_pathways == 0
  ) %>%

  pull(
    Comparison
  )


cat(
  "\nComparisons with no FDR-significant pathways:\n"
)


print(
  empty_comparisons
)


# ======================================================
# 15. Empty-panel annotation
# ======================================================

if (
  length(
    empty_comparisons
  ) > 0
) {


  middle_category <- category_order[
    ceiling(
      length(
        category_order
      ) /
        2
    )
  ]


  blank_label_df <- tibble(

    Comparison = factor(
      empty_comparisons,
      levels =
        comparison_levels
    ),

    Functional_category = factor(
      middle_category,
      levels =
        category_order
    ),

    x =
      0.5,

    label =
      "No FDR-significant pathways"
  )


} else {


  blank_label_df <- tibble(

    Comparison =
      factor(
        character(),
        levels =
          comparison_levels
      ),

    Functional_category =
      factor(
        character(),
        levels =
          category_order
      ),

    x =
      numeric(),

    label =
      character()
  )
}


# ======================================================
# 16. Export functional-category summary
# ======================================================

write.csv(

  category_summary_raw %>%

    arrange(
      Comparison,
      desc(
        N_pathways
      ),
      Functional_category
    ),

  file.path(
    output_dir,
    "stat_MetaCyc_functional_category_summary.csv"
  ),

  row.names = FALSE
)


write.csv(

  comparison_totals,

  file.path(
    output_dir,
    "stat_MetaCyc_significant_pathway_counts_by_comparison.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 17. Colors
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
# 18. Calculate appropriate fixed x-axis maximum
#
# Because bars are stacked by Enriched_group, the
# maximum must be based on the TOTAL bar length rather
# than on an individual colored segment.
# ======================================================

stacked_totals <- category_summary %>%

  group_by(
    Comparison,
    Functional_category
  ) %>%

  summarise(

    Total_N =
      sum(
        N_pathways
      ),

    .groups =
      "drop"
  )


max_stacked_n <- max(
  stacked_totals$Total_N,
  na.rm = TRUE
)


if (
  is.infinite(
    max_stacked_n
  ) ||
    is.na(
      max_stacked_n
    ) ||
    max_stacked_n < 1
) {

  max_stacked_n <- 1
}


# ======================================================
# 19. Functional-category bar plot
# ======================================================

p_category <- ggplot(

  category_summary,

  aes(

    x =
      N_pathways,

    y =
      Functional_category,

    fill =
      Enriched_group
  )

) +


  geom_col(

    width =
      0.75,

    color =
      "black",

    linewidth =
      0.2
  ) +


  geom_text(

    data =
      blank_label_df,

    aes(

      x =
        x,

      y =
        Functional_category,

      label =
        label
    ),

    inherit.aes =
      FALSE,

    size =
      4,

    fontface =
      "italic",

    color =
      "grey35"
  ) +


  facet_wrap(

    ~ Comparison,

    nrow =
      1,

    scales =
      if (
        axis_mode ==
          "fixed"
      ) {

        "fixed"

      } else {

        "free_x"
      }
  ) +


  scale_fill_manual(

    values =
      group_colors,

    drop =
      FALSE,

    name =
      NULL
  ) +


  scale_x_continuous(

    breaks =
      scales::pretty_breaks(),

    expand =
      expansion(
        mult = c(
          0,
          0.03
        )
      )
  ) +


  labs(

    x =
      "Number of differential pathways",

    y =
      "Functional category",

    title =
      "Functional categories of differential MetaCyc pathways",

    subtitle =
      paste0(
        "D30 + D60 + D90; Wilcoxon rank-sum test; ",
        "BH-adjusted FDR < 0.05"
      )
  ) +


  theme_bw(
    base_size = 12
  ) +


  theme(

    panel.grid.major.y =
      element_blank(),

    panel.grid.minor =
      element_blank(),


    axis.text =
      element_text(
        color = "black"
      ),


    axis.title =
      element_text(
        face = "bold"
      ),


    strip.background =
      element_rect(
        fill = "white",
        color = "black"
      ),


    strip.text =
      element_text(
        face = "bold",
        size = 12
      ),


    legend.position =
      "bottom",


    plot.title =
      element_text(
        face = "bold",
        hjust = 0.5
      ),


    plot.subtitle =
      element_text(
        hjust = 0.5
      )
  )


# ======================================================
# 20. Fixed x-axis if requested
# ======================================================

if (
  axis_mode ==
    "fixed"
) {


  p_category <- p_category +

    coord_cartesian(

      xlim = c(
        0,
        max_stacked_n *
          1.08
      ),

      clip =
        "off"
    )
}


# ======================================================
# 21. Display figure
# ======================================================

print(
  p_category
)


# ======================================================
# 22. Save figure
# ======================================================

ggsave(

  filename = file.path(
    output_dir,
    "Figure_MetaCyc_differential_pathway_functional_categories.pdf"
  ),

  plot =
    p_category,

  width =
    13,

  height =
    6
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_MetaCyc_differential_pathway_functional_categories.png"
  ),

  plot =
    p_category,

  width =
    13,

  height =
    6,

  dpi =
    600
)


# ======================================================
# 23. Print summary
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "Functional-category summary\n"
)

cat(
  "========================================\n"
)


print(

  category_summary_raw %>%

    arrange(
      Comparison,
      desc(
        N_pathways
      )
    )
)


# ======================================================
# 24. Finish
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "MetaCyc functional-category analysis completed.\n"
)

cat(
  "Only pathways with FDR < 0.05 were summarized.\n"
)

cat(
  "Functional categories were assigned by predefined ",
  "keyword matching against pathway descriptions.\n",
  sep = ""
)

cat(
  "Axis mode:",
  axis_mode,
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
  "Output directory:\n"
)

cat(
  output_dir,
  "\n"
)

cat(
  "========================================\n"
)