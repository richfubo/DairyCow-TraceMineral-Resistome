########################################################
## 16S genus-level stacked bar plot
##
## Figure 4
##
## Input:
##   metadata/metadata.csv
##
##   processed_data/16S/
##     feature-table.tsv
##     taxonomy.tsv
##
## Analysis:
##   1. Assign ASVs/features to genus level
##   2. Remove unclassified/uninterpretable genera
##   3. Calculate relative abundance using the total
##      abundance of the complete feature table
##   4. Select the Top 20 interpretable genera based on
##      mean relative abundance across all samples
##   5. Calculate mean relative abundance for each
##      Time × Group × Genus combination
##
## Plot:
##   x-axis = ITM / OTM1 / OTM2
##   facets = D0 / D30 / D60 / D90
##   fill = Top 20 interpretable genera
##
## Important:
##   No "Others" category is added.
##   Unclassified genera are not displayed.
##   Therefore, stacked bars are not expected to sum
##   to 100%.
##
## Output:
##   results/Fig4_16S_genus_stacked_bar/
########################################################


# ======================================================
# 1. Load packages
# ======================================================

library(tidyverse)
library(scales)
library(grid)


# ======================================================
# 2. Project paths
#
# Run this script from the repository root:
# DairyCow-TraceMineral-Resistome/
# ======================================================

PROJECT_DIR <- "."


feature_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "16S",
  "feature-table.tsv"
)


taxonomy_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "16S",
  "taxonomy.tsv"
)


meta_path <- file.path(
  PROJECT_DIR,
  "metadata",
  "metadata.csv"
)


output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig4_16S_genus_stacked_bar"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ======================================================
# 3. Read metadata
# ======================================================

metadata <- read.csv(
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
  colnames(metadata)
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


metadata <- metadata %>%

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
  unique(metadata$Group),
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
  unique(metadata$Time),
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


if (anyDuplicated(metadata$ID) > 0) {

  stop(
    "Duplicated sequencing sample IDs were found in metadata.csv."
  )
}


metadata$Group <- factor(
  metadata$Group,
  levels = expected_groups
)


metadata$Time <- factor(
  metadata$Time,
  levels = expected_times
)


cat(
  "\n========================================\n"
)

cat(
  "16S genus stacked-bar analysis\n"
)

cat(
  "========================================\n"
)


cat(
  "\nMetadata samples:",
  nrow(metadata),
  "\n"
)


cat(
  "\nSamples by Group × Time:\n"
)


print(
  table(
    metadata$Group,
    metadata$Time
  )
)


# ======================================================
# 5. Read QIIME2 feature table
#
# QIIME2-exported feature tables may contain:
#
# # Constructed from biom file
#
# as the first line. If present, it is skipped.
# ======================================================

first_feature_line <- readLines(
  feature_path,
  n = 1,
  warn = FALSE
)


feature_skip <- ifelse(

  grepl(
    "^# Constructed from biom file",
    first_feature_line
  ),

  1,

  0
)


feature <- read.delim(

  feature_path,

  header = TRUE,

  sep = "\t",

  skip = feature_skip,

  check.names = FALSE,

  comment.char = "",

  quote = "",

  stringsAsFactors = FALSE
)


# First column = feature/ASV ID
colnames(feature)[1] <- "FeatureID"


feature$FeatureID <- as.character(
  feature$FeatureID
)


# ======================================================
# 6. Read taxonomy
# ======================================================

taxonomy <- read.delim(

  taxonomy_path,

  header = TRUE,

  sep = "\t",

  check.names = FALSE,

  comment.char = "",

  quote = "",

  stringsAsFactors = FALSE
)


# Compatible with:
# Feature ID
# Feature.ID
# #OTU ID
colnames(taxonomy) <- gsub(
  "^Feature\\.ID$",
  "FeatureID",
  colnames(taxonomy)
)


colnames(taxonomy) <- gsub(
  "^Feature ID$",
  "FeatureID",
  colnames(taxonomy)
)


colnames(taxonomy) <- gsub(
  "^#OTU ID$",
  "FeatureID",
  colnames(taxonomy)
)


if (!"FeatureID" %in% colnames(taxonomy)) {

  colnames(taxonomy)[1] <- "FeatureID"
}


if (!"Taxon" %in% colnames(taxonomy)) {

  stop(
    "taxonomy.tsv must contain a Taxon column."
  )
}


taxonomy$FeatureID <- as.character(
  taxonomy$FeatureID
)


# ======================================================
# 7. Check feature IDs
# ======================================================

cat(
  "\nFeature table rows:",
  nrow(feature),
  "\n"
)


cat(
  "Taxonomy rows:",
  nrow(taxonomy),
  "\n"
)


matched_features <- intersect(
  feature$FeatureID,
  taxonomy$FeatureID
)


cat(
  "Features with taxonomy:",
  length(matched_features),
  "\n"
)


if (length(matched_features) == 0) {

  stop(
    paste0(
      "No FeatureID overlap was found between ",
      "feature-table.tsv and taxonomy.tsv."
    )
  )
}


# ======================================================
# 8. Check sequencing sample columns
# ======================================================

feature_sample_cols <- setdiff(
  colnames(feature),
  "FeatureID"
)


common_samples <- intersect(
  feature_sample_cols,
  metadata$ID
)


cat(
  "\nFeature-table samples:",
  length(feature_sample_cols),
  "\n"
)


cat(
  "Matched metadata samples:",
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


missing_in_feature <- setdiff(
  metadata$ID,
  feature_sample_cols
)


if (length(missing_in_feature) > 0) {

  cat(
    "\nSamples in metadata but missing from feature table:\n"
  )

  print(
    missing_in_feature
  )
}


extra_feature_samples <- setdiff(
  feature_sample_cols,
  metadata$ID
)


if (length(extra_feature_samples) > 0) {

  cat(
    "\nSamples in feature table but missing from metadata:\n"
  )

  print(
    extra_feature_samples
  )
}


if (length(common_samples) < 5) {

  stop(
    paste0(
      "Too few matched samples. ",
      "Please check sample IDs in feature-table.tsv ",
      "and metadata.csv."
    )
  )
}


# Retain only matched sequencing samples
feature <- feature %>%

  select(
    FeatureID,
    all_of(common_samples)
  )


metadata <- metadata[
  match(
    common_samples,
    metadata$ID
  ),
  ,
  drop = FALSE
]


# ======================================================
# 9. Convert feature abundance columns to numeric
# ======================================================

feature[
  common_samples
] <- lapply(

  feature[
    common_samples
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
# 10. Parse taxonomy to genus level
# ======================================================

taxonomy <- taxonomy %>%

  separate(

    Taxon,

    into = c(
      "Kingdom",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species"
    ),

    sep = ";",

    fill = "right",

    extra = "merge",

    remove = FALSE
  ) %>%

  mutate(

    across(

      c(
        Kingdom,
        Phylum,
        Class,
        Order,
        Family,
        Genus,
        Species
      ),

      ~ trimws(.x)
    ),


    Genus = gsub(
      "^g__",
      "",
      Genus
    ),


    Genus = gsub(
      "^D_5__",
      "",
      Genus
    ),


    Genus = ifelse(

      is.na(Genus) |

        Genus == "" |

        Genus == "__" |

        str_detect(
          tolower(Genus),
          "uncultured|unclassified|unknown|metagenome"
        ),

      "Unclassified",

      Genus
    )
  )


# ======================================================
# 11. Clean genus names
#
# Remove selected suffixes used in the original analysis.
# ======================================================

clean_genus_name <- function(x) {


  x <- as.character(x)


  # Numeric suffix:
  # Bifidobacterium_38875 -> Bifidobacterium
  x <- str_replace(
    x,
    "_[0-9]+$",
    ""
  )


  # Selected A-E suffixes
  x <- str_replace(
    x,
    "_A$",
    ""
  )

  x <- str_replace(
    x,
    "_B$",
    ""
  )

  x <- str_replace(
    x,
    "_C$",
    ""
  )

  x <- str_replace(
    x,
    "_D$",
    ""
  )

  x <- str_replace(
    x,
    "_E$",
    ""
  )


  return(
    x
  )
}


taxonomy <- taxonomy %>%

  mutate(

    Genus_clean =
      clean_genus_name(
        Genus
      )
  )


# ======================================================
# 12. Define interpretable genera
#
# The following filtering rules are retained from the
# original analysis.
# ======================================================

bad_genus_pattern <- paste(

  c(

    "^CAG",
    "^UBA",
    "^UCG",
    "^RUG",
    "^SFM",
    "^PeH",
    "^RF",
    "^NK",
    "^GCA",
    "^GCF",
    "^DTU",
    "^QAM",
    "^QEM",
    "^QHM",
    "^QYM",
    "^QZM",
    "^QGM",

    "^Ambiguous_taxa",

    "^Incertae_Sedis",

    "^Family_XIII",

    "^Lachnospiraceae_[A-Z0-9]",

    "^Ruminococcaceae_[A-Z0-9]",

    "^Eubacterium_",

    "^Blautia_A_",
    "^Blautia_B_",
    "^Blautia_C_",

    "^Clostridium_[A-Z]",

    "^Bacteroides_[A-Z]",

    "^Prevotella_[A-Z]",

    "^Cryptobacteroides$",

    "^Paramuribaculum$"
  ),

  collapse = "|"
)


taxonomy <- taxonomy %>%

  mutate(

    Is_unclassified =
      Genus_clean ==
      "Unclassified",


    Is_interpretable = case_when(

      Is_unclassified ~
        FALSE,

      str_detect(
        Genus,
        bad_genus_pattern
      ) ~
        FALSE,

      TRUE ~
        TRUE
    )
  )


cat(
  "\nTaxonomic genera before filtering:\n"
)


print(

  taxonomy %>%

    count(
      Genus_clean,
      sort = TRUE
    ) %>%

    slice_head(
      n = 30
    )
)


cat(
  "\nInterpretable genus candidates:\n"
)


print(

  taxonomy %>%

    filter(
      Is_interpretable
    ) %>%

    count(
      Genus_clean,
      sort = TRUE
    ) %>%

    slice_head(
      n = 30
    )
)


# ======================================================
# 13. Export genus classification/filtering table
# ======================================================

write.csv(

  taxonomy %>%

    select(
      FeatureID,
      Taxon,
      Genus,
      Genus_clean,
      Is_unclassified,
      Is_interpretable
    ),

  file.path(
    output_dir,
    "stat_16S_genus_taxonomy_filtering.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 14. Merge feature table with taxonomy
# ======================================================

feature_tax <- feature %>%

  left_join(

    taxonomy %>%

      select(
        FeatureID,
        Genus,
        Genus_clean,
        Is_interpretable,
        Is_unclassified
      ),

    by =
      "FeatureID"
  ) %>%

  mutate(

    Genus = ifelse(
      is.na(Genus),
      "Unclassified",
      Genus
    ),

    Genus_clean = ifelse(
      is.na(Genus_clean),
      "Unclassified",
      Genus_clean
    ),

    Is_unclassified = ifelse(
      is.na(Is_unclassified),
      TRUE,
      Is_unclassified
    ),

    Is_interpretable = ifelse(
      is.na(Is_interpretable),
      FALSE,
      Is_interpretable
    )
  )


# ======================================================
# 15. Calculate total abundance using ALL features
#
# Relative abundance denominator is based on the entire
# feature table, including genera that are subsequently
# excluded from display.
# ======================================================

sample_total <- feature_tax %>%

  summarise(

    across(
      all_of(common_samples),
      ~ sum(
        .x,
        na.rm = TRUE
      )
    )
  ) %>%

  pivot_longer(

    cols =
      everything(),

    names_to =
      "ID",

    values_to =
      "Total_abundance"
  )


if (
  any(
    sample_total$Total_abundance <= 0
  )
) {

  stop(
    "One or more samples have zero total feature abundance."
  )
}


# ======================================================
# 16. Aggregate all interpretable genera
# ======================================================

genus_interpretable <- feature_tax %>%

  filter(

    Is_interpretable,

    Genus_clean !=
      "Unclassified"
  ) %>%

  group_by(
    Genus_clean
  ) %>%

  summarise(

    across(

      all_of(
        common_samples
      ),

      ~ sum(
        .x,
        na.rm = TRUE
      )
    ),

    .groups =
      "drop"
  ) %>%

  rename(
    Genus =
      Genus_clean
  )


cat(
  "\nInterpretable genus-level abundance table:\n"
)


cat(
  "Genera:",
  nrow(
    genus_interpretable
  ),
  "\n"
)


cat(
  "Samples:",
  length(
    common_samples
  ),
  "\n"
)


if (
  nrow(
    genus_interpretable
  ) == 0
) {

  stop(
    paste0(
      "No interpretable genera were retained. ",
      "Please check taxonomy and filtering rules."
    )
  )
}


# ======================================================
# 17. Convert to relative abundance
# ======================================================

df_interpretable <- genus_interpretable %>%

  pivot_longer(

    cols =
      all_of(
        common_samples
      ),

    names_to =
      "ID",

    values_to =
      "Abundance"
  ) %>%

  mutate(
    ID =
      as.character(
        ID
      )
  ) %>%

  left_join(
    sample_total,
    by = "ID"
  ) %>%

  left_join(

    metadata %>%

      select(
        ID,
        Time,
        Group
      ),

    by =
      "ID"
  ) %>%

  mutate(

    RelAbund =
      Abundance /
      Total_abundance
  )


df_interpretable$RelAbund[
  is.na(
    df_interpretable$RelAbund
  )
] <- 0


# ======================================================
# 18. Select Top 20 interpretable genera
#
# Ranking is based on mean whole-community relative
# abundance across all matched samples.
# ======================================================

top20_table <- df_interpretable %>%

  filter(
    !is.na(Time),
    !is.na(Group)
  ) %>%

  group_by(
    Genus
  ) %>%

  summarise(

    mean_relative_abundance =
      mean(
        RelAbund,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  ) %>%

  arrange(
    desc(
      mean_relative_abundance
    )
  ) %>%

  slice_head(
    n = 20
  ) %>%

  mutate(
    Rank =
      row_number()
  )


top20_interpretable <- top20_table$Genus


cat(
  "\nTop 20 interpretable genera:\n"
)


print(
  top20_table
)


# Export Top20 list
write.csv(

  top20_table,

  file.path(
    output_dir,
    "stat_16S_Top20_interpretable_genera.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 19. Retain Top20 genera
# ======================================================

df_top20 <- df_interpretable %>%

  filter(
    Genus %in%
      top20_interpretable
  ) %>%

  mutate(

    Genus = factor(
      Genus,
      levels =
        top20_interpretable
    )
  )


# ======================================================
# 20. Summarize Time × Group × Genus
#
# Mean relative abundance across biological samples.
# ======================================================

df_sum <- df_top20 %>%

  filter(
    !is.na(Time),
    !is.na(Group)
  ) %>%

  group_by(
    Time,
    Group,
    Genus
  ) %>%

  summarise(

    RelAbund =
      mean(
        RelAbund,
        na.rm = TRUE
      ),

    n_samples =
      n(),

    .groups =
      "drop"
  )


# ======================================================
# 21. Check Top20 proportion represented in each bar
#
# Values are intentionally allowed to be < 1 because
# non-Top20 and unclassified genera are not displayed.
# ======================================================

bar_totals <- df_sum %>%

  group_by(
    Time,
    Group
  ) %>%

  summarise(

    total_top20 =
      sum(
        RelAbund,
        na.rm = TRUE
      ),

    .groups =
      "drop"
  )


cat(
  "\nTop20 interpretable genera total relative abundance per bar:\n"
)


print(
  bar_totals
)


# ======================================================
# 22. Export plotting data
# ======================================================

write.csv(

  df_sum,

  file.path(
    output_dir,
    "plot_16S_Top20_genus_mean_relative_abundance.csv"
  ),

  row.names = FALSE
)


write.csv(

  bar_totals,

  file.path(
    output_dir,
    "stat_16S_Top20_relative_abundance_per_bar.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 23. Colors
# ======================================================

top20_colors <- c(

  "#F8766D",
  "#FC8D62",
  "#FDB863",
  "#E6F598",
  "#ABDDA4",

  "#66C2A5",
  "#3288BD",
  "#5E4FA2",
  "#B2ABD2",
  "#F1A340",

  "#FF7F00",
  "#E31A1C",
  "#FB9A99",
  "#CAB2D6",
  "#6A3D9A",

  "#B3DE69",
  "#FFFF99",
  "#A6CEE3",
  "#1F78B4",
  "#33A02C"
)


names(
  top20_colors
) <- top20_interpretable


# ======================================================
# 24. Plot
# ======================================================

p <- ggplot(

  df_sum,

  aes(

    x =
      Group,

    y =
      RelAbund,

    fill =
      Genus
  )

) +


  geom_col(

    width =
      0.88,

    color =
      "grey30",

    linewidth =
      0.08
  ) +


  scale_fill_manual(

    values =
      top20_colors,

    drop =
      FALSE
  ) +


  scale_y_continuous(

    labels =
      scales::percent_format(
        accuracy = 1
      ),

    expand =
      c(
        0,
        0
      )
  ) +


  facet_wrap(

    ~ Time,

    nrow =
      1,

    scales =
      "free_x"
  ) +


  labs(

    x =
      "Group",

    y =
      "Relative abundance (%)",

    fill =
      "Genus"
  ) +


  guides(

    fill =
      guide_legend(

        ncol =
          1,

        byrow =
          TRUE
      )
  ) +


  theme_bw(
    base_size = 14
  ) +


  theme(

    panel.grid =
      element_blank(),


    strip.background =
      element_rect(
        fill = "grey90",
        color = "black"
      ),


    strip.text =
      element_text(
        face = "bold"
      ),


    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        color = "black"
      ),


    axis.text.y =
      element_text(
        color = "black"
      ),


    axis.title =
      element_text(
        face = "bold"
      ),


    panel.spacing.x =
      grid::unit(
        0.1,
        "lines"
      ),


    legend.title =
      element_text(
        face = "bold"
      ),


    legend.text =
      element_text(
        size = 9
      ),


    legend.position =
      "right"
  )


print(
  p
)


# ======================================================
# 25. Save figure
# ======================================================

ggsave(

  filename = file.path(
    output_dir,
    "Figure_16S_Top20_interpretable_genus_stacked_bar.pdf"
  ),

  plot =
    p,

  width =
    10,

  height =
    6
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_16S_Top20_interpretable_genus_stacked_bar.png"
  ),

  plot =
    p,

  width =
    10,

  height =
    6,

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
  "16S genus stacked-bar analysis completed.\n"
)

cat(
  "Top 20 interpretable genera only.\n"
)

cat(
  "Unclassified genera were excluded.\n"
)

cat(
  "No Others category was added.\n"
)

cat(
  "Therefore, bar heights are not expected to reach 100%.\n"
)

cat(
  "Feature-table input:\n"
)

cat(
  feature_path,
  "\n"
)

cat(
  "Taxonomy input:\n"
)

cat(
  taxonomy_path,
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