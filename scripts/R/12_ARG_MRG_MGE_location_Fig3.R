########################################################
## ARG / MRG / MGE genomic-location distribution
##
## Figure 3
##
## Panels:
##   C: ARG
##   D: MRG
##   E: MGE
##
## Input:
##   metadata/metadata.csv
##   processed_data/genomic_context/
##     ARG_MRG_MGE_sample_location_summary_wide_a31_a120.tsv
##
## Time points included:
##   D30 + D60 + D90
##
## Analysis:
##   Mean proportions of genes located on
##   plasmid-associated or chromosome-like/non-plasmid
##   contigs are calculated for each treatment group.
##
## Bars:
##   Stacked percentages normalized to 100%
##
## Labels:
##   Percentage labels are displayed within the
##   corresponding stacked regions.
##
## Output:
##   results/Fig3_genomic_location/
########################################################


# ======================================================
# 1. Load packages
# ======================================================

library(tidyverse)
library(cowplot)
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


location_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "genomic_context",
  "ARG_MRG_MGE_sample_location_summary_wide_a31_a120.tsv"
)


output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig3_genomic_location"
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


required_meta_cols <- c(
  "ID",
  "Group",
  "Time"
)


missing_meta_cols <- setdiff(
  required_meta_cols,
  colnames(meta_df)
)


if (length(missing_meta_cols) > 0) {

  stop(
    paste0(
      "Missing metadata columns: ",
      paste(
        missing_meta_cols,
        collapse = ", "
      )
    )
  )
}


# ======================================================
# 4. Read genomic-location summary
# ======================================================

loc_df <- read.delim(
  location_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


required_location_cols <- c(
  "sample",
  "category",
  "Plasmid_gene_proportion",
  "Chromosome_gene_proportion"
)


missing_location_cols <- setdiff(
  required_location_cols,
  colnames(loc_df)
)


if (length(missing_location_cols) > 0) {

  stop(
    paste0(
      "Missing genomic-location columns: ",
      paste(
        missing_location_cols,
        collapse = ", "
      )
    )
  )
}


# ======================================================
# 5. Prepare metadata
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


if (anyDuplicated(meta_df$ID) > 0) {

  stop(
    "Duplicated sequencing sample IDs were found."
  )
}


# ======================================================
# 6. Prepare genomic-location data
# ======================================================

loc_df <- loc_df %>%

  mutate(

    sample = as.character(sample),

    category = as.character(category),

    Plasmid_gene_proportion =
      as.numeric(
        Plasmid_gene_proportion
      ),

    Chromosome_gene_proportion =
      as.numeric(
        Chromosome_gene_proportion
      )
  )


# Check expected categories
expected_categories <- c(
  "ARG",
  "MRG",
  "MGE"
)


unexpected_categories <- setdiff(
  unique(loc_df$category),
  expected_categories
)


if (length(unexpected_categories) > 0) {

  warning(
    paste0(
      "Unexpected category value(s): ",
      paste(
        unexpected_categories,
        collapse = ", "
      )
    )
  )
}


loc_df$category <- factor(
  loc_df$category,
  levels = c(
    "ARG",
    "MRG",
    "MGE"
  )
)


# ======================================================
# 7. Match metadata and retain D30-D90
# ======================================================

df <- loc_df %>%

  left_join(

    meta_df %>%
      select(
        ID,
        Group,
        Time
      ),

    by = c(
      "sample" = "ID"
    )
  ) %>%

  filter(

    !is.na(Group),

    !is.na(Time),

    Time %in% c(
      "D30",
      "D60",
      "D90"
    )
  ) %>%

  droplevels()


# ======================================================
# 8. Sample checks
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "Genomic-location sample check\n"
)

cat(
  "========================================\n"
)


cat(
  "\nSamples by category × Group × Time:\n"
)


print(
  table(
    df$category,
    df$Group,
    df$Time
  )
)


cat(
  "\nUnique samples used:",
  length(
    unique(
      df$sample
    )
  ),
  "\n"
)


# ======================================================
# 9. Location colors
# ======================================================

location_cols <- c(

  Plasmid =
    "#F4B183",

  Chromosome =
    "#9DC3E6"
)


# ======================================================
# 10. Calculate mean genomic-location proportions
#
# D30, D60, and D90 are pooled within each treatment
# group before calculating the mean proportions.
# ======================================================

summary_df <- df %>%

  group_by(
    category,
    Group
  ) %>%

  summarise(

    Plasmid_percent =
      mean(
        Plasmid_gene_proportion,
        na.rm = TRUE
      ) *
      100,

    Chromosome_percent =
      mean(
        Chromosome_gene_proportion,
        na.rm = TRUE
      ) *
      100,

    n_samples =
      n(),

    .groups =
      "drop"
  ) %>%

  mutate(

    Total_percent =
      Plasmid_percent +
      Chromosome_percent,


    Plasmid_percent_norm =
      ifelse(

        Total_percent > 0,

        Plasmid_percent /
          Total_percent *
          100,

        0
      ),


    Chromosome_percent_norm =
      ifelse(

        Total_percent > 0,

        Chromosome_percent /
          Total_percent *
          100,

        0
      )
  )


# ======================================================
# 11. Convert to long format for plotting
#
# Chromosome is placed first and Plasmid second so that
# chromosome-associated proportions appear at the bottom
# of the stacked bars.
# ======================================================

plot_df <- summary_df %>%

  select(

    category,

    Group,

    n_samples,

    Plasmid_percent_norm,

    Chromosome_percent_norm
  ) %>%

  pivot_longer(

    cols = c(
      Chromosome_percent_norm,
      Plasmid_percent_norm
    ),

    names_to =
      "Location",

    values_to =
      "Percent"
  ) %>%

  mutate(

    Location = recode(

      Location,

      "Plasmid_percent_norm" =
        "Plasmid",

      "Chromosome_percent_norm" =
        "Chromosome"
    ),


    Location = factor(

      Location,

      levels = c(
        "Chromosome",
        "Plasmid"
      )
    ),


    category = factor(

      category,

      levels = c(
        "ARG",
        "MRG",
        "MGE"
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


    label_percent = case_when(

      Percent == 0 ~
        "",

      Percent < 3 ~
        "",

      TRUE ~
        sprintf(
          "%.1f%%",
          Percent
        )
    )
  )


cat(
  "\nPlot data:\n"
)


print(
  plot_df
)


# ======================================================
# 12. Export genomic-location summary
# ======================================================

write.csv(

  plot_df,

  file.path(
    output_dir,
    "ARG_MRG_MGE_location_percentage_summary.csv"
  ),

  row.names = FALSE
)


# Also export the wide-format summary
write.csv(

  summary_df,

  file.path(
    output_dir,
    "ARG_MRG_MGE_location_percentage_summary_wide.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 13. Plot theme
# ======================================================

theme_single <- theme_classic(
  base_size = 12
) +

  theme(

    legend.position =
      "top",

    legend.title =
      element_blank(),

    legend.text =
      element_text(
        color = "black",
        size = 8.5
      ),

    legend.key.size =
      grid::unit(
        0.35,
        "cm"
      ),

    legend.margin =
      margin(
        0,
        0,
        0,
        0
      ),

    legend.box.margin =
      margin(
        0,
        0,
        -3,
        0
      ),


    axis.title.x =
      element_blank(),

    axis.title.y =
      element_text(
        face = "bold",
        size = 10.5
      ),


    axis.text =
      element_text(
        color = "black",
        size = 8.8
      ),

    axis.text.x =
      element_text(
        color = "black",
        size = 9,
        angle = 0,
        hjust = 0.5
      ),

    axis.text.y =
      element_text(
        color = "black",
        size = 8.5
      ),


    axis.line =
      element_line(
        linewidth = 0.50
      ),

    axis.ticks =
      element_line(
        linewidth = 0.50
      ),

    axis.ticks.length =
      grid::unit(
        0.10,
        "cm"
      ),


    panel.border =
      element_rect(
        fill = NA,
        color = "black",
        linewidth = 0.55
      ),


    plot.margin =
      ggplot2::margin(
        4,
        4,
        4,
        2
      )
  )


# ======================================================
# 14. Stacked genomic-location bar plot
# ======================================================

plot_bar_only <- function(
    plot_df,
    cat_name
) {


  df_sub <- plot_df %>%

    filter(
      category == cat_name
    )


  p <- ggplot(

    df_sub,

    aes(
      x = Group,
      y = Percent,
      fill = Location
    )

  ) +


    geom_col(

      width = 0.62,

      color = "black",

      linewidth = 0.30,

      position = "stack"
    ) +


    geom_text(

      aes(
        label = label_percent
      ),

      position =
        position_stack(
          vjust = 0.5
        ),

      color = "black",

      size = 3.1,

      fontface = "bold"
    ) +


    scale_fill_manual(

      values =
        location_cols,

      breaks = c(
        "Plasmid",
        "Chromosome"
      ),

      labels = c(
        "Plasmid",
        "Chromosome"
      )
    ) +


    scale_y_continuous(

      name =
        "Proportion (%)",

      limits = c(
        0,
        100
      ),

      breaks = c(
        0,
        25,
        50,
        75,
        100
      ),

      expand =
        expansion(
          mult = c(
            0,
            0
          )
        )
    ) +


    labs(
      x = NULL
    ) +


    theme_single


  return(
    p
  )
}


# ======================================================
# 15. Vertical category strip
# ======================================================

make_vertical_strip <- function(
    label
) {


  ggplot() +


    annotate(

      "rect",

      xmin = 0,

      xmax = 1,

      ymin = 0,

      ymax = 1,

      fill = "grey90",

      color = NA
    ) +


    annotate(

      "text",

      x = 0.5,

      y = 0.5,

      label = label,

      angle = 90,

      fontface = "bold",

      size = 4
    ) +


    theme_void() +


    coord_cartesian(

      xlim = c(
        0,
        1
      ),

      ylim = c(
        0,
        1
      ),

      expand = FALSE
    )
}


# ======================================================
# 16. Construct individual panel
# ======================================================

make_panel <- function(
    plot_df,
    cat_name,
    panel_label
) {


  strip <- make_vertical_strip(
    cat_name
  )


  bar <- plot_bar_only(
    plot_df,
    cat_name
  )


  panel_core <- cowplot::plot_grid(

    strip,

    bar,

    nrow = 1,

    rel_widths = c(
      0.13,
      1
    )
  )


  panel_final <- cowplot::ggdraw(
    panel_core
  ) +

    cowplot::draw_label(

      panel_label,

      x = 0.00,

      y = 1.00,

      hjust = 0,

      vjust = 1,

      fontface = "bold",

      size = 14
    )


  return(
    panel_final
  )
}


# ======================================================
# 17. Generate panels C-E
# ======================================================

p_ARG <- make_panel(

  plot_df =
    plot_df,

  cat_name =
    "ARG",

  panel_label =
    "C"
)


p_MRG <- make_panel(

  plot_df =
    plot_df,

  cat_name =
    "MRG",

  panel_label =
    "D"
)


p_MGE <- make_panel(

  plot_df =
    plot_df,

  cat_name =
    "MGE",

  panel_label =
    "E"
)


print(
  p_ARG
)


print(
  p_MRG
)


print(
  p_MGE
)


# ======================================================
# 18. Combine panels C-E
# ======================================================

p_combined <- cowplot::plot_grid(

  p_ARG,

  p_MRG,

  p_MGE,

  nrow = 1,

  align = "h",

  rel_widths = c(
    1,
    1,
    1
  )
)


print(
  p_combined
)


# ======================================================
# 19. Save combined figure
# ======================================================

ggsave(

  filename = file.path(
    output_dir,
    "Figure_ARG_MRG_MGE_location_percentage_CDE.pdf"
  ),

  plot =
    p_combined,

  width =
    12,

  height =
    4.2
)


ggsave(

  filename = file.path(
    output_dir,
    "Figure_ARG_MRG_MGE_location_percentage_CDE.png"
  ),

  plot =
    p_combined,

  width =
    12,

  height =
    4.2,

  dpi =
    600
)


# ======================================================
# 20. Finish
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "ARG/MRG/MGE genomic-location analysis completed.\n"
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