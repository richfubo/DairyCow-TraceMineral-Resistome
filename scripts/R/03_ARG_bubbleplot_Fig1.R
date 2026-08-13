########################################################
## ARG abundance bubble plot
##
## Bubble size:
## Median normalized abundance
##
## Bubble color:
## log10(normalized abundance + 1e-6)
##
## Metadata columns:
## ID | Time | Group | Sample
##
## Group:
## ITM, OTM1, OTM2
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
  "Fig1_ARG_bubbleplot"
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
  "Group",
  "Sample"
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


# Rename ID for downstream merging
meta_df <- meta_df %>%
  rename(
    SampleID = ID
  ) %>%
  mutate(
    SampleID = as.character(SampleID),
    Time = as.character(Time),
    Group = as.character(Group),
    Sample = as.character(Sample)
  )


# ======================================================
# 4. Set treatment and time order
# ======================================================

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
  "Metadata rows:",
  nrow(meta_df),
  "\n"
)

cat(
  "Unique sequencing samples:",
  length(
    unique(meta_df$SampleID)
  ),
  "\n"
)

cat(
  "Unique cows:",
  length(
    unique(meta_df$Sample)
  ),
  "\n"
)

print(
  table(
    meta_df$Group,
    meta_df$Time
  )
)


if (anyDuplicated(meta_df$SampleID) > 0) {

  stop(
    "Duplicated sequencing sample IDs were found."
  )
}


# ======================================================
# 6. Read normalized ARG abundance matrix
#
# First column:
# ARG category
#
# Remaining columns:
# a1-a120
# ======================================================

feature_df <- read.delim(
  feature_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ======================================================
# 7. Convert abundance matrix to long format
# ======================================================

feature_col <- colnames(feature_df)[1]

arg_long <- feature_df %>%

  rename(
    ARG_subtype = all_of(feature_col)
  ) %>%

  mutate(
    ARG_subtype = as.character(ARG_subtype)
  ) %>%

  pivot_longer(

    cols = -ARG_subtype,

    names_to = "SampleID",

    values_to = "Abundance"

  ) %>%

  mutate(

    SampleID = as.character(SampleID),

    Abundance = as.numeric(Abundance)
  )


# Replace missing abundance with zero
arg_long$Abundance[
  is.na(arg_long$Abundance)
] <- 0


# Negative abundance values are not allowed
arg_long$Abundance[
  arg_long$Abundance < 0
] <- 0


# ======================================================
# 8. Check sample matching
# ======================================================

matrix_samples <- unique(
  arg_long$SampleID
)

common_samples <- intersect(
  meta_df$SampleID,
  matrix_samples
)


cat(
  "Samples matched between metadata and ARG table:",
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
  meta_df$SampleID,
  matrix_samples
)


if (length(missing_samples) > 0) {

  cat(
    "Samples missing from ARG abundance table:\n"
  )

  print(
    missing_samples
  )
}


# ======================================================
# 9. Merge metadata
# ======================================================

arg_long <- arg_long %>%

  inner_join(

    meta_df %>%
      select(
        SampleID,
        Time,
        Group
      ),

    by = "SampleID"
  )


# ======================================================
# 10. Calculate median abundance
#
# Median abundance is calculated for each:
# Time × Group × ARG category
# ======================================================

df_sum <- arg_long %>%

  group_by(
    Time,
    Group,
    ARG_subtype
  ) %>%

  summarise(

    Abundance = median(
      Abundance,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%

  mutate(

    logAbundance = log10(
      Abundance + 1e-6
    ),

    Time_Group = paste(
      Time,
      Group,
      sep = "_"
    )
  )


# ======================================================
# 11. Set ARG category order
# ======================================================

df_sum$ARG_subtype <- factor(

  df_sum$ARG_subtype,

  levels = sort(
    unique(
      df_sum$ARG_subtype
    )
  )
)


# ======================================================
# 12. Set time and treatment order
# ======================================================

df_sum$Time <- factor(
  df_sum$Time,
  levels = c(
    "D0",
    "D30",
    "D60",
    "D90"
  )
)


df_sum$Group <- factor(
  df_sum$Group,
  levels = c(
    "ITM",
    "OTM1",
    "OTM2"
  )
)


x_levels <- df_sum %>%

  distinct(
    Time,
    Group,
    Time_Group
  ) %>%

  arrange(
    Time,
    Group
  ) %>%

  pull(
    Time_Group
  )


df_sum$Time_Group <- factor(
  df_sum$Time_Group,
  levels = x_levels
)


# ======================================================
# 13. Bubble plot 1
#
# x = Time × Group
# y = ARG category
# ======================================================

p1 <- ggplot(

  df_sum,

  aes(
    x = Time_Group,
    y = ARG_subtype
  )

) +

  geom_point(

    aes(
      size = Abundance,
      color = logAbundance
    ),

    alpha = 0.85
  ) +

  scale_size_area(
    max_size = 12
  ) +

  labs(

    x = "Time × Group",

    y = "ARG subtype",

    size = "Median abundance",

    color = "log10(abundance)"
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),

    panel.grid.major = element_line(
      linewidth = 0.2
    ),

    panel.grid.minor = element_blank()
  )


print(
  p1
)


# ======================================================
# 14. Bubble plot 2
#
# Faceted by sampling time
# x = treatment group
# ======================================================

p2 <- ggplot(

  df_sum,

  aes(
    x = Group,
    y = ARG_subtype
  )

) +

  geom_point(

    aes(
      size = Abundance,
      color = logAbundance
    ),

    alpha = 0.85
  ) +

  facet_grid(
    . ~ Time,
    scales = "free_x",
    space = "free_x"
  ) +

  scale_size_area(
    max_size = 12
  ) +

  scale_x_discrete(
    drop = FALSE
  ) +

  labs(

    x = "Group",

    y = "ARG subtype",

    size = "Median abundance",

    color = "log10(abundance)"
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(

    panel.grid.major = element_line(
      linewidth = 0.2
    ),

    panel.grid.minor = element_blank()
  )


print(
  p2
)


# ======================================================
# 15. Save figures
# ======================================================

ggsave(

  filename = file.path(
    output_dir,
    "ARG_bubbleplot_TimeGroup.pdf"
  ),

  plot = p1,

  width = 12,

  height = 8
)


ggsave(

  filename = file.path(
    output_dir,
    "ARG_bubbleplot_facetTime.pdf"
  ),

  plot = p2,

  width = 14,

  height = 7
)


# ======================================================
# 16. Export summarized data
# ======================================================

write.csv(

  df_sum,

  file.path(
    output_dir,
    "ARG_bubbleplot_median_abundance.csv"
  ),

  row.names = FALSE
)


# ======================================================
# 17. Finish
# ======================================================

cat(
  "\nARG abundance bubble plot completed.\n"
)

cat(
  "Results saved to:\n",
  output_dir,
  "\n"
)