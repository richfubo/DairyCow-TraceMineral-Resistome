########################################################
## D90 MRG Random Forest biomarker analysis
##
## Figure 2
##
## Input:
##   metadata/metadata.csv
##   processed_data/MRG/
##     MRG_merged_subtype_abundance.tsv
##
## Analysis:
##   D90 samples only
##
## Pairwise comparisons:
##   ITM vs OTM1
##   ITM vs OTM2
##   OTM1 vs OTM2
##
## Metal categories:
##   Copper
##   Zinc
##   Iron
##
## Random Forest:
##   ntree = 1000
##   importance = TRUE
##   top_n = 6 genes per metal category
##
## Feature importance:
##   MeanDecreaseAccuracy is used when available;
##   otherwise MeanDecreaseGini is used.
##
## Point color:
##   Treatment group with higher mean abundance
##
## Output:
##   results/Fig2_MRG_random_forest/
########################################################


# ======================================================
# 1. Load packages
# ======================================================

library(tidyverse)
library(randomForest)
library(stringr)
library(patchwork)
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


mrg_path <- file.path(
  PROJECT_DIR,
  "processed_data",
  "MRG",
  "MRG_merged_subtype_abundance.tsv"
)


output_dir <- file.path(
  PROJECT_DIR,
  "results",
  "Fig2_MRG_random_forest"
)


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Number of top genes retained for each metal category
top_n <- 6


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
# 4. Read normalized MRG subtype matrix
# ======================================================

mrg_df <- read.delim(
  mrg_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


colnames(mrg_df)[1] <- "feature"


mrg_df$feature <- as.character(
  mrg_df$feature
)


# ======================================================
# 5. Prepare metadata
#
# Expected:
# ID | Time | Group | Sample
#
# Group:
# ITM, OTM1, OTM2
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


if (anyDuplicated(meta_df$ID) > 0) {

  stop(
    "Duplicated sequencing sample IDs were found."
  )
}


# ======================================================
# 6. Retain D90 samples
# ======================================================

meta_d90 <- meta_df %>%

  filter(
    Time == "D90"
  ) %>%

  droplevels()


cat(
  "\n========================================\n"
)

cat(
  "D90 Random Forest analysis\n"
)

cat(
  "========================================\n"
)


cat(
  "\nD90 sample size per group:\n"
)


print(
  table(
    meta_d90$Group
  )
)


cat(
  "\nTotal D90 samples:",
  nrow(meta_d90),
  "\n"
)


# ======================================================
# 7. Metal-related feature keywords
# ======================================================

metal_keywords <- list(


  Copper =
    "copper|\\bcop|\\bpco|\\bcus|\\bcue|\\btcrb|\\bsil",


  Zinc =
    "zinc|\\bznt|\\bzit|\\bznu|\\bczc|\\bzur",


  Iron =
    paste0(
      "iron|\\bfeo|\\bfep|\\bfhu|\\bfie|",
      "\\btonb|\\bfur|\\bdps|ferritin|",
      "bacterioferritin|\\bbfr"
    )
)


metal_order <- c(
  "Copper",
  "Zinc",
  "Iron"
)


# ======================================================
# 8. Gene-name cleaner
#
# Examples:
# Copper_dnaK -> dnaK
# Zinc_zntA   -> zntA
# Iron_fur    -> fur
# ======================================================

clean_gene_name <- function(x) {


  x0 <- as.character(
    x
  )


  x0 <- basename(
    x0
  )


  # Standardize delimiters
  x1 <- x0 %>%

    str_replace_all(
      "\\|",
      "_"
    ) %>%

    str_replace_all(
      ":",
      "_"
    ) %>%

    str_replace_all(
      ";",
      "_"
    ) %>%

    str_replace_all(
      "\\s+",
      "_"
    )


  # Remove metal prefixes
  x1 <- str_replace(

    x1,

    regex(
      "^(Copper|Cooper|Cu|Zinc|Zn|Iron|Fe)[_\\-\\.]+",
      ignore_case = TRUE
    ),

    ""
  )


  x_low <- tolower(
    x1
  )


  gene <- case_when(


    # --------------------------
    # Copper-related genes
    # --------------------------

    str_detect(x_low, "tcrb") ~ "tcrB",


    str_detect(x_low, "copa") ~ "copA",
    str_detect(x_low, "copb") ~ "copB",
    str_detect(x_low, "copc") ~ "copC",
    str_detect(x_low, "copd") ~ "copD",
    str_detect(x_low, "cope") ~ "copE",
    str_detect(x_low, "copr") ~ "copR",
    str_detect(x_low, "cops") ~ "copS",


    str_detect(x_low, "pcoa") ~ "pcoA",
    str_detect(x_low, "pcob") ~ "pcoB",
    str_detect(x_low, "pcoc") ~ "pcoC",
    str_detect(x_low, "pcod") ~ "pcoD",
    str_detect(x_low, "pcoe") ~ "pcoE",
    str_detect(x_low, "pcor") ~ "pcoR",
    str_detect(x_low, "pcos") ~ "pcoS",


    str_detect(x_low, "cusa") ~ "cusA",
    str_detect(x_low, "cusb") ~ "cusB",
    str_detect(x_low, "cusc") ~ "cusC",
    str_detect(x_low, "cusf") ~ "cusF",
    str_detect(x_low, "cusr") ~ "cusR",
    str_detect(x_low, "cuss") ~ "cusS",


    str_detect(x_low, "cueo") ~ "cueO",
    str_detect(x_low, "cuer") ~ "cueR",


    str_detect(x_low, "sila") ~ "silA",
    str_detect(x_low, "silb") ~ "silB",
    str_detect(x_low, "silc") ~ "silC",
    str_detect(x_low, "sile") ~ "silE",
    str_detect(x_low, "silf") ~ "silF",
    str_detect(x_low, "silp") ~ "silP",
    str_detect(x_low, "silr") ~ "silR",
    str_detect(x_low, "sils") ~ "silS",


    # --------------------------
    # Zinc-related genes
    # --------------------------

    str_detect(x_low, "znta") ~ "zntA",
    str_detect(x_low, "zntr") ~ "zntR",
    str_detect(x_low, "zitb") ~ "zitB",


    str_detect(x_low, "znua") ~ "znuA",
    str_detect(x_low, "znub") ~ "znuB",
    str_detect(x_low, "znuc") ~ "znuC",


    str_detect(x_low, "czca") ~ "czcA",
    str_detect(x_low, "czcb") ~ "czcB",
    str_detect(x_low, "czcc") ~ "czcC",
    str_detect(x_low, "czcd") ~ "czcD",
    str_detect(x_low, "czcr") ~ "czcR",
    str_detect(x_low, "czcs") ~ "czcS",


    str_detect(x_low, "zur") ~ "zur",


    # --------------------------
    # Iron-related genes
    # --------------------------

    str_detect(x_low, "acn") ~ "acn",
    str_detect(x_low, "pmrg") ~ "pmrG",

    str_detect(x_low, "ybto") ~ "ybtO",
    str_detect(x_low, "ybtp") ~ "ybtP",

    str_detect(x_low, "feta") ~ "fetA",
    str_detect(x_low, "fetb") ~ "fetB",


    str_detect(x_low, "feoa") ~ "feoA",
    str_detect(x_low, "feob") ~ "feoB",
    str_detect(x_low, "feoc") ~ "feoC",


    str_detect(x_low, "fief") ~ "fieF",


    str_detect(x_low, "fepa") ~ "fepA",
    str_detect(x_low, "fepb") ~ "fepB",
    str_detect(x_low, "fepc") ~ "fepC",
    str_detect(x_low, "fepd") ~ "fepD",
    str_detect(x_low, "fepg") ~ "fepG",


    str_detect(x_low, "fhua") ~ "fhuA",
    str_detect(x_low, "fhub") ~ "fhuB",
    str_detect(x_low, "fhuc") ~ "fhuC",
    str_detect(x_low, "fhud") ~ "fhuD",


    str_detect(x_low, "tonb") ~ "tonB",
    str_detect(x_low, "fur") ~ "fur",
    str_detect(x_low, "dps") ~ "dps",


    str_detect(
      x_low,
      "ferritin"
    ) ~
      "ferritin",


    str_detect(
      x_low,
      "bacterioferritin|\\bbfr\\b|bfr"
    ) ~
      "bfr",


    TRUE ~
      NA_character_
  )


  # If no known gene family is detected,
  # retain the final non-empty field.
  fallback <- sapply(

    str_split(
      x1,
      "_"
    ),

    function(v) {

      v <- v[
        v != ""
      ]


      if (length(v) == 0) {

        return(
          NA_character_
        )

      } else {

        return(
          v[
            length(v)
          ]
        )
      }
    }
  )


  gene <- ifelse(

    is.na(gene) |
      gene == "",

    fallback,

    gene
  )


  gene <- ifelse(

    is.na(gene) |
      gene == "",

    x1,

    gene
  )


  gene <- str_replace_all(
    gene,
    "\\s+",
    "_"
  )


  return(
    gene
  )
}


# ======================================================
# 9. Build D90 abundance matrix for one metal category
# ======================================================

build_metal_matrix_D90 <- function(
    metal_name,
    keyword_pattern
) {


  metal_df <- mrg_df %>%

    filter(

      str_detect(
        feature,
        regex(
          keyword_pattern,
          ignore_case = TRUE
        )
      )
    )


  if (nrow(metal_df) == 0) {

    warning(
      paste0(
        "No ",
        metal_name,
        "-related subtype was found."
      )
    )

    return(
      NULL
    )
  }


  x <- metal_df[
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
    metal_df
  )[-1]


  colnames(mat) <- make.unique(
    as.character(
      metal_df$feature
    )
  )


  mat[is.na(mat)] <- 0

  mat[mat < 0] <- 0


  # Match D90 samples
  common_ids <- intersect(
    meta_d90$ID,
    rownames(mat)
  )


  if (length(common_ids) == 0) {

    warning(
      paste0(
        metal_name,
        ": no matched D90 samples."
      )
    )

    return(
      NULL
    )
  }


  mat <- mat[
    common_ids,
    ,
    drop = FALSE
  ]


  meta_sub <- meta_d90[
    match(
      common_ids,
      meta_d90$ID
    ),
    ,
    drop = FALSE
  ]


  stopifnot(
    all(
      rownames(mat) ==
        meta_sub$ID
    )
  )


  # Remove features with zero total abundance
  # or zero variance
  valid_cols <- apply(

    mat,

    2,

    function(v) {

      sum(
        v,
        na.rm = TRUE
      ) > 0 &&

        var(
          v,
          na.rm = TRUE
        ) > 0
    }
  )


  mat <- mat[
    ,
    valid_cols,
    drop = FALSE
  ]


  if (ncol(mat) == 0) {

    warning(
      paste0(
        metal_name,
        ": no variable subtype available for Random Forest."
      )
    )

    return(
      NULL
    )
  }


  return(

    list(

      mat = mat,

      meta = meta_sub,

      metal = metal_name
    )
  )
}


# ======================================================
# 10. Run Random Forest for one metal and one pair
# ======================================================

run_rf_one_metal_pair <- function(
    metal_name,
    keyword_pattern,
    g1,
    g2,
    top_n = 6
) {


  metal_obj <- build_metal_matrix_D90(

    metal_name =
      metal_name,

    keyword_pattern =
      keyword_pattern
  )


  if (is.null(metal_obj)) {

    return(
      NULL
    )
  }


  mat <- metal_obj$mat

  meta <- metal_obj$meta


  # Select the two treatment groups
  sub_meta <- meta %>%

    filter(
      Group %in% c(
        g1,
        g2
      )
    ) %>%

    droplevels()


  sub_mat <- mat[
    sub_meta$ID,
    ,
    drop = FALSE
  ]


  # Remove features without variation
  valid_cols <- apply(

    sub_mat,

    2,

    function(v) {

      var(
        v,
        na.rm = TRUE
      ) > 0
    }
  )


  sub_mat <- sub_mat[
    ,
    valid_cols,
    drop = FALSE
  ]


  if (ncol(sub_mat) == 0) {

    warning(
      paste0(
        metal_name,
        " | ",
        g1,
        " vs ",
        g2,
        ": no variable subtype."
      )
    )

    return(
      NULL
    )
  }


  sub_mat <- as.data.frame(
    sub_mat
  )


  sub_mat[] <- lapply(
    sub_mat,
    as.numeric
  )


  y <- factor(
    sub_meta$Group
  )


  # ----------------------------------
  # Random Forest
  # ----------------------------------

  set.seed(
    123
  )


  rf_model <- randomForest::randomForest(

    x =
      sub_mat,

    y =
      y,

    ntree =
      1000,

    importance =
      TRUE
  )


  # Out-of-bag classification error
  oob_error <- rf_model$err.rate[
    nrow(
      rf_model$err.rate
    ),
    "OOB"
  ] *
    100


  # ----------------------------------
  # Feature importance
  # ----------------------------------

  imp_df <- as.data.frame(

    randomForest::importance(
      rf_model
    )
  ) %>%

    tibble::rownames_to_column(
      "Feature"
    )


  if (
    "MeanDecreaseAccuracy" %in%
      colnames(imp_df)
  ) {

    imp_df <- imp_df %>%

      rename(
        Importance =
          MeanDecreaseAccuracy
      )

  } else if (
    "MeanDecreaseGini" %in%
      colnames(imp_df)
  ) {

    imp_df <- imp_df %>%

      rename(
        Importance =
          MeanDecreaseGini
      )

  } else {

    stop(
      paste0(
        "No MeanDecreaseAccuracy or ",
        "MeanDecreaseGini was found in ",
        "the Random Forest importance output."
      )
    )
  }


  imp_df <- imp_df %>%

    mutate(

      Gene =
        clean_gene_name(
          Feature
        )
    )


  # Retain the maximum RF importance for each
  # cleaned gene name.
  imp_gene <- imp_df %>%

    group_by(
      Gene
    ) %>%

    slice_max(
      order_by =
        Importance,
      n = 1,
      with_ties = FALSE
    ) %>%

    ungroup() %>%

    transmute(

      Gene,

      Importance,

      Representative_feature =
        Feature
    ) %>%

    arrange(
      desc(
        Importance
      )
    ) %>%

    slice_head(
      n = top_n
    )


  # ----------------------------------
  # Determine the group with the
  # higher mean abundance
  # ----------------------------------

  direction_df <- sub_mat %>%

    tibble::rownames_to_column(
      "ID"
    ) %>%

    left_join(

      sub_meta %>%
        select(
          ID,
          Group
        ),

      by = "ID"
    ) %>%

    pivot_longer(

      cols =
        -c(
          ID,
          Group
        ),

      names_to =
        "Feature",

      values_to =
        "Abundance"
    ) %>%

    mutate(

      Gene =
        clean_gene_name(
          Feature
        )
    ) %>%

    group_by(
      Gene,
      Group
    ) %>%

    summarise(

      mean_abundance =
        mean(
          Abundance,
          na.rm = TRUE
        ),

      .groups =
        "drop"
    ) %>%

    pivot_wider(

      names_from =
        Group,

      values_from =
        mean_abundance,

      values_fill =
        0
    ) %>%

    mutate(

      Higher_group =
        case_when(

          .data[[g2]] >
            .data[[g1]] ~
            g2,

          .data[[g1]] >
            .data[[g2]] ~
            g1,

          TRUE ~
            "Equal"
        )
    ) %>%

    select(
      Gene,
      Higher_group
    )


  result <- imp_gene %>%

    left_join(
      direction_df,
      by = "Gene"
    ) %>%

    mutate(

      Metal =
        metal_name,

      Metal = factor(
        Metal,
        levels = c(
          "Copper",
          "Zinc",
          "Iron"
        )
      ),

      Comparison =
        paste0(
          g1,
          " vs ",
          g2
        ),

      OOB_error_percent =
        oob_error
    )


  return(
    result
  )
}


# ======================================================
# 11. Run three-metal analysis for one pairwise comparison
# ======================================================

run_pairwise_three_metals_plot <- function(
    g1,
    g2,
    top_n = 6
) {


  res_list <- lapply(

    metal_order,

    function(
        metal_name
    ) {


      run_rf_one_metal_pair(

        metal_name =
          metal_name,

        keyword_pattern =
          metal_keywords[[
            metal_name
          ]],

        g1 =
          g1,

        g2 =
          g2,

        top_n =
          top_n
      )
    }
  )


  names(
    res_list
  ) <- metal_order


  res_list <- purrr::compact(
    res_list
  )


  if (length(res_list) == 0) {

    warning(
      paste0(
        g1,
        " vs ",
        g2,
        ": no Random Forest result generated."
      )
    )

    return(
      NULL
    )
  }


  plot_df <- bind_rows(
    res_list
  ) %>%

    mutate(

      Metal =
        factor(
          Metal,
          levels = c(
            "Copper",
            "Zinc",
            "Iron"
          )
        ),

      Higher_group =
        factor(
          Higher_group,
          levels = c(
            "ITM",
            "OTM1",
            "OTM2",
            "Equal"
          )
        ),

      # Unique internal label prevents duplicated
      # gene names across metal facets.
      Gene_plot =
        paste(
          Metal,
          Gene,
          sep = "___"
        )
    )


  # Set gene order separately within metal categories
  gene_levels <- plot_df %>%

    group_by(
      Metal
    ) %>%

    arrange(
      Importance,
      .by_group = TRUE
    ) %>%

    pull(
      Gene_plot
    ) %>%

    unique()


  plot_df <- plot_df %>%

    mutate(

      Gene_plot =
        factor(
          Gene_plot,
          levels = gene_levels
        )
    )


  group_cols_pair <- c(

    ITM =
      "#8ECFC9",

    OTM1 =
      "#FFBE7A",

    OTM2 =
      "#FA7F6F",

    Equal =
      "grey60"
  )


  # ----------------------------------
  # Plot
  # ----------------------------------

  p <- ggplot(

    plot_df,

    aes(
      x = Gene_plot,
      y = Importance
    )

  ) +


    geom_segment(

      aes(
        xend = Gene_plot,
        y = 0,
        yend = Importance
      ),

      color = "grey72",

      linewidth = 0.65
    ) +


    geom_point(

      aes(
        color = Higher_group
      ),

      size = 3.2
    ) +


    coord_flip(
      clip = "off"
    ) +


    facet_grid(

      Metal ~ .,

      scales = "free_y",

      space = "free_y"
    ) +


    scale_x_discrete(

      labels = function(x) {

        str_replace(
          x,
          "^.*___",
          ""
        )
      }
    ) +


    scale_color_manual(

      values =
        group_cols_pair,

      drop =
        FALSE
    ) +


    labs(

      title =
        paste0(
          g1,
          " vs ",
          g2
        ),

      x =
        NULL,

      y =
        "RF importance",

      color =
        "Higher"
    ) +


    theme_bw(
      base_size = 10
    ) +


    theme(

      panel.grid.minor =
        element_blank(),

      panel.grid.major.y =
        element_blank(),

      panel.grid.major.x =
        element_line(
          linewidth = 0.2,
          color = "grey88"
        ),


      strip.background =
        element_rect(
          fill = "white",
          color = "black",
          linewidth = 0.45
        ),

      strip.text.y =
        element_text(
          face = "bold",
          size = 9
        ),


      plot.title =
        element_text(
          face = "bold",
          hjust = 0.5,
          size = 11
        ),


      axis.title.x =
        element_text(
          face = "bold",
          size = 9
        ),

      axis.text.y =
        element_text(
          face = "bold",
          color = "black",
          size = 7.5,
          margin =
            ggplot2::margin(
              r = 2
            )
        ),

      axis.text.x =
        element_text(
          color = "black",
          size = 7.5
        ),


      legend.position =
        "bottom",

      legend.title =
        element_text(
          size = 8,
          face = "bold"
        ),

      legend.text =
        element_text(
          size = 8
        ),

      legend.key.size =
        grid::unit(
          0.35,
          "cm"
        ),


      plot.margin =
        ggplot2::margin(
          3,
          8,
          3,
          8
        )
    )


  return(

    list(

      plot = p,

      data = plot_df
    )
  )
}


# ======================================================
# 12. Generate three pairwise Random Forest analyses
# ======================================================

rf_ITM_vs_OTM1 <- run_pairwise_three_metals_plot(

  g1 = "ITM",

  g2 = "OTM1",

  top_n = top_n
)


rf_ITM_vs_OTM2 <- run_pairwise_three_metals_plot(

  g1 = "ITM",

  g2 = "OTM2",

  top_n = top_n
)


rf_OTM1_vs_OTM2 <- run_pairwise_three_metals_plot(

  g1 = "OTM1",

  g2 = "OTM2",

  top_n = top_n
)


# ======================================================
# 13. Display individual pairwise plots
# ======================================================

if (!is.null(
  rf_ITM_vs_OTM1
)) {

  print(
    rf_ITM_vs_OTM1$plot
  )
}


if (!is.null(
  rf_ITM_vs_OTM2
)) {

  print(
    rf_ITM_vs_OTM2$plot
  )
}


if (!is.null(
  rf_OTM1_vs_OTM2
)) {

  print(
    rf_OTM1_vs_OTM2$plot
  )
}


# ======================================================
# 14. Combine into Figure 2 Random Forest panel
# ======================================================

plot_list <- list(

  rf_ITM_vs_OTM1,

  rf_ITM_vs_OTM2,

  rf_OTM1_vs_OTM2

) %>%

  purrr::compact()


if (length(plot_list) > 0) {


  p_all <- patchwork::wrap_plots(

    lapply(
      plot_list,
      function(x) {
        x$plot
      }
    ),

    nrow = 1
  ) +


    patchwork::plot_annotation(

      title =
        "D90 Fe/Cu/Zn-related MRG biomarkers",

      subtitle =
        paste0(
          "Random Forest importance; colors indicate ",
          "the group with higher mean abundance"
        ),

      theme =
        theme(

          plot.title =
            element_text(
              face = "bold",
              hjust = 0.5,
              size = 13
            ),

          plot.subtitle =
            element_text(
              hjust = 0.5,
              size = 9
            )
        )
    )


  print(
    p_all
  )
}


# ======================================================
# 15. Combine Random Forest results
# ======================================================

rf_results_list <- list()


if (!is.null(
  rf_ITM_vs_OTM1
)) {

  rf_results_list[[
    "ITM_vs_OTM1"
  ]] <-
    rf_ITM_vs_OTM1$data
}


if (!is.null(
  rf_ITM_vs_OTM2
)) {

  rf_results_list[[
    "ITM_vs_OTM2"
  ]] <-
    rf_ITM_vs_OTM2$data
}


if (!is.null(
  rf_OTM1_vs_OTM2
)) {

  rf_results_list[[
    "OTM1_vs_OTM2"
  ]] <-
    rf_OTM1_vs_OTM2$data
}


if (length(rf_results_list) > 0) {

  rf_all_results <- bind_rows(
    rf_results_list
  )

} else {

  rf_all_results <- tibble()
}


# ======================================================
# 16. Export Random Forest results
# ======================================================

if (nrow(rf_all_results) > 0) {


  write.csv(

    rf_all_results %>%

      select(
        Comparison,
        Metal,
        Gene,
        Representative_feature,
        Importance,
        Higher_group,
        OOB_error_percent
      ),

    file.path(
      output_dir,
      "D90_MRG_randomForest_top6_all_comparisons.csv"
    ),

    row.names = FALSE
  )


  # OOB error summary
  oob_summary <- rf_all_results %>%

    distinct(
      Comparison,
      Metal,
      OOB_error_percent
    ) %>%

    arrange(
      Comparison,
      Metal
    )


  write.csv(

    oob_summary,

    file.path(
      output_dir,
      "D90_MRG_randomForest_OOB_error_summary.csv"
    ),

    row.names = FALSE
  )
}


# ======================================================
# 17. Save combined figure
# ======================================================

if (
  exists("p_all") &&
    length(plot_list) > 0
) {


  ggsave(

    filename = file.path(
      output_dir,
      "Figure_D90_MRG_randomForest_top6.pdf"
    ),

    plot = p_all,

    width = 15,

    height = 6.5
  )


  ggsave(

    filename = file.path(
      output_dir,
      "Figure_D90_MRG_randomForest_top6.png"
    ),

    plot = p_all,

    width = 15,

    height = 6.5,

    dpi = 600
  )
}


# ======================================================
# 18. Print selected biomarkers
# ======================================================

cat(
  "\n================ ITM vs OTM1 =================\n"
)


if (!is.null(
  rf_ITM_vs_OTM1
)) {

  print(

    rf_ITM_vs_OTM1$data %>%

      select(
        Metal,
        Gene,
        Importance,
        Higher_group,
        OOB_error_percent
      )
  )
}


cat(
  "\n================ ITM vs OTM2 =================\n"
)


if (!is.null(
  rf_ITM_vs_OTM2
)) {

  print(

    rf_ITM_vs_OTM2$data %>%

      select(
        Metal,
        Gene,
        Importance,
        Higher_group,
        OOB_error_percent
      )
  )
}


cat(
  "\n================ OTM1 vs OTM2 =================\n"
)


if (!is.null(
  rf_OTM1_vs_OTM2
)) {

  print(

    rf_OTM1_vs_OTM2$data %>%

      select(
        Metal,
        Gene,
        Importance,
        Higher_group,
        OOB_error_percent
      )
  )
}


# ======================================================
# 19. Finish
# ======================================================

cat(
  "\n========================================\n"
)

cat(
  "D90 MRG Random Forest analysis completed.\n"
)

cat(
  "top_n = ",
  top_n,
  "\n",
  sep = ""
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