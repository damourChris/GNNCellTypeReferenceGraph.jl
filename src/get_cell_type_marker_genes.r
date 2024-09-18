get_cell_type_marker_genes <- function(eset) {
  # Filter out genes with p-values bigger than threshold
  anova_significant_threshold <- 0.0001
  tuxey_significant_threshold <- 10e-9

  expression_matrix <- Biobase::exprs(eset)

  # Log transform the expression matrix
  expression_matrix <- log2(expression_matrix + 1)
  rownames(expression_matrix) <- Biobase::fData(eset)$ensembl_id

  sample_group <- factor(eset[["cell ontology:ch1"]])

  print("Running ANOVA test...")
  anova_results <- apply(expression_matrix, 1, function(gene_expr) {
    anova(aov(gene_expr ~ sample_group))
  })

  p_values <- lapply(anova_results, function(x) {
    x[["Pr(>F)"]][1]
  })

  adjusted_p_values <- p.adjust(p_values, method = "fdr")
  names(adjusted_p_values) <- rownames(expression_matrix)

  significant_genes <-
    names(adjusted_p_values)[adjusted_p_values < anova_significant_threshold]

  # Reduce the original expression matrix to only the significant genes
  expression_matrix <- expression_matrix[significant_genes, ]

  print("Running Tukey test...")
  tukey_results <- apply(expression_matrix, 1, function(gene_expr) {
    model <- aov(gene_expr ~ sample_group)
    TukeyHSD(model)
  })

  names(tukey_results) <- significant_genes


  print("Identifying cell types...")
  cell_types_gene_pairings <- lapply(tukey_results, function(res) {
    tukey_pvals <- res$sample_group[, 4]
    sig_comps <- res$sample_group[tukey_pvals < tuxey_significant_threshold, ]

    types <- lapply(rownames(sig_comps), function(row) {
      strsplit(row, "-")
    })

    # Count frequencies
    frequency_table <- table(unlist(types))

    # Identify the most frequent string
    most_frequent_string <- names(which.max(frequency_table))

    # Make sure that the most frequent string is not empty or NULL
    if (most_frequent_string == "" || is.null(most_frequent_string)) {
      return(NULL)
    }

    most_frequent_string

    # This should be recorded only if the most frequent string count is
    # bigger than 1 and all the other strings counts are 1
    if (frequency_table[most_frequent_string] > 1 &&
      all(frequency_table[
        -which(names(frequency_table) == most_frequent_string)
      ]) == 1) {
      return(most_frequent_string)
    } else {
      return(NULL)
    }
  })

  # Filter out the NULL values
  cell_types_gene_pairings_filter <-
    cell_types_gene_pairings[!sapply(cell_types_gene_pairings, is.null)]

  known_cell_types <- unique(unlist(cell_types_gene_pairings))

  cell_type_marker_genes <- list()


  for (value in known_cell_types) {
    indices <- which(unlist(cell_types_gene_pairings_filter) == value)


    if (length(indices) == 0) {
      next
    }


    new_list <- names(cell_types_gene_pairings_filter[indices])


    cell_type_marker_genes[[value]] <- new_list
  }

  return(cell_type_marker_genes)
}
