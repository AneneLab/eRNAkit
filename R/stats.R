#' Compute Tissue/Cell Type Specificity Scores
#'
#' This function calculates specificity scores for eRNA or genes based on their
#' expression profiles across multiple tissues or cell types. The specificity
#' score is computed using Shannon entropy and a logarithmic correction.
#'
#' @param x A data frame or matrix where rows represent eRNAs and columns
#'          represent different tissues or cell types. Only numeric columns
#'          (expression values) are used for calculations.
#'
#' @return A data frame with three columns:
#' \describe{
#'   \item{ID}{The identifier (row names of input data).}
#'   \item{Specificity_Score}{The calculated specificity score for each row.}
#'   \item{Is_Specific}{The tissue/cell type where the transcript is specific, or FALSE if not specific.}
#' }
#'
#' @details
#' The function follows these steps:
#' \enumerate{
#'   \item Compute expression ratios by normalizing each eRNA's expression across tissues/cell types.
#'   \item Calculate Shannon entropy for each eRNA.
#'   \item Compute the specificity score as \eqn{log_2(N) - H}, where \(N\) is the number of cell types
#'         and \(H\) is Shannon entropy.
#'   \item Define specificity by checking if the highest expression ratio is at least 2× the second highest
#'         and if the specificity score is greater than 1.
#' }
#'
#' @export
entropyTS <- function(x){

  numX <- x[, sapply(x, is.numeric),
            drop = FALSE]

  # S1 - expression ratios
  to_exp <- rowSums(numX)
  eR <- sweep(numX, 1, to_exp, "/")

  # S2 - Shannon entropy
  shannon <- apply(eR, 1, function(x) {
    -sum(x[x > 0] * log2(x[x > 0]))
  })

  # S3 - Specificity
  N <- ncol(numX)
  SS <- log2(N) - shannon

  # S4 - Define specificity label
  is_specific <- sapply(1:nrow(eR), function(i) {
    x2 <- setNames(as.numeric(eR[i, ]), colnames(eR))

    if(all(is.na(x2))) {
      return(FALSE)
    }

    sR <- sort(x2, decreasing = TRUE)
    rName <- rownames(eR)[i]
    score <- SS[rName == rownames(numX)]

    if(sR[1] / sR[2] > 2 && score > 1) {
      return(names(sR[1]))
    }
    return(FALSE)
  })

  # S5 - Bind results
  res <- data.frame(
    ID = rownames(numX),
    Specificity_Score = SS,
    Is_Specific = is_specific
  )

  return(res)

}

#' Count to log2 count per million "CPM"
#'
#' Convert count to log2 CPM on matrix
#'
#' @param df: Data frame
#' @param side: Operate on row-1 or column-2
#'
#' @return Data frame of converted counts
#' @export
#'
#' @examples
#' log2_count2cpm(df, 2)
log2_count2cpm <- function(df="data", side="r:1.c:2") {
  # This two offsets were got from limma:voom paper
  # counts_offset = 0.5
  # library_offset = 1
  apply(df, side, function(x) log2((x + 0.5 /(sum(x)+1))*1000000))
}


#' Detect outliers in a matrix.
#'
#' Function to apply the two outlier function to a matrix.
#' This is specific to eRNAkit usage
#'
#' @param m matrix of expression values gene/eRNA in column and sample in rows.
#'          It expects the row.names to set to gene IDs.
#' @param type String indicating the p-value estimation method.
#'             "pnorm" or "permutation"
#' @param name The name to use for the new gene column.
#'             This option is to differentiate eRNA/mRNA runs.
#'
#' @returns Table of results based on the output of the outlier functions.
#' @export
#'
#' @examples
#' outlier_matrix(matrix(1:9, nrow = 3, ncol = 3), "pnorm")
outlier_matrix <- function(m, type="pnorm", name="eRNA") {

  res <- lapply(rownames(m), function(id) {
    if(type == "pnorm"){
      df <- eRNAkit::pnorm_outliers(m[id, ])
    } else {
      df <- eRNAkit::permutation_outliers(m[id, ])
    }

    df[[name]] <- id
    return(df)
  })

  # NOTE: Adjusted p-values here
  res <- do.call(rbind, res)
  res$p_adj <- p.adjust(res$P, method = "BH")
  res$outlier <- (res$p_adj <= 0.05)

  return(res)
}

#' Detect outliers in a named vector using normalized z-score.
#'
#' Function to detect outliers in a named vector of values.
#' It uses CDF of pnorm to get p-values
#'
#' @param x A named vector of values.
#'
#' @returns A table of all the results including intermediate tables.
#' @export
#'
#' @examples
#' pnorm_outliers(c(2, 4, 7, 8,9))
pnorm_outliers <- function(x) {

  median <- median(x, na.rm = TRUE)
  mad <- eRNAkit::compute_mad(x)
  # by 0.6745 normalizes to standard normal z-scores.
  z <- 0.6745 * (x - median) / mad

  # Compute p-values using normal approximation
  # 2-tailed not suitable for this case.
  p_norm <- 1 - stats::pnorm(z)
  maxs <- eRNAkit::max_cols(x)

  return(data.frame(
    Sample = names(x),
    Value = x,
    Z = z,
    P = p_norm,
    max = maxs[1],
    max2 = maxs[2]
  ))
}


#' Detect outliers in a named vector using permutations.
#'
#' Function to detect outliers in a named vector of values.
#' It uses permutations to get empirical p-values.
#'
#' @param x A named vector of values.
#' @param n Number of permutations to run.
#'
#' @returns A table of all the results including intermediate tables.
#' @export
#'
#' @examples
#' permutation_outliers(c(2, 4, 7, 8,9), n=100)
permutation_outliers <- function(x, n = 1000) {

  set.seed(856)
  median <- median(x, na.rm = TRUE)
  mad <- eRNAkit::compute_mad(x)
  z <- 0.6745 * (x - median) / mad
  p_perm <- numeric(length(x))

  # Speed optimized permutation
  perm_samp <- replicate(n, sample(x))
  perm_median <- apply(perm_samp, 2, median, na.rm = TRUE)
  perm_mad <- apply(perm_samp, 2, eRNAkit::compute_mad)
  perm_z <- 0.6745 * sweep(perm_samp, 2, perm_median,"-") / perm_mad


  # Compute p-values using normal approximation
  # 2-tailed not suitable for this case.
  p_perm <- rowMeans(sweep(perm_z, 1, z, `>=`))
  p_perm <- matrix(p_perm, ncol = 1)
  rownames(p_perm) <- names(x)

  maxs <- eRNAkit::max_cols(x)

  return(data.frame(
    Sample = names(x),
    Value = x,
    Z = z,
    P = p_perm,
    max = maxs[1],
    max2 = maxs[2]
  ))
}

