#' Compute ssGSEA profiles and scores for all samples in a matrix
#'
#' This function applies single-sample GSEA (ssGSEA) across all columns (samples) in an expression matrix.
#' It returns both the full enrichment score profiles and the final ssGSEA scores.
#'
#' @param mat A numeric matrix of gene expression data with genes as rows and samples as columns. Row names must be gene identifiers.
#' @param gene.set A character vector of gene identifiers representing the gene set of interest.
#' @param alpha Exponent parameter used for the weighting of ranks (default: 0.25). Higher values give more weight to highly expressed genes.
#' @param norm Logical; if TRUE (default), the enrichment profile is normalized using min-max scaling as in GSVA.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{profiles}{A matrix of full ssGSEA enrichment profiles, with one column per sample.}
#'   \item{scores}{A numeric vector of ssGSEA scores (the final enrichment value per sample).}
#' }
#'
#' @examples
#' mat <- matrix(rnorm(1000), nrow = 100)
#' rownames(mat) <- paste0("Gene", 1:100)
#' gene.set <- paste0("Gene", sample(1:100, 10))
#' result <- ssGSEA_matrix(mat, gene.set)
#' head(result$scores)
#'
#' @export
ssGSEA_matrix <- function(mat, gene.set, alpha = 0.25, norm = TRUE) {

  profiles <- apply(mat, 2, function(s) {
    ssGSEA(s, gene.set, alpha, norm)
  })

  scores <- apply(mat, 2, function(s) {
    ssgsea_es_final(s, gene.set, alpha, norm)
  })

  return(list(
    profiles = profiles,
    scores = scores
  ))
}




#' Single-sample Gene Set Enrichment Analysis (ssGSEA)
#'
#' Computes the ssGSEA profile and final enrichment score for a single gene expression vector
#' and a given gene set.
#'
#' @param e A named numeric vector of gene expression values. Names must correspond to gene identifiers.
#' @param gs A character vector of genes defining the gene set.
#' @param a Numeric value indicating the exponent alpha used in the weighting scheme (default = 0.25).
#' @param scale Logical indicating whether to normalize the enrichment profile using min-max normalization (default = TRUE).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{fp}{The full running enrichment score profile (numeric vector).}
#'   \item{fp_tail}{The final ssGSEA enrichment score (last value of the profile).}
#' }
#'
#' @examples
#' expr <- c(GENE1 = 3.2, GENE2 = 1.5, GENE3 = 2.8, GENE4 = 0.9)
#' gs <- c("GENE1", "GENE3")
#' result <- ssGSEA(expr, gs)
#' result$fp_tail  # Final ssGSEA score
#'
#' @export
ssGSEA <- function(e, gs, a = 0.25, scale = TRUE){

  helper <- function(exp, gene.set, alpha = 0.25, norm = TRUE) {

    exp <- exp[!is.na(exp)]

    rgenes <- sort(exp, decreasing = TRUE)
    gnames <- names(rgenes)

    N <- length(rgenes)
    Nh <- length(intersect(gene.set, gnames))

    if (Nh == 0) return(rep(0, N))

    # Indicator
    indicator <- gnames %in% gene.set

    # Weighting
    ranks <- seq_len(N)
    weights <- ranks^alpha

    # Compute Phit and Pmiss vectors
    sum_hit <- sum(weights[indicator])
    sum_miss <- N - Nh

    Phit <- ifelse(indicator, weights / sum_hit, 0)
    Pmiss <- ifelse(!indicator, 1 / sum_miss, 0)

    # Running enrichment score profile
    running_ES <- cumsum(Phit - Pmiss)

    # Optional normalization (GSVA min-max style)
    if (norm) {
      min_ES <- min(running_ES)
      max_ES <- max(running_ES)
      if (max_ES != min_ES) {
        running_ES <- (running_ES - min_ES) / (max_ES - min_ES)
      } else {
        running_ES <- rep(0, length(running_ES))
      }
    }

    return(running_ES)
  }

  ##
  fp <- helper(exp=e, gene.set = gs, a=a, norm=scale)
  fp2 <- tail(fp, 1)

  return(list("fp"=fp, "fp_tail"=fp2))
}


