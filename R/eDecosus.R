#' Make eRNA based signature for fine cell typoe deconvolution.
#'
#' This function creates the signature for eDecosus from an \code{eRNAkitDB} database object
#' based on the specified \code{type} ("cell" or "organ"), minimum specificity,
#' and CPM thresholds. It returns the top \code{n} entries per group (cell or organ)
#' by two ranking methods:
#' 1. Specificity Score then CPM
#' 2. CPM then Specificity Score
#' The results from both methods are combined and deduplicated.
#'
#' @param eRNAkitDB A list containing at least the \code{EOC} data frame with
#'   eRNA signatures and their metadata.
#' @param type Character; one of \code{"cell"} or \code{"organ"} specifying which
#'   signature type to filter.
#' @param specificity Numeric; minimum Specificity Score to retain rows
#'   (default: 2).
#' @param cpm Numeric; minimum CPM (counts per million) to retain rows
#'   (default: 20).
#' @param n Integer; number of top rows to return per group after filtering and sorting
#'   (default: 10).
#'
#' @return A data frame containing filtered and combined top eRNA signatures
#'   grouped by the specified \code{type} (cell or organ).
#'
#' @examples
#' \dontrun{
#' eRNAkitDB <- readRDS("path/to/eRNAkitDB.rds")
#' eSig <- make_eRNASig(eRNAkitDB, type = "cell",
#'                                                  specificity = 2,
#'                                                  cpm = 20,
#'                                                  n = 10) }
#' @export
make_eRNASig <- function(db=eRNAkitDB, type = "cell",
                         specificity = 2, cpm = 20,
                         n = 10) {
  #
  spec <- db$EOC
  spec <- spec[spec$type %in% type & spec$Is_Specific != "FALSE" & spec$Is_Specific == spec$space, ]
  spec$pair <- paste0(spec$E, "_", spec$Is_Specific)
  spec <- spec[spec$Specificity_Score >= specificity & spec$CPM >= cpm, ]
  spec <- split(spec, spec$space)

  # Helper
  topS <- function(df, side = c("sc", "cs"), n = 10) {
    method <- match.arg(side)

    if (method == "sc") {
      # specificity first and then cpm
      df <- df[order(-df$Specificity_Score, -df$CPM), ]
    } else if (method == "cs") {
      # cpm first and the specificity
      df <- df[order(-df$CPM, -df$Specificity_Score), ]
    }

    head(df, n)
  }

  # Apply both filters
  sig1 <- lapply(spec, topS, side="sc", n=n)
  sig2 <- lapply(spec, topS, side="cs", n=n)

  eRNASig <- do.call(rbind, c(sig1, sig2))
  eRNASig <- eRNASig[!duplicated(eRNASig$pair), ]

  return(eRNASig)
}



#' Compute ssGSEA scores for multiple gene sets across all samples
#'
#' @param mat Numeric matrix of expression data (genes in rows, samples in columns).
#' @param gene.sets Named list of character vectors (each vector = genes in a set).
#' @param alpha Weight exponent for ranks (default 0.25).
#' @param norm Logical; normalize enrichment profile before taking score (default TRUE).
#'
#' @return Data frame with columns:
#'   \describe{
#'     \item{gene_set}{Name of the gene set.}
#'     \item{sample}{Sample name (column in \code{mat}).}
#'     \item{score}{Numeric ssGSEA score.}
#'   }
#' @export
GSEAmm <- function(mat, gene.sets, alpha = 0.25, norm = TRUE) {

  scores <- vector("list", length(gene.sets))
  names(scores) <- names(gene.sets)

  for (set in names(gene.sets)) {
    gs <- gene.sets[[set]]

    # One score per sample (final ES value)
    scr <- apply(mat, 2, function(s) {
      ssGSEA(s, gs, a = alpha, scale = norm)$fp_tail
    })

    scores[[set]] <- data.frame(
      gene_set = set,
      sample = colnames(mat),
      score = scr,
      stringsAsFactors = FALSE
    )
  }

  scores <- do.call(rbind, scores)
  rownames(scores) <- NULL
  return(scores)
}



#' Compute ssGSEA profiles and scores for all samples in a matrix for a single gene set
#'
#' This function applies single-sample GSEA (ssGSEA) across all columns (samples) in an expression matrix.
#' It returns JUST final ssGSEA scores. If you need the full profile then run the ssGSEA itslf with apply.
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
GSEAms <- function(mat, gene.set, alpha = 0.25, norm = TRUE) {

  res <- apply(mat, 2, function(s) {
    ssGSEA(s, gene.set, a = alpha, scale = norm)
  })

  profiles <- do.call(cbind, lapply(res, `[[`, "fp"))
  scores   <- vapply(res, `[[`, numeric(1), "fp_tail")

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


