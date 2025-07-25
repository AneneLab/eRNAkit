#' Prioritize Functional Regulatory Variants Using eRNA–mRNA Interactions.
#'
#' This function identifies and prioritizes putative functional regulatory non-coding variants
#' by integrating variant data with enhancer RNA (eRNA) to mRNA interaction networks.
#' It processes a VCF file, maps variants to eRNAs, identifies downstream regulatory targets,
#' and links variants to potential gene regulatory effects using the `eRNAkit` database.
#'
#' @param vcf Path to a VCF file containing variant calls. The file names must match the column names in E and G
#' @param E Data frame of eRNA expression. Row names must match the eRNA IDs from eRNAkitDB and column names must match vcf file names.
#' @param G Data frame of eRNA expression. Row names must match the Ensemble gene IDs and column names must match vcf file names.
#' @param n Minimum number of samples to contain the a variant for it to be tested. (default: `2`).
#' @param fdr False discovery rate threshold for filtering results (default: `0.1`).
#' @param db eRNAkitDB object.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{var}{Preprocessed variant data.}
#'   \item{var2erna}{Mapped variants to eRNAs.}
#'   \item{einteraction}{eRNA-mRNA interactions involving those variants.}
#'   \item{varReg}{Final regulatory variant prioritization results. With hit - significant hit, E and G full test for checks}
#' }
#'
#' @export
cisR <- function(vcf="path/to/vcf", E="eRNA", G="mRNA",
                 n=2, fdr=0.1, db=eRNAkitDB){
  cisRdb <- list()
  cisRdb[["var"]] <- eRNAkit::process_vcf(vcf, n=n)
  cisRdb[["var2erna"]] <- eRNAkit::var2erna(cisRdb$var, db$core)
  cisRdb[["einteraction"]] <- eRNAkit::link2rr(cisRdb$var2erna, db$R2R)
  cisRdb[["varReg"]] <- eRNAkit::varReg(cisRdb$var, cisRdb$var2erna,
                               cisRdb$einteraction, E=E, G=G, t=fdr)
  return(cisRdb)
}



#' Variant Regulatory Impact Analysis via eRNA and mRNA Expression
#'
#' This function prioritizes regulatory variants by assessing their impact on
#' enhancer RNA (eRNA) and downstream mRNA expression. It integrates variant
#' carrier sample information, eRNA annotations, and eRNA–mRNA interaction data
#' to identify variants associated with significant expression changes in both
#' eRNAs and their putative mRNA targets.
#'
#' NOTE: This is mainly an internal function for cisR.
#'
#' @param output A \code{data.frame} of variant-sample mappings,
#' from \code{process_vcfs()}, containing at least \code{id} (variant
#'   IDs) and \code{sample} (carrier sample IDs).
#' @param output2 A \code{data.frame} linking variants to eRNAs, from
#'   \code{var2erna()}, containing columns \code{id} (variant IDs) and \code{E}
#'   (eRNA IDs).
#' @param output3 A \code{data.frame} representing eRNA–mRNA interactions, from
#'   \code{linkE2RD()}, with at least columns \code{E} (eRNA IDs) and \code{G}
#'   (mRNA IDs).
#' @param E Numeric matrix or \code{data.frame} of eRNA expression (rows = eRNAs,
#'   columns = samples). Row names should be eRNAkit IDs "enxxx".
#' @param G Numeric matrix or \code{data.frame} of mRNA expression (rows = mRNAs,
#'   columns = samples). Row names should be Ensemble IDs to match output1-3.
#'   @param t False discovery rate threshold for filtering results (default: `0.1`).
#'
#' @return A named \code{list} with components:
#' \describe{
#'   \item{hits}{\code{data.frame} of variants with significant differential
#'     expression in both eRNAs and linked mRNAs, including log fold change,
#'     moderated t-statistics, adjusted p-values, and direction of effect.}
#'   \item{E}{Full \code{data.frame} of eRNA-level differential expression
#'     statistics for all tested variants.}
#'   \item{G}{Full \code{data.frame} of mRNA-level differential expression
#'     statistics for only significant eRNAs in E.}
#' }
#'
#' @details
#' The function performs the following key steps:
#' \enumerate{
#'   \item Filters eRNA and mRNA expression data for features with sufficient
#'     expression (>1 count in more than 2 samples), then log2-transforms the
#'     data.
#'   \item Estimates global empirical Bayes variance priors for eRNA and mRNA
#'     expression using \code{gprior()}.
#'   \item For each variant in \code{output2}, tests differential expression of
#'     linked eRNAs between variant carriers and non-carriers using
#'     \code{limma::lmFit} and variance shrinkage via \code{shrinkFit()}.
#'   \item Extracts significantly altered eRNAs (FDR ≤ 0.1) and, for each,
#'     tests linked mRNAs (from \code{output3}) for differential expression
#'     under the same model.
#'   \item Compiles final hits as variants with significant effects on both
#'     eRNAs and their downstream mRNA targets, applying a combined FDR
#'     threshold (~1%).
#' }
#'
#' The use of an FDR cutoff of 0.1 at both eRNA and mRNA levels results in an
#' overall approximate FDR of 0.01 for the final hits.
#'
#' @seealso \code{\link{process_vcfs}}, \code{\link{var2erna}},
#'   \code{\link{linkE2RD}}, \code{\link{gprior}}, \code{\link{shrinkFit}}
#'
#' @importFrom limma lmFit
#' @export
varReg <- function(output, output2,
                   output3, E="eRNA", G="mRNA", t=0.1){

  # Check
  eRNAkit::checkN(E, "E")
  eRNAkit::checkN(G, "G")

  # Filter and scale
  log_E <- log2(E[rowSums(E > 1) > 2, ] + 1)
  log_G <- log2(G[rowSums(G > 1) > 2, ] + 1)

  # Get global priors
  E.prior <- eRNAkit::ebayes_prior(log_E)
  G.prior <- eRNAkit::ebayes_prior(log_G)

  ## Analyse E
  bout <- lapply(output2$id, function(v){
    varSamples <- output[output$id %in% v, ][["sample"]]
    er <- output2[output2$id %in% v, ][["E"]]

    if (length(er) == 0 ) {
      return(NULL)
    }

    e_exp <- log_E[rownames(log_E) %in% er, ]

    if (nrow(e_exp) == 0 ) {
      return(NULL)
    }

    group <- ifelse(names(e_exp) %in% varSamples, 1, 0)
    design <- model.matrix(~factor(group))
    ##
    fit <- limma::lmFit(e_exp, design)
    fit <- eRNAkit::shrink_lmfit(fit,E.prior["df"], E.prior["s2"])

    fit$var <- v
    fit$ID <- rownames(e_exp)

    return(fit)
  })
  bout <- do.call(rbind, bout)
  bout$sig <- ifelse(bout$local.FDR <= t, "sig", "")
  bout$is <- "E"
  rownames(bout) <- NULL

  # Restrict
  output2b <- bout[bout$sig %in% "sig", ][["var"]]

  ## Analyse linked G
  bout2 <- lapply(output2b, function(v){
    varSamples <- output[output$id %in% v, ][["sample"]]
    er <- output2[output2$id %in% v, ][["E"]]
    eInter <- output3[output3$E %in% er, ]

    if (nrow(eInter) == 0 ) {
      return(NULL)
    }
    idmap <- unique(eInter$G)
    ## Get expression
    e_exp <- log_G[rownames(log_G) %in% idmap, ]

    if (nrow(e_exp) == 0 ) {
      return(NULL)
    }

    group <- ifelse(names(e_exp) %in% varSamples, 1, 0)
    design <- model.matrix(~factor(group))
    ##
    fit <- limma::lmFit(e_exp, design)
    fit <- eRNAkit::shrink_lmfit(fit,E.prior["df"], E.prior["s2"])

    fit$var <- v
    fit$ID <- rownames(e_exp)

    return(fit)
  })
  bout2 <- do.call(rbind, bout2)
  bout2$sig <- ifelse(bout2$local.FDR <= t, "sig", "")
  bout2$is <- "G"
  rownames(bout2) <- NULL

  # NOTE: The use of FDR 0.1 at both level is justified.
  # Final FDR will be 0.1 * 0.1 = 1%
  ghit <- bout2[bout2$sig %in% "sig", ]
  ehit <- bout[bout$var %in% ghit$var, ]
  message(paste("Identifed ", nrow(ehit), "E hits", " with",
                nrow(ghit), " associated G."))

  hits <- rbind(ehit, ghit)
  hits$dir.var <- ifelse(hits$logFC > 0, "up", "down")
  hits <- hits[order(hits$var), ]

  if (nrow(hits) == 0) {
    message("No significant hits.")
  }

  return(list("hits"=hits, "E"=bout, "G"=bout2))
  }


#' Link eRNAs to their Regulatory Targets
#'
#' This function filters the R2R interaction table to return only entries
#' matching a given set of eRNAs (from `output2`).
#'
#' @param output2 A data frame with a column `E` containing eRNA IDs. Out of var2erna
#' @param R2R eRNAkitDB R2R object.
#'
#' @return A data frame with filtered R2R interactions and added `"interaction"` column set to `"R2R"`.
#' @export
#'
#' @examples
#' \dontrun{
#' output2 <- var2erna(output)
#' result <- linkE2RD(output2)
#' }
link2rr <- function(output2, R2R="dd"){
  cols <- c("pair", "E", "G", "Gs", "biotype",
            "confidence", "D2D", "interaction")

  rr <- within(subset(R2R, E %in% output2$E),
               interaction <- "R2R")[cols]

  return(rr) }


#' Map Variants to Overlapping eRNAs
#'
#' This function takes a set of variants (from process_vcf() function)
#' and returns overlaps with enhancer RNAs (eRNAs) using genomic coordinates.
#'
#' @param output A data.frame with a column named `id` in the format "chr:pos:ref:alt".
#' Typically the output of \code{\link{process_vcfs}}.
#' @param eRNA eRNAkitDB core object.
#'
#' @return A data.frame with columns \code{id} (variant ID) and \code{E} (eRNA ID).
#' Returns \code{NULL} if no overlaps are found.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' v <- process_vcf("vcf_folder/")
#' hits <- var2erna(v)
#' head(hits)
#' }
var2erna <- function(output, eRNA="db") {

  vcf_output <- output[!duplicated(output$id), ]

  coord_matrix <- do.call(rbind, strsplit(vcf_output$id, ":", fixed = TRUE))
  tab <- data.frame(
    chr = coord_matrix[, 1],
    start = as.integer(coord_matrix[, 2]),
    end = as.integer(coord_matrix[, 2]),
    strand = "*",
    id = vcf_output$id,
    stringsAsFactors = FALSE
  )

  tab <- eRNAkit::getGRange(tab)
  etab <- eRNAkit::getGRange(eRNA)
  hits <- GenomicRanges::findOverlaps(tab, etab)

  # Match
  matched <- cbind(
    as.data.frame(S4Vectors::mcols(tab)[S4Vectors::queryHits(hits), ]),
    as.data.frame(S4Vectors::mcols(etab)[S4Vectors::subjectHits(hits), ])
  )


  if (nrow(matched) == 0) {
    message("No variant–eRNA intersect found")
    return(NULL)
  }

  # 1:7 here is intentional.
  # It allows the code to fail if wrong input
  matched <- matched[c(1, 7)]
  names(matched) <- c("id", "E")

  return(matched)
}



#' Process VCF files from a folder
#'
#' This function reads all `.vcf` files in a given folder, filters for PASS SNPs (single nucleotide variants),
#' and compiles them into a long-format data frame where each row is a variant found in a sample.
#' Variants are retained only if they are found in at least `n` samples. This function can handle .gz
#'
#' @param folder A character string giving the path to the folder containing `.vcf` files.
#' @param n Minimum number of samples in which a variant must appear to be retained. Default is 2.
#'
#' @return A data frame with columns `id` (variant identifier) and `sample` (sample name).
#' @export
#'
#' @examples
#' \dontrun{
#' df <- process_vcfs("path/to/vcf_folder", n = 3)
#' head(df)
#' }
process_vcf <- function(folder, n=2) {
  vcfs <- list.files(folder,
                     pattern = "\\.vcf$",
                     full.names = TRUE)

  if (length(vcfs) == 0) {
    stop("No VCF files found.")
  }

  ## Do check here for vcfs not empty

  out <- vector("list", length(vcfs))
  names(out) <- sub("\\.vcf$", "", basename(vcfs))


  for (i in seq_along(vcfs)) {
    df <- read.table(vcfs[i], header = FALSE,
                     quote = "\"", comment.char = "#")

    df <- df[df$V7 == "PASS", ]
    df <- df[nchar(as.character(df$V4)) == 1 & nchar(as.character(df$V5)) == 1, ]
    df$id <- paste(df$V1, df$V2, df$V4, df$V5, sep = ":")
    df$sample <- sub("\\.vcf(\\.gz)?$", "", basename(vcfs[i]))
    df <- df[c("id", "sample")]
    out[[i]] <- df
  }

  out <- do.call(rbind, out)
  id_counts <- table(out$id)
  out <- out[out$id %in% names(id_counts[id_counts >= n]), ]
  out$id <- sub("chr", "", out$id)
  return(out)
}



#' Compute Moderated t-Statistics with Empirical Bayes Variance Shrinkage
#'
#' Computes the posterior variance estimate by combining
#' prior variance estimates with residual variance from a linear model fit.
#' It then calculates moderated t-statistics, p-values, and BH correction.
#'
#' @param fit A \code{lmFit} object returned by \code{limma::lmFit} containing
#'   coefficients, standard errors, residual degrees of freedom, and residual variances.
#' @param df Numeric scalar representing the prior degrees of freedom obtained from
#'   empirical Bayes estimation.
#' @param s2 Numeric scalar representing the prior variance estimate obtained from
#'   empirical Bayes estimation.
#'
#' @return A \code{data.frame} with the following columns:
#' \describe{
#'   \item{logFC}{Estimated log fold change of coefficients (2nd column of \code{fit$coefficients}).}
#'   \item{AveExpr}{Average expression values (\code{fit$Amean}).}
#'   \item{t}{Moderated t-statistics computed with shrunk variance.}
#'   \item{P.Value}{Two-sided p-values associated with moderated t-statistics.}
#'   \item{local.FDR}{Adjusted p-values controlling false discovery rate using BH method.}
#'   \item{B}{Log-odds of differential expression (set to NA here).}
#' }
#'
#' @details
#' The posterior variance estimate is computed as a weighted average of
#' the prior variance and the residual variance from the model fit, weighted
#' by their respective degrees of freedom. The moderated t-statistics
#' improve variance estimation stability especially when sample sizes are small.
#'
#' @examples
#' \dontrun{
#' library(limma)
#' dm <- model.matrix(~factor(c(1,1,0,0)))
#' fit <- lmFit(matrix(rnorm(100), ncol=4), dm)
#' fit <- eBayes(fit)
#' prior <- c(s2=fit$s2.prior, df=fit$df.prior)
#' shrinkFit(fit, prior["df"], prior["s2"])
#' }
#'
#' @importFrom stats pt p.adjust
#' @export
shrink_lmfit <- function(fit, df, s2) {
  # Compute posterior variance estimate
  s2.post <- (df * s2 + fit$df.residual * fit$sigma^2) /
    (df + fit$df.residual)

  # Compute moderated t-statistics
  t.mod <- fit$coef / sqrt(fit$stdev.unscaled^2 * s2.post)
  df.total <- fit$df.residual + df
  p.val <- 2 * pt(-abs(t.mod), df = df.total)

  # Table
  out <- data.frame(
    logFC     = fit$coefficients[, 2],
    AveExpr   = fit$Amean,
    t         = t.mod[, 2],
    P.Value   = p.val[, 2],
    local.FDR = p.adjust(p.val[, 2], method = "BH"),
    B         = NA
  )

  return(out)
}


#' Estimate Empirical Bayes Prior Parameters
#'
#' Computes the prior variance (\eqn{s_0^2}) and prior degrees of freedom (\eqn{d_0})
#' for a given expression matrix using empirical Bayes shrinkage, as implemented in the \pkg{limma} package.
#'
#' This function fits a simple intercept-only linear model to the input data and estimates
#' hyper parameters required for moderated t-statistics.
#'
#' Works best with log2 transformed matrix log2(x)
#'
#' @param x A numeric matrix of expression values with features as rows and samples as columns.
#' @return A named numeric vector with elements:
#' \describe{
#'   \item{\code{s2}}{Prior variance (\eqn{s_0^2})}
#'   \item{\code{df}}{Prior degrees of freedom (\eqn{d_0})}
#' }
#' @examples
#' mat <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' gprior(mat)
#'
#' @export
#' @importFrom limma lmFit eBayes
ebayes_prior <- function(x){
  dm <- matrix(1, ncol = 1, nrow = ncol(x))
  fit <- limma::lmFit(x, dm)
  fit <- limma::eBayes(fit)
  return(c("s2" = fit$s2.prior,
           "df" = fit$df.prior))

  }
