#' Identify Potential Regulatory Variants
#'
#' This function predicts potential regulatory non-coding variants based on target mRNA decay and translation efficiency (TE).
#' It is tailored for exploratory use in **single-sample scenarios**, such as variant interpretation, rather than full expression-based
#' prioritisation across a cohort.
#'
#' The input `input` must be a data frame containing only the variants. It must include
#' **three columns** named **chr**, **start**, and **end**. Start and end must be the same value.
#'
#' @param input Data frame containing variants. Must have exactly three columns: `chr`,
#'   `start`, and `end`.
#' @param fdr Numeric. FDR threshold for significance (default is 0.1).
#' @param db eRNAkitDB object.
#'
#' @return A named list with result objects.
#'
#' @examples
#' # Assuming input_df has columns chr, start, end
#' # and eRNAkitDB is loaded with required tables:
#' result <- to_dete(input_df, fdr = 0.05, db = eRNAkitDB)
#'
#' @export
TranCi_explore <- function(input, fdr=0.1, db=eRNAkitDB){

  tab <- input
  tab$chr <- sub("chr", "", tab$chr)
  tab$strand <- "*"
  tab$id <- paste0("row", rownames(tab))

  ## Collect
  cisRdb <- list()
  cisRdb[["var"]] <- tab

  ## Intersect
  tab <- eRNAkit::getGRange(tab)
  etab <- eRNAkit::getGRange(db$core)
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

  matched <- matched[c("mcols.id", "mcols.E")]
  names(matched) <- c("id", "E")


  ## Add to collection
  cisRdb[["var2erna"]] <- matched
  cisRdb[["einteraction"]] <- eRNAkit::link2rr(cisRdb$var2erna, db$R2R)

  dete <- eRNAkit::to_dete(cisRdb$einteraction$E, db, t=fdr)
  cisRdb[["e2decay"]] <- dete$de
  cisRdb[["e2TE"]] <- dete$te

  ## need to do ranking here

  return(cisRdb)
}




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
TranCi <- function(vcf="path/to/vcf", E="eRNA", G="mRNA",
                 n=2, fdr=0.1, db=eRNAkitDB){
  cisRdb <- list()
  cisRdb[["var"]] <- eRNAkit::process_vcf(vcf, n=n)
  cisRdb[["var2erna"]] <- eRNAkit::var2erna(cisRdb$var, db$core)
  cisRdb[["einteraction"]] <- eRNAkit::link2rr(cisRdb$var2erna, db$R2R)

  dete <- eRNAkit::to_dete(cisRdb$einteraction$E, db, t=fdr)
  cisRdb[["e2decay"]] <- dete$de
  cisRdb[["e2TE"]] <- dete$te

  var <- eRNAkit::varReg(cisRdb$var, cisRdb$var2erna,
                         cisRdb$einteraction, E=E, G=G, t=fdr)
  cisRdb$E <- var$E
  cisRdb$G <- var$G

  result <- eRNAkit::cisR_Score(cisRdb)
  result$RAW <- cisRdb


  return(result)
}




#' Score Evidence of Regulatory Variants
#'
#' This function assigns quantitative scores to variants in a the cisRdb interaction database
#' based on multiple layers of evidence, including eRNA impact,eRNA-gene regulation,
#' decay impact, and translation. It outputs a final classification of each variant into
#' one of three evidence classes: "empirical", "predictive", or "weak".
#'
#' @param cisRdb A list containing regulatory interaction data from cisR run, with the following named components:
#' \describe{
#'   \item{einteraction}{Data frame of enhancer-variant-gene interactions with columns \code{id}, \code{E}, and \code{G}.}
#'   \item{var2erna}{Data frame mapping variant IDs to enhancer IDs with columns \code{id} and \code{E}.}
#'   \item{E}{Data frame of enhancer-only analyses, including a \code{sig} column and a \code{trend} column.}
#'   \item{G}{Data frame of enhancer-to-gene analysis, with a \code{sig} column and a \code{trend} column.}
#'   \item{e2TE}{Data frame with translation and decay evidence, including \code{pair}, \code{sig}, and \code{trend} columns.}
#' }
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{variants}}{A data frame with unique variant identifiers and the highest evidence score observed across interactions. Contains columns \code{chr}, \code{start}, \code{ref}, \code{alt}, \code{id}, \code{score}, and \code{s.class} (evidence class).}
#'   \item{\code{evidence}}{The full interaction-level data frame annotated with component evidence and intermediate scores.}
#' }
#'
#' @details
#' The scoring is additive, with each layer contributing as follows:
#' \itemize{
#'   \item Enhancer-only (E.impact): 5 points
#'   \item Enhancer-gene (E2G.impact): 2.5 points
#'   \item Decay impact: 1 point
#'   \item Translation impact: 1 point
#' }
#' All variants start with a base score of 0.5. Final classifications are assigned based on the score:
#' \itemize{
#'   \item \code{empirical}: score ≥ 6.5
#'   \item \code{predictive}: 1.5 ≤ score < 6.5
#'   \item \code{weak}: score < 1.5
#' }
#'
#' @examples
#' # Assuming cisRdb is a properly formatted list as described above
#' result <- cisR_Score(cisRdb)
#' head(result$variants)
#' head(result$evidence)
#'
#' @export
cisR_Score <- function(cisRdb){
  # Extract impact
  extract_cisRdb <- function(x, id, name) {
    x <- x[x$sig %in% "sig", ]
    x <- x[, c(id, "trend")]
    names(x) <- c(id, name)
    return(x)
  }

  base <- merge(cisRdb$einteraction, cisRdb$var2erna, by="E")
  base$varE <- paste0(base$id, "_", base$E)

  empE <- extract_cisRdb(cisRdb$E, "varE", "E.impact")
  emp <- extract_cisRdb(cisRdb$G, "varEG", "E2G.impact")
  emp1 <- extract_cisRdb(cisRdb$e2TE, "pair", "TE.impact")
  emp2 <- extract_cisRdb(cisRdb$e2TE, "pair", "decay.impact")

  ## Merge
  base <- merge(base, empE, by="varE", all=T)
  base$varEG <- paste0(base$varE, "_", base$G)
  base <- merge(base, emp, by="varEG", all=T)
  dete_impact <- merge(emp2, emp1, by="pair", all=T)
  base <- merge(base, dete_impact, by="pair", all=T)

  #
  sp_id <- as.data.frame(do.call(rbind,
                                 strsplit(base$id, split = "[:_]",
                                          fixed = FALSE)))
  colnames(sp_id) <- c("chr", "start", "ref", "alt")
  sp_id$id <- base$id
  base <- base[, !names(base) %in% c("pair", "varEG", "varE", "id")]
  base <- cbind(sp_id, base)

  base$score <- 0.5
  base$score <- base$score +
    ifelse(!is.na(base$E.impact), 5, 0) +
    ifelse(!is.na(base$E2G.impact), 2.5, 0) +
    ifelse(!is.na(base$decay.impact),        1, 0) +
    ifelse(!is.na(base$TE.impact),     1, 0)


  ###
  res <- aggregate(score ~ chr + start + ref + alt + id, data = base, FUN = max)
  res$s.class <- ifelse(res$score >= 6.5, "empirical",
                        ifelse(res$score >= 1.5, "predictive", "weak") )

  return(list("variants"=res, "evidence"=base))
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
#'   @param t False discovery rate threshold for filtering G results (default: `0.1`).
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

  # Helper
  pp_out <- function(df, is, fdrt = 0.05) {
    df <- do.call(rbind, df)
    df$sig <- ifelse(round(df$local.FDR, 3) <= fdrt, "sig", "")
    df$trend <- ifelse(df$logFC > 0, "Increasing", "Decreasing")
    df$is <- is
    rownames(df) <- NULL
    return(df)
  }

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
    fit$n.var <- length(varSamples)

    return(fit)
  })
  ## The FDR here is the same as p.value due to single E tested
  bout <- pp_out(bout, "E", 0.05)
  bout$varE <- paste0(bout$var, "_", bout$ID)
  sigvar <- bout[bout$sig %in% "sig", ][["var"]]

  ## Analyse linked G
  bout2 <- lapply(sigvar, function(v){
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
    # This is valid because 1:1 E to var mapping
    fit$E <- er
    fit$n.var <- length(varSamples)

    return(fit)
  })
  bout2 <- pp_out(bout2, "G", 0.1)
  bout2$varEG <- paste0(bout2$var, "_", bout2$E, "_", bout2$ID)

  return(list("E"=bout, "G"=bout2))
  }




#' Link a given eRNA to significant impact on mRNA decay or translation efficiency
#'
#' This function processes the "eDecay2mRNA" and "TE2mRNA" slots on eRNAkitDB to get their impact on target mRNA.
#' It requires the "G" object in the DB map the gene symbol of those slots to the Ensemble id.
#'
#' @param E A vector containing eRNAs of interest. eRNA IDs must match the naming in eRNAkitDB.
#' @param DB The eRNAkitDB object.
#' @param t Numeric threshold for significance on adjusted p-values (False Discovery Rate). Default is 0.1.
#'
#' @return A named list with two elements:
#' \itemize{
#'   \item \code{de}: data.frame of significant mRNA decay impact.
#'   \item \code{te}: data.frame of significant translational efficiency impact.
#' }
#'
#' @details
#' This expects the exact slot naming convention from the eRNAkitDB.
#' The output retains the columns in the eRNAkitDB slots but replaces pair for easy mapping.
#'
#' @examples
#' \dontrun{
#' results <- to_dete(E = eRNAs,
#'                   DB = eRNAkitDB,
#'                   t = 0.1)
#' head(results$de)
#' head(results$te)
#' }
#'
#' @export
to_dete <- function(E, DB, t=0.1){
  de <- DB$eDecay2mRNA[DB$eDecay2mRNA$E %in% E, ]
  te <- DB$TE2mRNA[DB$TE2mRNA$E %in% E, ]

  # Helper
  processTAB <- function(df) {
    # FDR per E
    out <- lapply(split(df, df$E), function(subdf) {
      subdf$local.FDR <- p.adjust(subdf$pvalue, method = "BH")
      return(subdf)
    })
    out <- do.call(rbind, out)
    rownames(df) <- NULL

    # Get significant
    out$sig <- ifelse(out$local.FDR <= t, "sig", "")
    out <- out[out$sig %in% "sig", ]

    out <- merge(DB$G[c("ID", "symbol")], out,
                 by.x = "symbol", by.y = "G", all.y = TRUE)

    if (nrow(out) == 0) {
      return(cbind(out, pair = character()))
    }

    out$pair <- paste0(out$E, "_", out$ID)

    return(out)
  }

  de <- processTAB(de)
  te <- processTAB(te)

 return(list("de"=de, "te"=te))
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

    #### You can now add something here if you want
    ### This will require a different logic because model for > 1bp need thinking
    df$keep <- ifelse(nchar(df$V4) == 1 & grepl("^([ACGT])(,[ACGT])*$", df$V5),
                      "yes", "")
    df <- df[df$keep %in% "yes", ]

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




