# Last DB build
# eRNAkitDB$version <- "0.2.2"
# eRNAkitDB$genome <- "GRCh38"
# eRNAkitDB$date <- "11-08-2025"
# BioML::DBsave(eRNAkitDB, "eRNAkitDB")



#' Write eRNAkit Database to Standard Bioinformatics File
#'
#' Export the core data in an eRNAkit database to a file in one of the
#' common bioinformatics formats: BED, GTF, or FASTA.
#'
#' @param db A list representing the eRNAkit database object. Must contain
#'   \code{core} data.frame with columns \code{chr}, \code{start}, and \code{end},
#'   and \code{E_sequence} data.frame with \code{label} and \code{sequence} columns.
#'   Defaults to \code{eRNAkitDB}.
#' @param output Character string specifying the output file path or name.
#'   Defaults to \code{"eRNADBv0.2.2.gtf"}.
#' @param format Character string indicating the output format, one of
#'   \code{"bed"}, \code{"gtf"}, or \code{"fasta"}. Default is \code{"gtf"}.
#' @param chr Logical; if \code{TRUE}, ensures chromosome names start with \code{"chr"},
#'   otherwise removes \code{"chr"} prefix. Default is \code{FALSE}.
#' @param mito Optional character string to rename mitochondrial chromosome
#'   labels (e.g., \code{"chrM"} to \code{"MT"}). Default is \code{NULL}.
#'
#' @details
#' This function converts the core annotation and sequence data stored in the eRNAkit
#' database to a chosen file format, writing to disk with optional metadata appended
#' as header comments (for BED/GTF) or embedded in FASTA headers.
#'
#' @return Invisibly returns \code{NULL}. Writes output to the specified file.
#'
#' @examples
#' \dontrun{
#' exportDB(db=eRNAkitDB, mito="MT") # For ensemble
#' exportDB(output="eRNADB_chrv0.2.2.gtf", chr=T, mito="MT") # for ensemble
#' exportDB(db=eRNAkitDB, output="eRNADBv0.2.2.bed", format="bed", mito="MT")
#' exportDB(db=eRNAkitDB, output="eRNADBv0.2.2.fa", format="fasta")
#' }
#'
#' @export
exportDB <- function(db=eRNAkitDB, output="eRNADBv0.2.2.gtf", format="gtf",
                          chr = F, mito = NULL){

  ##### Check db is sane
  if (!all(c("chr", "start", "end") %in% names(db$core))) {
    stop("db$core must have: 'chr', 'start', and 'end'")
  }

  if (!all(c("label", "sequence") %in% names(db$E_sequence))) {
    stop("db$core must have: 'lable' and 'sequence'")
  }


  # bed and gtf
  io1 <- function(input, file.name = "file.bed", meta=NULL) {
    options(scipen = 999)
    con <- file(file.name, open = "wt")

    if (!is.null(meta)) {
      writeLines(paste0("#!meta ", meta), con)
    }

    write.table(input, con, sep = "\t",
                col.names = FALSE,
                row.names = FALSE, quote = FALSE)

    close(con)
  }

  # fasta
  io2 <- function(input, file.name="file.fa", meta = NULL) {
    options(scipen = 999)
    input$label <- paste0(">", input$label)
    if (!is.null(meta)) {
      input$label <- paste(input$label, meta)
    }

    con <- file(file.name, open = "wt")

    for (i in seq_len(nrow(input))) {
      writeLines(input$label[i], con)
      writeLines(input$sequence[i], con)
    }

    close(con)
  }

  # Process
  df <- db$core

  #####
  df$chr <- if (chr) {
    ifelse(grepl("^chr", df$chr),
           df$chr, paste0("chr", df$chr))
  } else {
    gsub("^chr", "", df$chr)
  }

  ####
  if (!is.null(mito)) {
    df$chr <- gsub("^chrM$|^chrMT$|^MT$|^M$",
                   mito, df$chr)
  }


  out <- switch(format,
                # bed6
                bed = {
                  data.frame(
                    chr = df$chr,
                    start = df$start,
                    end = df$end,
                    name = df$E,
                    score = df$score,
                    strand = ".",
                    stringsAsFactors = FALSE
                  )
                },
                # gtf
                gtf = {
                  data.frame(
                    chr = df$chr,
                    source = paste0("eRNAkit_", db$version),
                    type = "eRNA",
                    start = df$start,
                    end = df$end,
                    score = df$score,
                    strand = ".",
                    frame = ".",
                    attribute = df$attribute,
                    stringsAsFactors = FALSE
                  )
                },
                # fasta
                fasta = {
                  fasta = {
                    as.data.frame(db$E_sequence)

                  }
                }
  )


  if (format == "fasta") {
    io2(input = out, file.name = output,
        meta = paste(db$version, db$genome, db$date, sep="_"))
  } else {
    io1(input = out, file.name = output,
        meta = paste(db$version, db$genome, db$date, sep="_"))
  }

  moutput <- normalizePath(output, mustWork = FALSE)
  message("File written to: ", moutput)

}




#' Create n.bp Genomic Windows or Retain Original Interval for eRNAkitDB$core
#'
#' Given a genomic interval, this function splits it into non-overlapping windows
#' of n base pairs if the interval length is greater than n. If the interval
#' is shorter than or equal to n, it is returned as-is.
#'
#' @param db A list representing the eRNAkit database object. Must contain
#'   \code{core} data.frame with columns \code{chr}, \code{start}, and \code{end}
#' @param size An integer indicating the size of the windows "n".
#'
#' @return A \code{data.frame} with columns \code{chrom}, \code{start}, \code{end}, and \code{name}.
#' If the interval is longer than size, the result contains multiple rows representing the windows.
#' For shorter intervals, a single row is returned.
#'
#' @export
make_window <- function(db = eRNAkitDB, size = 100) {
  df <- db$core

  # Apply to each row and combine results
  out <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    start <- df$start[i]
    end <- df$end[i]
    chrom <- df$chr[i]
    name <- df$E[i]
    len <- df$score[i]

    if (len <= size) {
      data.frame(
        chrom = chrom,
        start = start,
        end = end,
        name = name,
        iname = name,
        stringsAsFactors = FALSE
      )
    } else {
      seq_starts <- seq(start, end, by = size)
      seq_ends <- pmin(seq_starts + (size - 1), end)
      data.frame(
        chrom = chrom,
        start = seq_starts,
        end = seq_ends,
        name = paste0(name, "_", seq_along(seq_starts)),
        iname = name,
        stringsAsFactors = FALSE
      )
    }
  }))

  rownames(out) <- NULL
  return(out)
}





