#' Load a database from the eRNAkit.
#'
#' Function to load DB .rds DB into environment.
#'
#' @param name Name of the db set to load (without .rds)
#' @return The loaded R object
#' @export
#'
#' @examples
#' ribo <- loadDB("ribo")
loadDB <- function(name) {
  path <- system.file("extdata", paste0(name, ".rds"),
                      package = "eRNAkit")

  if (path == "") stop("Database not found : ", name)
  readRDS(path)
}




# Class
base = "data.frame"

# Set the class
setClass("eRNAkitDB", slots = list(core = base, E = base, G = base,
                                  EOC = base, loc = base,
                                  R2R = base,
                                  R2R_extend = base,
                                  R2R_location = base,
                                  D2D = base,
                                  PAcy = base,
                                  eRNA2TF = base,
                                  eRNA2RBP = base,
                                  eDecay = base,
                                  eDecay2mRNA = base,
                                  eTE = base,
                                  eTE2mRNA = base))

# Generics
setGeneric("getCoordinates", function(object, ...) standardGeneric("getCoordinates"))

setMethod("getCoordinates", "eRNAkitDB", function(object,
                                                  format = c("bed", "gtf", "fasta"),
                                                  chr.prefix = TRUE,
                                                  rename.mito = NULL,
                                                  genome = NULL) {
  format <- match.arg(format)
  df <- object@core

  required_cols <- c("chr", "start", "end")
  if (!all(required_cols %in% names(df))) {
    stop("core slot must contain columns: 'chr', 'start', and 'end'")
  }

  # Add chr prefix if needed
  df$chr <- if (chr.prefix) {
    ifelse(grepl("^chr", df$chr), df$chr, paste0("chr", df$chr))
  } else {
    gsub("^chr", "", df$chr)
  }

  # Rename mitochondrial chromosome if requested
  if (!is.null(rename.mito)) {
    df$chr <- gsub("^chrM$|^chrMT$|^MT$|^M$", rename.mito, df$chr)
  }

  row_ids <- if (!is.null(rownames(df))) rownames(df) else paste0("eRNA_", seq_len(nrow(df)))

  out <- switch(format,
                bed = {
                  data.frame(
                    chr = df$chr,
                    start = df$start,
                    end = df$end,
                    name = row_ids,
                    score = 0,
                    strand = if ("strand" %in% names(df)) df$strand else ".",
                    stringsAsFactors = FALSE
                  )
                },
                gtf = {
                  data.frame(
                    chr = df$chr,
                    source = if ("source" %in% names(df)) df$source else "eRNAkit",
                    type = if ("type" %in% names(df)) df$type else "exon",
                    start = df$start,
                    end = df$end,
                    score = ".",
                    strand = if ("strand" %in% names(df)) df$strand else ".",
                    frame = ".",
                    attribute = paste0("gene_id \"", row_ids, "\";"),
                    stringsAsFactors = FALSE
                  )
                },
                fasta = {
                  if (is.null(genome)) stop("To export FASTA, provide a valid BSgenome object.")
                  requireNamespace("GenomicRanges", quietly = TRUE)
                  requireNamespace("Biostrings", quietly = TRUE)
                  gr <- GenomicRanges::GRanges(
                    seqnames = df$chr,
                    ranges = IRanges::IRanges(start = df$start + 1, end = df$end),
                    strand = if ("strand" %in% names(df)) df$strand else "*"
                  )
                  seqs <- Biostrings::getSeq(genome, gr)
                  names(seqs) <- row_ids
                  seqs
                }
  )

  return(out)
})


