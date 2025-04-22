#' Create n.bp Genomic Windows or Retain Original Interval
#'
#' Given a genomic interval, this function splits it into non-overlapping windows
#' of n base pairs if the interval length is greater than n. If the interval
#' is shorter than or equal to n, it is returned as-is.
#'
#' @param chrom A character or numeric value representing the chromosome name or number.
#' @param start An integer representing the start coordinate of the genomic interval (1-based).
#' @param end An integer representing the end coordinate of the genomic interval (inclusive).
#' @param name A character string providing an identifier for the interval.
#' @param size An integer indicating the size of the windows "n".
#'
#' @return A \code{data.frame} with columns \code{chrom}, \code{start}, \code{end}, and \code{name}.
#' If the interval is longer than size, the result contains multiple rows representing the windows.
#' For shorter intervals, a single row is returned.
#'
#' @export
make_windows <- function(chrom, start, end, name, size=100) {
  len <- end - start + 1
  if (len <= size) {
    return(data.frame(chrom = chrom, start = start, end = end, name = name))
  } else {
    seq_starts <- seq(start, end, by = size)
    seq_ends <- pmin(seq_starts + (size-1), end)
    return(data.frame(
      chrom = chrom,
      start = seq_starts,
      end = seq_ends,
      name = paste0(name, "_", seq_along(seq_starts)),
      iname = name
    ))
  }
}


#' Extract eRNA sequences from a FASTA file using BED coordinates
#'
#' This function extracts DNA sequences for eRNAs from a reference FASTA file based on
#' coordinates provided in a BED file. It optionally supports strand-specific
#' extraction and writes the result to a new FASTA file.
#'
#' @param b Data frame with intervals. It must "chr", "start", "end", "strand" as names. With any other mcols.
#' @param fasta Character. Path to the reference genome in FASTA format.
#'
#' @return A \code{DNAStringSet} object containing the extracted sequences.
#'
#' @import GenomicRanges
#' @export
eRNASeq <- function(b, fasta) {
  print("Needs implementing")

  return(seqs)
}
