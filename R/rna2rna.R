#' Extract extract locations of eRNA:mRNA interactions.
#'
#' Function to extract locations of eRNA:mRNA interactions from GRanges object.
#'
#' @param chimGR Chimeric junction GRanges object. It should be made from ChimPRO output.
#' @param geneGR Gene annotation GRanges object. It should be made from annGE$gene.
#' @param enhGR Enhancer annotation GRanges object. It should be made from annGE$erna.
#'
#' @returns bedpe type table. xx.x and xx.y are the target pairs
#' @export
#'
rna2rnaR <- function(chimGR, geneGR, enhGR){
  ovgene <- findOverlaps(chimGR, geneGR)
  overna <- findOverlaps(chimGR, enhGR)

  in_ge <- cbind(as.data.frame(chimGR[queryHits(ovgene)]),
                 as.data.frame(geneGR[subjectHits(ovgene)]))
  in_er <- cbind(as.data.frame(chimGR[queryHits(overna)]),
                 as.data.frame(enhGR[subjectHits(overna)]))

  out <- merge(in_er, in_ge, by = "mcols.name")

  return(out)
}


#' Extract and count the eRNA:mRNA interactions.
#'
#' Function to extract and count eRNA:mRNA interactions from GRanges object.
#'
#' @param chimGR Chimeric junction GRanges object. It should be made from ChimPRO output.
#' @param geneGR Gene annotation GRanges object. It should be made from annGE$gene.
#' @param enhGR Enhancer annotation GRanges object. It should be made from annGE$erna.
#'
#' @returns Interactions counts. Uses IDs for naming enhancer and genes.
#' @export
#'
rna2rna <- function(chimGR, geneGR, enhGR){
  ovgene <- findOverlaps(chimGR, geneGR)
  overna <- findOverlaps(chimGR, enhGR)

  in_ge <- cbind(as.data.frame(chimGR[queryHits(ovgene)]),
                 as.data.frame(geneGR[subjectHits(ovgene)]))
  in_er <- cbind(as.data.frame(chimGR[queryHits(overna)]),
                 as.data.frame(enhGR[subjectHits(overna)]))

  get_col <- c("mcols.name", "mcols.side",
               "mcols.ID", "mcols.symbol",
               "mcols.biotype")

  in_ge <- in_ge[get_col]
  in_er <- in_er[get_col]

  out <- merge(in_er, in_ge, by = "mcols.name")

  out <- aggregate(rep(1, nrow(out)),
                   by = list(out$mcols.ID.x, out$mcols.ID.y),
                   FUN = sum)
  colnames(out) <- c("enh", "ID", "count")

  return(out)
}

#' Get genomic range object for annGE or chimericPRO.
#'
#' Make GRanges object from dataframe.
#' This function expects a specific input so make sure you have that.
#' It must have columns named "chr", "start", "end", "strand"
#'
#' @param tab
#'
#' @returns Returns GRanges
#' @export
getGRange <- function(tab=exmp){
  gobj <- GRanges(seqnames = tab$chr,
                  ranges = IRanges(start = tab$start,
                                   end = tab$end),
                  strand = tab$strand,
                  mcols = tab[!names(tab) %in% c("chr", "start", "end", "strand")])
  return(gobj)
}


#' Process chimeric junction file from star-alignment.
#'
#' Function to process chimeric junction files into a bedpe format.
#' It splits the pairs to make it easy to find overlaps.
#'
#' @param file The file path for the chimeric.junction from star.
#'
#' @returns A list with n1 and n2 pairs. The name column match the ID.
#' @export
chimericPRO <- function(file) {

  chim <- read.table(file, header = FALSE,
                     stringsAsFactors = FALSE)

  # STAR's Chimeric.out.junction format
  colnames(chim) <- c("chr1", "pos1", "strand1", "chr2", "pos2", "strand2", "jtype",
                      "repeat_left_len1", "repeat_right_len2", "read_name", "start_aln1",
                      "cigar_aln1", "start_aln2", "cigar_aln2")

  chim$start_aln1 <- as.numeric(chim$start_aln1)
  chim$start_aln2 <- as.numeric(chim$start_aln2)

  # This set to strand + because the star ouput is respect to +
  chim$end1 <- mapply(eRNAkit::eeCigar, chim$start_aln1,
                      "+", chim$cigar_aln1)

  chim$end2 <- mapply(eRNAkit::eeCigar, chim$start_aln2,
                      "+", chim$cigar_aln2)

  bedpe <- chim[c("chr1", "start_aln1", "end1", "strand1",
                  "chr2", "start_aln2", "end2", "strand2")]
  bedpe$name <- rownames(bedpe)

  out <- list()
  out[["n1"]] <- bedpe[c("chr1", "start_aln1", "end1", "strand1", "name")]
  out[["n2"]] <- bedpe[c("chr2", "start_aln2", "end2", "strand2", "name")]

  names(out$n1) <- c("chr", "start", "end", "strand", "name")
  out$n1$side <- "n1"
  names(out$n2) <- c("chr", "start", "end", "strand", "name")
  out$n2$side <- "n2"
  #
  out <- rbind(out$n1, out$n2)

  rm(chim, bedpe)

  return(out)
}


#' Find Overlaps Between emi$core and a Genomic Coordinate String
#'
#' Parses a coordinate string and finds overlapping entries in a emi$core.
#' Returns a data.frame including all metadata columns.
#'
#' @param erob A `GRanges` object to search within.
#' @param s A character string in the format "chr1:100-200" or "1:100-200".
#'
#' @return A data.frame of overlapping GRanges entries with metadata columns.
#' Returns an empty data.frame if no overlaps are found.
#'
#' @examples
#' find_overlap_df(my_gr, "chr1:2000-3000")
#'
#' @import GenomicRanges
#' @export
byRange <- function(erob, s) {
  query <- eRNAkit::parseGR(s)
  hits <- subsetByOverlaps(erob, query)
  hits <- as.data.frame(hits)

  # Return as data.frame, even if empty
  if (nrow(hits) == 0) {
    return(data.frame())
  } else {
    hits$ID <- sub("ID\\s*([^;]+);.*", "\\1",
                   hits$mcols.attribute)
    hits <- hits[c("ID", "seqnames", "start", "end")]

    return(hits)
  }
}


#' Parse a Genomic Coordinate String into a GRanges Object
#'
#' Takes a genomic coordinate string in the format `"chr1:100-200"` or `"1:100-200"`
#' (case-insensitive) and returns a `GRanges` object representing that region.
#'
#' @param c A character string specifying the genomic region.
#'   Valid formats include `"chr1:100-200"`, `"CHR1:100-200"`, or `"1:100-200"`.
#'   Whitespace is ignored. Chromosome prefix is normalized to lowercase `"chr"` format.
#'
#' @return A `GRanges` object corresponding to the parsed genomic range.
#'
#' @examples
#' parse_genomic_range("chr1:100-200")
#' parse_genomic_range("CHR2:1000-2000")
#' parse_genomic_range("3:500-1000")
#'
#' @import GenomicRanges
#' @export
parseGR <- function(c) {
  # rm spaces and case
  c <- gsub("\\s+", "", c)

  # Match chr1:2345-5678 or 1:2345-5678
  regec <- regexec("^([cC]?[hH]?[rR]?\\w+):(\\d+)-(\\d+)$", c)
  me <- regmatches(c, regec)

  if (length(me[[1]]) != 4) {
    stop("Invalid: Use chr1:23-56 or 1:23-56")
  }

  seq <- me[[1]][2]
  start <- as.integer(me[[1]][3])
  end <- as.integer(me[[1]][4])

  # Fix chr
  if (grepl("^chr", tolower(seq))) {
    nam <- sub("^chr", "", tolower(seq))
  } else {
    nam <- seq
  }

  return(GRanges(nam, IRanges(start, end)))
}
