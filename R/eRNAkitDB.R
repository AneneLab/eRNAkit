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








setGeneric("load.inper", function(x, evi, info, gene.col="gene_id", info.id="ID",
                                  cell.col="cell", pen="no",
                                  method.col ="info.method") standardGeneric("load.inper"))

setGeneric("integrate.inper", function(x, L1="bonfe", L2="fishe",
                                       p=0.05, penalise=F) standardGeneric("integrate.inper"))








emi <- readRDS("/Volumes/Atlas/Github/eRNAkit/inst/extdata/emi.rds")
emi_ext <- readRDS("/Volumes/Atlas/Github/eRNAkit/inst/extdata/emi_ext.rds")
emi_ext$interaction_rna_extended <- emi_ext$interaction_rna_extended[emi_ext$interaction_rna_extended$enh %in% emi$core$E, ]
emi_ext$interaction_rna_locations <- emi_ext$interaction_rna_locations[emi_ext$interaction_rna_locations$ID.x %in% emi$core$E, ]

emi_ext$E <- emi_ext$E[emi_ext$E$ID %in% emi$core$E, ]
emi_ext$D2D <- emi_ext$D2D[emi_ext$D2D$E %in% emi$core$E, ]

names(emi_ext$interaction_rna_extended) <- c("file","G", "Gs","biotype", "E","count","cpm","pair","mol","method","cell","library",
                                             "geo", "loc","treatment")

emi_ext$interaction_rna_locations <- emi_ext$interaction_rna_locations[c(1:10, 13)]
emi_ext2 <- list()
emi_ext2[["R2R_extend"]] <- emi_ext$interaction_rna_extended
  emi_ext2[["R2R_location"]] <- emi_ext$interaction_rna_locations
  emi_ext2[["D2D"]] <- emi_ext$D2D
  emi_ext2[["E"]] <- emi_ext$E
  emi_ext2[["G"]] <- emi_ext$G

  remap <- read.delim("~/Desktop/remap_eRNA_overlap.bed", header=FALSE)
  remap <- remap[c(1:4, 13)]

  split_cols <- strsplit(as.character(remap$V4), ":")

  # Create new columns by extracting from the list
  remap$TF <- sapply(split_cols, `[`, 1)
  remap$cell <- sapply(split_cols, `[`, 2)
  remap$V1 <- sub("chr", "", remap$V1)
  remap <- remap[c(1:3, 5:7)]
  names(remap) <- c("chr", "start", "end", "E", "TF", "cell")

  emi_ext2[["eRNA2TF"]] <- remap
  emi_ext2[["eRNA2RBP"]] <- ovmat
  emiH <- emi
  emiH$meta <- NULL
  eRNAkitDB <- c(emiH, emi_ext2)




decayy <- readRDS("/Volumes/Atlas/Github/eRNAkit/decayDB/core.rds")
eRNAd <- decayy[["eRNA_decay"]]
eRNAd <- eRNAd[!is.na(eRNAd$ntime), ]
eRNAd <- eRNAd[!is.na(eRNAd$HF_m1), ]
names(eRNAd)[17] <- "E"
eRNAkitDB[["eDecay"]] <- eRNAd

sigpre <- decayy[["significant_eRNA_predict_mRNA_stability"]]
sigpre <- sigpre[c(1:7, 10)]
names(sigpre) <- c("A0","High","G","trend","pvalue","E","FDR","pair" )
eRNAkitDB[["eDecay2mRNA"]] <- sigpre
DBsave(eRNAkitDB, "eRNAkitDB")


coret <- eRNAkitDB$core
coret$chr <- paste0("chr", coret$chr)
BioML::write.bed(coret[c(1, 4, 5, 10)], "eRNA_chr.bed")


# Generics
setGeneric("load.inper", function(x, evi, info, gene.col="gene_id", info.id="ID",
                                  cell.col="cell", pen="no",
                                  method.col ="info.method") standardGeneric("load.inper"))

setGeneric("integrate.inper", function(x, L1="bonfe", L2="fishe",
                                       p=0.05, penalise=F) standardGeneric("integrate.inper"))



# Methods
setMethod("load.inper", "RBPInper", function(x, evi, info, gene.col, info.id,
                                             cell.col, pen, method.col) {
  # Sanity check
  if(!is(object = x, class2 = "RBPInper")) {
    stop('Object must be of class "RBPInper"\n')
  }

  # Load the evidence
  pvals <- names(evi)[sapply(evi, function(x) is.numeric(x) || is.integer(x))]

  x@meta <- evi
  rownames(x@meta) <- paste(rownames(evi), evi[[gene.col]], sep = "_")
  x@evidence <- x@meta[pvals]

  # Load the info
  x@cells <- unique(info[[cell.col]])
  x@information <- lapply(x@cells, function(cc){
    subinfo <- info[info[[cell.col]] %in% cc, ]
  })

  names(x@information) <- x@cells
  x@id <- info.id
  x@gene_id <- gene.col
  x@cell.col <- cell.col

  # Load the method to penalize
  if(pen != "no"){
    x@method.column <- info[info[[method.col]] %in% pen, ]
    x@pen.group <- unique(x@method.column[[x@cell.col]])
    x@penalised <- lapply(x@pen.group, function(cp){
      subinfo <- x@method.column[x@method.column[[x@cell.col]] %in% cp, ]
    })
    names(x@penalised) <- x@pen.group
  }

  x
})


setMethod("integrate.inper", "RBPInper", function(x, L1, L2, p, penalise) {
  # Sanity check
  if(!is(object = x, class2 = "RBPInper")) {
    stop('Object must be of class "RBPInper"\n')
  }

  if(penalise){

    x@information <- lapply(x@information, function(tab){
      tabb <- tab[!tab[[x@id]] %in% x@method.column[[x@id]], ]
      if(nrow(tabb) >= 1){
        return(tabb)
      }else{
      }
    })

    x@information <- x@information[sapply(x@information, function(x) !is.null(x))]

    x@cells <- unique(unlist(sapply(x@information, function(tab2){
      tab2[[x@cell.col]] })))
  }

  ##
  inteed <- lapply(x@cells, function(ci){
    idinteed <- x@evidence[x@information[[ci]][[x@id]]]
    idinteed[[ci]] <- apply(idinteed, 1, RBPInper::fusion, alpha=p, type=L1)
    idinteed <- idinteed[ci]
  })

  x@L1.result <- do.call(cbind, inteed)

  # Do penalty
  if(penalise){
    print("Runing the penalised version")
    inteedP <- lapply(x@pen.group, function(ci){
      idinteedP <- x@evidence[x@penalised[[ci]][[x@id]]]

      idinteedP[[ci]] <- apply(idinteedP, 1, RBPInper::fusion, alpha=p, type=L1)
      idinteedP <- idinteedP[ci]
    })
    L1.pen <- do.call(cbind, inteedP)
    L1.pen$pen <- apply(L1.pen, 1, RBPInper::fusion,
                        alpha=p, type=L2)
    x@L1.result$pen <- L1.pen$pen
  }



  # res <- x@L1.result
  x@L2.result <- x@L1.result
  x@L2.result$global <- apply(x@L2.result, 1, RBPInper::fusion,
                              alpha=p, type=L2)

  x@L2.result <- cbind(x@meta[x@gene_id], x@L2.result["global"])
  x@L2.result$adjP <- p.adjust(x@L2.result$global, method = "BH")
  x@L2.result$call <- ifelse(x@L2.result$adjP <= p, "Hit", "")

  x
})


#' @title Integrate RBP interaction profiles
#'
#' @description Meta-analysis of RNA binding protein interaction profiles
#'
#' @param evi Data frame of p-value evidence, genes in row and data in columns
#'
#' @param info Data frame of information file, at least data ID and cell type include
#'
#' @param gene.col Column name in "evi" p-value matrix with gene IDs
#'
#' @param info.id Column name in "info" information matrix with matching sample IDs
#'
#' @param cell.col Column name in "info" information matrix with cell type name
#'
#' @param L1 Method for p-value merging at cell type level, default "bonf"
#'
#' @param L2 Method for p-value merging at global level, default "fisher"
#'
#' @param p P-value threshold for calls
#'
#' @param penalise Indicate if a method should be penalized. Sub-process for DE RNA-Seq.
#'
#' @param pen If penalise TRUE, then indicate the method to be penalized found in method.col below.
#'            This must match information in info file.
#'
#' @param method.col If penalise TRUE and pen given, then indicate the column in the info file that match pen.
#'
#' @return RBPInper object with cell type level results and global results
#'
#' @keywords RBPInper, p-value merging, RNA binding protein
#'
#' @examples
#'
#' @export
#'
rbpinper.run <- function(evi, info, gene.col="gene_id",
                         info.id="ID", cell.col="cell",
                         L1="bonfe", L2="fishe", p=0.05, penalise=F,
                         pen="no", method.col ="info.method"){

  res <- new("RBPInper")
  res <- load.inper(res, evi=evi, info=info, gene.col=gene.col,
                    info.id=info.id, cell.col=cell.col, pen=pen, method.col=method.col)

  res <- integrate.inper(res, L1=L1, L2=L2, p=p, penalise=penalise)

  print("Integration complete")
  return(res)
}
