#' Pivot a wide db table to long
#'
#' Converts wide emi db table to long format for processing
#'
#' @param df A data frame in wide format.
#' @param keep_cols A character vector of column names to retain (not pivot).
#'
#' @return A data frame in long format with columns: the retained identifiers,
#'         a `variable` column for the former column names, and a `value` column for the values.
#'
#' @export
#' @importFrom tidyr pivot_longer
wide2long <- function(df, keep_cols) {
  tidyr::pivot_longer(
    df,
    cols = -all_of(keep_cols),
    names_to = "variable",
    values_to = "value"
  )
}


#' Check if all elements of a db "emi" is NULL.
#'
#' This is specific for Dr. Anene's DBs.
#' It is likely useful after subsetting the DB.
#'
#' @param db The DB object to check for NULL
#'
#' @returns TRUE if all elements are NULL or FALSE.
#' @export
check_null <- function(db) {
  res <- sapply(db, is.null)
  return(all(res))
}


#' Filter a list of data frames by column and value
#'
#' This function filters each data frame in a list by a specific value in a specified column.
#' It skips any data frames that do not contain the column.
#'
#' @param db A list of data frames.
#' @param col The name of the column to filter by (e.g., "E" or "G").
#' @param value The value to filter for.
#'
#' @return A subset of the db
#' @examples subDB(emi, "E", "en100035")
#' @export
subDB <- function(db, col, value) {
  sdb <- lapply(db, function(df) {
    if (!col %in% names(df)) return(NULL)
    out <- df[df[[col]] %in% value, ]
    if (nrow(out) == 0) return(NULL)
    return(out)
  })
}


#' Extract end coordinates from cigar string.
#'
#' Function to extract end position from cigar string.
#' The function only adds M, D and S to the count.
#'
#' @param start Start position.
#' @param strand The strand. Defualt is + standard.
#' @param cigar The cigar string.
#'
#' @returns The end coordinate.
#' @export
eeCigar <- function(start=2, strand="+", cigar="10M5D20M2I") {

  ops <- unlist(strsplit(cigar, "(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)",
                         perl = TRUE))

  tl <- 0

  #
  for (i in seq(1, length(ops), by = 2)) {
    l <- as.numeric(ops[i])
    op <- ops[i + 1]

    if (op %in% c("M", "D", "N")) {
      # Only M (Match), D (Deletion), and N (Skipped)
      tl <- tl + l
    }
  }

  # strand
  if (strand == "+") {
    endp <- start + tl - 1
  } else {
    endp <- start - tl + 1
  }

  return(max(1, endp))
}


#' Get id of max column by rows.
#'
#' Function to get ids of max columns in a vector.
#' Second max column is include for future tooling.
#' 0 values are set to NA to avoid return 0 as max.
#'
#' @param x Named vector to get max values.
#'
#' @returns Returns the max and second max ids.
#' @export
#'
#' @examples
#' max_cols(c(A = 10, B = 20, C = 30))
max_cols <- function(x) {
  x[x == 0] <- NA
  sorted <- order(x, decreasing = TRUE, na.last = TRUE)

  max1 <- ifelse(is.na(x[sorted[1]]), NA,
                 names(x)[sorted[1]])

  max2 <- ifelse(is.na(x[sorted[2]]), NA,
                 names(x)[sorted[2]])

  return(c(max1, max2))
}

#' Get id of max column by rows.
#'
#' Function to get ids of max columns in a data frame by row.
#' Second max column is include for future tooling.
#' 0 values are set to NA to avoid return 0 as max.
#'
#' @param df Data frame or matrix to get max values by row.
#'
#' @returns Returns a data frame of max and second max column ids.
#' @export
#'
#' @examples
#' max_cols(matrix(1:9, nrow = 3, ncol = 3))
max_cols_matrix <- function(df) {
  out <- apply(df, 1, function(row) {
    row[row == 0] <- NA
    sorted <- order(row, decreasing = TRUE, na.last = TRUE)

    max1 <- ifelse(is.na(row[sorted[1]]), NA,
                   colnames(df)[sorted[1]])

    max2 <- ifelse(is.na(row[sorted[2]]), NA,
                   colnames(df)[sorted[2]])
    c(max1, max2)
  })

  out <- as.data.frame(t(out))
  names(out) <- c("max1", "max2")
  return(out)
}


#' Compute median absolute deviation with constant.
#'
#' @param x A vector of values to compute.
#' @param constant The constant value to add mad equal to 0.
#' @param na.rm Option to remove "TRUE" or not "FALSE" NA values.
#'
#' @returns MAD for each value in the vector
#' @export
#'
#' @examples
#' compute_mad(c(2, 4, 7, 8,9))
compute_mad <- function(x) {
  mad <- stats::mad(x, constant = 1.4826,
                    na.rm = TRUE)

  # Add constant to zero values to avoid 0 division
  # Adding to just the 0 preserves dispersion
  ifelse(mad == 0, 1e-6, mad)
}


#' Filter matrix to remove rows with sum 0.
#'
#' @param matrix
#'
#' @returns The filtered matrix
#' @export
#'
#' @examples
#' remove_zeros(matrix(1:9, nrow = 3, ncol = 3))
remove_zeros <- function(matrix) {
  matrix[rowSums(matrix) != 0, , drop = FALSE]
}


#' Process Gene String
#'
#' Splits a comma-separated gene string into a character vector and trims any surrounding white spaces.
#'
#' @param ss A single character string with gene names separated by commas (e.g. `"gene1, gene2 , gene3"`).
#'
#' @return A character vector of cleaned gene names.
#' @export
stringS <- function(ss) {
  gS <- unlist(strsplit(ss, ","))
  trimmedS <- trimws(gS)
  return(trimmedS)
}
