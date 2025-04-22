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


#' List all database from the eRNAkit.
#'
#' Function to list DB .rds included with the eRNAkit.
#'
#' @return A character list of all the dbs included
#' @export
#'
#' @examples
#' listDB()
listDB <- function() {
  path <- system.file("extdata", package = "eRNAkit")
  if (path == "") return(character(0))

  list.files(path, pattern = "\\.rds$",
             full.names = FALSE)
}

