#' Lunch eRNAkit Web application.
#'
#' This lunches the eRNAkit UI for interactive exploration of the emi db.
#'
#' @returns The Web UI
#' @export
eRNkitApp <- function(){
  shiny::runApp(system.file("app", "app.R", package = "eRNAkit"),
                launch.browser = T)
}


