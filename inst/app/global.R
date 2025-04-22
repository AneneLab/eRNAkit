# DATA TRANSFORMATION AND NEW VARIABLES -----------------------------------
emi_path <- system.file("extdata", "emi.rds", package = "eRNAkit")
emi <- readRDS(emi_path)

# HELP & INTRO DATA ---------------------------------------------------------------
db <- emi$meta

# FLUID DESIGN FUNCTION ---------------------------------------------------
fluid_design <- function(id, w, x, y, z) {
  fluidRow(
    div(
      id = id,
      column(
        width = 6,
        uiOutput(w),
        uiOutput(y)
      ),
      column(
        width = 6,
        uiOutput(x),
        uiOutput(z)
      )
    )
  )
}
