# LIST OF REQUIRED PACKAGES -----------------------------------------------
required_packages <- c(
  "GenomicRanges", "igraph", "ggplot2", "data.table", "rintrojs", "shiny",
  "shinyBS", "shinycssloaders", "shinydashboard", "shinyjs",
  "shinyWidgets", "tidyverse", "DT"
)

# Install missing packages
install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Installing missing package: %s", pkg))
      if (pkg %in% rownames(installed.packages()) == FALSE) {
        if (pkg %in% c("GenomicRanges")) {
          BiocManager::install(pkg, ask = FALSE)
        } else {
          install.packages(pkg)
        }
      }
    }
  }
}

# Install if needed
install_if_missing(required_packages)

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))
