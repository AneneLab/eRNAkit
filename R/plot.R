
#' Plot Expression Across Tissue or Cells
#'
#' This function generates a bar plot showing the expression levels of eRNAs across different organs or cells.
#' It expects the data from emi db.
#' The data is transformed into a long format using `pivot_longer`, excluding certain columns (`E`, `Specificity_Score`, `Is_Specific`, `expressed`).
#' The plot is ordered by expression values in descending order, and a color scale is used based on the `E` column.
#'
#' @param x A data frame containing the expression data with columns representing organs or cells and an `E` column indicating the eRNA expression values.
#' @param t A string specifying the column name to be used for the x-axis. Defaults to `"Organ"`. It can be any column name in `x` (e.g., `"Cells"`).
#'
#' @return A `ggplot` object showing the expression data as a bar plot.
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr any_of
#'
#' @examples
#' # Example 2: Plot expression across cells
#' plotOC(x2, "Cells")
#' @export
plotOC <- function(x, t="Organ"){


  return(ggplot(x, aes(x = reorder(space, -CPM), y = log10(CPM+1e-06), fill = E)) +
           geom_bar(stat = "identity", position = "dodge") +
           theme_minimal() +
           theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
                 legend.position = "bottom",
                 legend.text = element_text(size = 2),
                 legend.title = element_text(size = 5),
                 legend.key.size = unit(0.6, "lines"),
                 legend.key.width = unit(0.2, "lines"),
                 legend.key.height = unit(0.2, "lines")) +
           labs(title = paste0("Expression across ", t),
                x = t,
                y = "CPM (Log10)",
                fill = "eRNA"))

}

#' Plot the figures for localisation data.
#'
#' This function expects tables pre-processed by the winloc function.
#'
#' @param df Data frame to plot from the winloc list of tables.
#' @param t The type of plot p1, p2, or p3, matching the winloc table
#'
#' @returns GGplot2 plot
#' @import ggplot2
#' @export
plotLoc <- function(df, t=c("p1", "p2", "p3")){
  if(tolower(trimws(t))=="p1") {
    ggplot(df, aes(x = E, y = logCPM, fill = C1)) +
      geom_boxplot(outlier.shape = NA, width = 0.6) +
      scale_fill_manual(values = c("cytosol"="#66C2A5","nucleus"="#8DA0CB")) +
      facet_grid(rows = vars(RNA)) +
      labs(x = "E", y = "CPM (log10)",
           fill = "Key", title = "Expression CN") +
      theme_minimal() +
      theme(
        strip.text = element_text(size = 8, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1,
                                   size=6),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.3, "lines"),
        legend.key.width = unit(0.5, "lines"),
        legend.key.height = unit(0.5, "lines"))

  } else if(tolower(trimws(t))=="p2"){
    ggplot(df, aes(x = E, y = Cell, fill = Value)) +
      geom_tile(width = 0.5, height = 0.3) +
      scale_fill_gradient2(low = "darkgreen", high = "darkred",
                           mid = "white", midpoint = 0) +
      labs(x = "", y = "Cell",fill = "REI",
        title = "polyA - vs +") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
        panel.spacing = unit(1, "lines"),
        panel.grid = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.3, "lines"),
        legend.key.width = unit(0.5, "lines"),
        legend.key.height = unit(0.5, "lines"))

  } else if(tolower(trimws(t))=="p3"){
    ggplot(df, aes(x = E, y = interaction(C1, Cell, sep = " | "),
                         fill = Value)) +
      geom_tile(color = "white", width = 0.5, height = 0.2) +
      scale_fill_gradient2(
        low = "#4575b4", mid = "white", high = "#d73027",
        midpoint = 0, name = "REI" ) +
      labs( x = "E", y = "Fraction | Cell",
        title = "subcellular fraction  vs cytoplasm" ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        panel.grid = element_blank() ,
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 7),
        legend.key.size = unit(0.3, "lines"),
        legend.key.width = unit(0.5, "lines"),
        legend.key.height = unit(0.5, "lines"))

  } else {
    print("Give t p1, p2, p3")
  }

}


#' Generate location-based subsets from emi loc table
#'
#' This function processesthe emi db loc table.
#'
#' It returns:
#' - p1: CPM values split by RNA type (only cytosol and nucleus)
#' - p2: REI values for `cytosol` polyA+ RNA split by Cell line
#' - p3: REI values comparing subcellular cytoplasmic compartments to cytosol (total RNA), split by Enhancer (E)
#' - d:  The original filtered input dataframe
#'
#' @param subloc emi$loc subset prefered.
#'
#' @return A named list with components: `p1`, `p2`, `p3`, and `d`.
#'
#' @importFrom dplyr filter %>%
#' @export
winLoc <- function(subloc){
  # CPM cytoplasm and nucleus
  p1 <- subloc %>%
    filter(Value.is == "CPM") %>%
    filter(C1 %in% c("cytosol", "nucleus")) %>%
    filter(RNA != "total")
  p1$logCPM <- log10(p1$Value + 1e-06)

  # REI cyto ployA+ and polyA+
  p2 <- subloc %>%
    filter(Value.is == "REI") %>%
    filter(RNA == "cytosol")

  # REI sub cellular and cytoplasm
  p3 <- subloc %>%
    filter(Value.is == "REI") %>%
    filter(C2 == "cytosol") %>%
    filter(RNA == "total")

  out <- list()
  out[["p1"]] <- p1
  out[["p2"]] <- p2
  out[["p3"]] <- p3
  out[["d"]] <- subloc

  return(out)

}
