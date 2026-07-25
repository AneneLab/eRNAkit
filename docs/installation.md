← [Back to eRNAkit](../README.md)

# Installation

### eRNAkit R package
The recommended way to install the **eRNAkit** package (including **eRNAkitApp**, small embeded data and associated workflows) is:

- `devtools::install_github("AneneLab/eRNAkit")`

Most dependencies will be installed automatically. If any are missing, install them manually:

- `install.packages(c("tidyr", "dplyr", "ggplot2", "igraph", "data.table", "rintrojs", "shiny"))`
- `BiocManager::install("GenomicRanges")`

The first launch of **eRNAkitApp** will also attempt to install any missing packages automatically.

Documentation for all implemented functions is available through the package, or [here](../man).

<br><br>
### eRNAkitDB
The complete reference database (see [Core](docs/core.md)) is provided in the repository as [`eRNAkitDB.rds`](../eRNAkitDB.rds).

A high-confidence subset required by **eRNAkitApp** is bundled with the package installation.

<br><br>
### Integration into other pipelines
Standard bioinformatics file formats, including **BED**, **GTF**, and **FASTA**, can be generated directly from **eRNAkitDB** using the package's utility functions.
These utilities automatically export annotations for use with tools such as **IGV**, **HTSeq-count**, and **bedtools**.

See [Core](docs/core.md) for additional details.
