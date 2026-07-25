### Installation
The recommended way to install eRNAkit and its associated database is by running:

- `devtools::install_github("AneneLab/eRNAkit")`

Most dependencies should install automatically. If not, install them manually using: 

- `install.packages(c("tidyr", "dplyr", "ggplot2", "igraph", "data.table", "rintrojs", "shiny"))`
- `BiocManager::install("GenomicRanges")`

The first run of the eRNAkitApp should also install any missing packages automatically.

### **Integration to other pipelines**
Standard bioinformatics files such as .bed, .gtf and .fa can easily be extracted using utility operations included with the package.
These utility functions can automatically parse the eRNAkitDB content into files that can be used directly with IGV, HTSeq-count, bedtools etc.

The eRNAkitDB database file includes a `core` table that can be used to recreate key resources.
For windowed analysis, use the make_window() function in the package.

Description of the implemented functions are through the package.

