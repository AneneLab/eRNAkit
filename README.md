# eRNAkit

**Expanding the Functional Atlas of human Enhancer RNAs Beyond the Nucleus**

**eRNAkit** is the first comprehensive resource specifically designed to investigate the cytoplasmic functions of enhancer RNAs (eRNAs).
It integrates a suite of analytical workflows for processing and analysing eRNA-related data, facilitating tasks such as data
manipulation, statistical analysis, and visualisation.  

To assist users, eRNAkit provides a transcriptome-wide map of:

- eRNA subcellular localisation
- eRNA-mRNA interactions
- Expression speificity across major tissues and primary cell types
- eRNA-ribosome associations

These insights help elucidate the functional roles of eRNAs beyond the nucleus.
We are currently working on expanding the analytics workflows and database to include - 
additional utility functions and profiles of eRNA under various stress conditions "STRESS eRNA".

### Installation
The recommended way to install eRNAkit and its associated database is by running:

- `devtools::install_github("username/eRNAkit")`

Most dependencies should install automatically. If not, install them manually using: 

- `install.packages(c("tidyr", "dplyr", "ggplot2", "igraph", "data.table", "rintrojs", "shiny"))`
- `BiocManager::install("GenomicRanges")`

The first run of the eRNAkitApp should also install any missing packages automatically.

### **Downloads**
Standard bioinformatics files such as .bed, .gtf and .fa are available in the `downloads` folder in the root directory.
Combinations of files needed for integration into other workflows are also provided.

The emi.rds database file includes a `core` table that can be used to recreate key resources.

For windowed analysis, 100bp windows of sufficiently long eRNAs are included in the `donwloads` folder. 
To generate windows of custom lengths, use the make_window() function in eRNAkit.

Description of the implemented functions are described in `eRNAkit_0.2.1.pdf` file in the root directory.

## Citation
**eRNAkit: Expanding the Functional Atlas of human Enhancer RNAs Beyond the Nucleus**
Natalia Benova, Rene Kuklinkova, Mahmoud Eldahshoury, Chinedu A. Anene
