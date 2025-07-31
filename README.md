# eRNAkit

**Expanding the Functional Atlas of human Enhancer RNAs Beyond the Nucleus**

**eRNAkit** is the first comprehensive resource specifically designed to investigate the post-transcriptional functions of human enhancer RNAs (eRNAs).

The resource has three layers including:
- **eRNAkitDB**
- **eRNAkitApp**
- **Workflows**

We are currently working on expanding the analytics workflows and database to include - 
additional utility functions and profiles of eRNA under various stress conditions "STRESS eRNA".


&nbsp;
&nbsp;


## 📦 eRNAkitDB
The eRNAkitDB object `eRNAkitDB.rds` contains harmonized eRNA annotations and interaction datasets to support integrative regulatory analysis.
It includes post-transcriptional regulatory information across human cell types, experimental techniques and molecules. 
Below is a summary of each slot.

#### Core Components
| Slot Name | Description |
|-----------|-------------|
| `core` | Consensus eRNA annotation (hg38), with genomic coordinates and identifiers. Serves as the central reference set. |
| `E` | BED formatted companion for the `core` annotation.|
| `G` | Gene-level BED formatted annotation for potential mRNA targets. |
| `EOC` | Tissue- and organ-level eRNA expression quantification and specificty. |
| `loc` | Sub cellular localisation expression and index (cytoplasm, nucleus, fractions) for eRNAs across cell types. |
| `PAcy` | Polysome and Ribosome profile infromation for eRNAs. |
| `w100bp_n25` | ±100bp window expression of eRNAs across 25 samples of polyA+/-. Serves as valuable resource for designing siRNAs, primers and cloning constructs. |
| `E_sequence` | FASTA-format eRNA sequences for oligo designs and motif/structure analysis. |

### Interaction Data
| Slot Name | Description |
|-----------|-------------|
| `R2R` | RNA–RNA interactions between eRNAs and mRNAs (from KARR-Seq, RIC-seq and PARIS). |
| `R2R_location` | Exact location of the R2R interactions for feature mapping and motif/structure analysis. |
| `D2D` | 3D genome interactions linking eRNAs to genes (via Hi-C, derives from `TADS`). |
| `eRNA2TF` | eRNA–TF associations from ChIP-seq overlaps (via remap2022). |
| `eRNA2RBP` | eRNA-RNA-binding protein interactions from eCLIP (via Encode datasets). |

### Decay and Translation
| Slot Name | Description |
|-----------|-------------|
| `eDecay` | eRNA RNA decay kinetics from time-resolved transcriptions (via Actinomycin D, set92). |
| `eDecay2mRNA` | Impact of eRNA half-lives on target mRNA stability (derives from `eDecay` and `R2R`). |
| `TE2mRNA` | Impact of eRNA expression on target mRNA translation efficiency (via RPFdb). |

### Sequence and Genomic Context
| Slot Name | Description |
|-----------|-------------|
| `TADS` | Genome wide topologically associating domains (TADS, via Encode). Serves as template for upstream regulations. |
| `Regulome` | Genome wide regulatory element status across multiple cell types (via Ensemble regulatory build). |

&nbsp;
**A column mapping for the information in these slots can be found inside eRNAkitDB.meta**

&nbsp;
&nbsp;


## 🔁 Workflows
Several upstream and downstream processing functions are available in the `eRNAkit` package to support common eRNA analysis workflows.
Examples include:

| Function             | Description                                                   |
|----------------------|---------------------------------------------------------------|
| `rna2rna()`         | Call RNA-RNA interactions in STAR alignment junction files.   |
| `listDB()`           | List available database components within eRNAkitDB.          |
| `fit_decay2()`       | Model decay kinetics of RNA expression in decay experiments. |
| `eRNkitApp()`        | Launch the interactive `eRNAkit` Shiny app.                   |
| Wet lab support      | Utility functions to extract information for desiging wet lab experiments such as siRNA, primers and cloning. |
| DB operations        | Utility functions to extract standard bioinformatics formats including, gtf, BED4, windows etc |

&nbsp;
*...and many more.* 

📘 For full documentation and usage examples is avilable through the standard R manual files.


&nbsp;
&nbsp;

## 🖥️ eRNAkitApp

The eRNAkit Shiny app provides an interactive dashboard for exploring the eRNAkitDB.
The UI provides multiple search options, including by eRNA, gene and coordinates.

#### Main programms
The dashboard features buttons to switch views and analyses, including:

| Label        | Content                                |
|--------------|-----------------------------------|
| LOCALISATION | eRNA localisation information. |
| TARGET       | eRNA-target interactions. |
| EXPRESSION   | eRNA expression across primary cells and major tissues. |
| RIBOSOME     | Ribosome and polysome association.                      |

### Notes
We will add new UI features as we develop the resource. 
Please, get in touch if you will want to see something implemented.
Currently, the UI is most suitable for exploring high confidence interactions for a few eRNAs. 
For large scale analysis, we recommend using the eRNAkitDB directly which contains extended datasets and annotations.

&nbsp;
&nbsp;

### Installation
The recommended way to install eRNAkit and its associated database is by running:

- `devtools::install_github("AneneLab/eRNAkit")`

Most dependencies should install automatically. If not, install them manually using: 

- `install.packages(c("tidyr", "dplyr", "ggplot2", "igraph", "data.table", "rintrojs", "shiny"))`
- `BiocManager::install("GenomicRanges")`

The first run of the eRNAkitApp should also install any missing packages automatically.

### **Integration to other pipelines**
Standard bioinformatics files such as .bed, .gtf and .fa can easily be extracted using utility operations including with package.
These utility functions can automatically parse the eRNAkitDB content into files that can be used directly with IGV, HTSeq-count, bedtools etc.

The emi.rds database file includes a `core` table that can be used to recreate key resources.
For windowed analysis, use the make_window() function in eRNAkit.

Description of the implemented functions are described in `eRNAkit_0.2.1.pdf` file in the root directory.

## Citation

&nbsp;

**For core components of the package, please cite:** <br>
eRNAkit: Expanding the Functional Atlas of human Enhancer RNAs Beyond the Nucleus <br>
Natalia Benova, Rene Kuklinkova, Mahmoud K Eldahshoury, Chinedu Anthony Anene <br>
doi: https://doi.org/10.1101/2025.04.25.650683

&nbsp;

**For eRNA's link to target mRNA stability and translation efficiency, please cite:** <br>
eRNAs Modulate mRNA Stability and Translation Efficiency to Bridge Transcriptional and Post-transcriptional Gene Regulation <br>
Rene Kuklinkova, Natalia Benova, Chinedu Anthony Anene <br>
doi: https://doi.org/10.1101/2025.07.06.663389 

&nbsp;

## Contacts
c.a.anene{at}leedsbeckett.ac.uk

## License
This software is licensed under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).

- **Free for academic and non-commercial research.**
- **Commercial use requires a license. Please contact caanenedr{@}outlook.com for details.**
