# **Core**       ← [Back to eRNAkit](../README.md)

## Overview

The Core is the foundational framework of eRNAkit, providing a unified
platform for studying human enhancer RNAs (eRNAs) as functional
regulatory molecules. It integrates multiple orthogonal genomic and
transcriptomic resources into **eRNAkitDB**, a harmonised annotation
database centred on experimentally defined enhancer coordinates. The
database combines complementary evidence spanning expression,
subcellular localisation, RNA–RNA interactions, translation, RNA
stability and regulatory context to support systems-level, RNA-centric
investigation of eRNA biology. Alongside eRNAkitDB, the Core module
includes the interactive **eRNAkitApp** for data exploration and a
comprehensive collection of analytical **Workflows** supporting
bioinformatic analysis, data extraction and wet-lab experimental design.

Although several annotations are provided as harmonised reference
resources, the accompanying workflows expose the intermediate datasets
and processing steps used to generate them. This allows users to
reproduce, customise or extend individual analyses to meet the
requirements of specific biological questions or analytical pipelines.



## 📦 eRNAkitDB

The eRNAkitDB object `eRNAkitDB.rds` contains harmonized eRNA
annotations and interaction datasets to support integrative regulatory
analysis. It includes post-transcriptional regulatory information across
human cell types, experimental techniques and molecules. You can
download the eRNAkitDB.rds from the root of the repository.

Below is a summary of each slot.

#### Foundation

| Slot Name | Description |
|---------------------------------|---------------------------------------|
| `core` | Consensus eRNA annotation (hg38), with genomic coordinates and identifiers. Serves as the central reference set. |
| `E` | BED formatted companion for the `core` annotation. |
| `G` | Gene-level BED formatted annotation for potential mRNA targets. |
| `EOC` | Tissue- and organ-level eRNA expression quantification and specificty. |
| `loc` | Sub cellular localisation expression and index (cytoplasm, nucleus, fractions) for eRNAs across cell types. |
| `PAcy` | Polysome and Ribosome profile infromation for eRNAs. |
| `w100bp_n25` | ±100bp window expression of eRNAs across 25 samples of polyA+/-. Serves as valuable resource for designing siRNAs, primers and cloning constructs. |
| `E_sequence` | FASTA-format eRNA sequences for oligo designs and motif/structure analysis. |

### Interaction Data

| Slot Name | Description |
|---------------------------------|---------------------------------------|
| `R2R` | RNA–RNA interactions between eRNAs and mRNAs (from KARR-Seq, RIC-seq and PARIS). |
| `R2R_location` | Exact location of the R2R interactions for feature mapping and motif/structure analysis. |
| `D2D` | 3D genome interactions linking eRNAs to genes (via Hi-C, derives from `TADS`). |
| `eRNA2TF` | eRNA–TF associations from ChIP-seq overlaps (via remap2022). |
| `eRNA2RBP` | eRNA-RNA-binding protein interactions from eCLIP (via Encode datasets). |

### Decay and Translation

| Slot Name | Description |
|---------------------------------|---------------------------------------|
| `eDecay` | eRNA RNA decay kinetics from time-resolved transcriptions (via Actinomycin D, set92). |
| `eDecay2mRNA` | Impact of eRNA half-lives on target mRNA stability (derives from `eDecay` and `R2R`). |
| `TE2mRNA` | Impact of eRNA expression on target mRNA translation efficiency (via RPFdb). |

### Sequence and Genomic Context

| Slot Name | Description |
|---------------------------------|---------------------------------------|
| `TADS` | Genome wide topologically associating domains (TADS, via Encode). Serves as template for upstream regulations. |
| `Regulome` | Genome wide regulatory element status across multiple cell types (via Ensemble regulatory build). |

**A column mapping for the information in these slots can be found inside eRNAkitDB.meta**



## 🔁 Workflows

Several upstream and downstream processing functions are available in
the `eRNAkit` package to support common eRNA analysis workflows.
Examples include:

| Function | Description |
|---------------------|---------------------------------------------------|
| `R2R()` | Call RNA-RNA interactions in STAR alignment junction files. |
| `fit_decay2()` | Model decay kinetics of RNA expression in decay experiments. |
| `eRNkitApp()` | Launch the interactive `eRNAkit` Shiny app. |
| Wet lab support | Utility functions to extract information for desiging wet lab experiments such as siRNA, primers and cloning. |
| DB operations | Utility functions to extract standard bioinformatics formats including, gtf, BED4, windows etc |

 *...and many more.*

📘 Full documentation and usage examples is available through the
standard R manual files.



## 🖥️ eRNAkitApp

The eRNAkit Shiny app provides an interactive dashboard for exploring
the eRNAkitDB. The UI provides multiple search options, including by
eRNA, gene and coordinates.

#### Main programms

The dashboard features buttons to switch views and analyses, including:

| Label        | Content                                                 |
|--------------|---------------------------------------------------------|
| LOCALISATION | eRNA localisation information.                          |
| TARGET       | eRNA-target interactions.                               |
| EXPRESSION   | eRNA expression across primary cells and major tissues. |
| RIBOSOME     | Ribosome and polysome association.                      |

### Notes

We will add new UI features as we develop the resource. Please, get in
touch if you will want to see something implemented. Currently, the UI
is most suitable for exploring high confidence interactions for a few
eRNAs. For large scale analysis, we recommend using the eRNAkitDB
directly which contains extended datasets and annotations.



## Citation

**For the Foundation components of the package, please cite:** 
<br> eRNAkit: Expanding the Functional Atlas of human Enhancer RNAs Beyond
the Nucleus <br> Natalia Benova, Rene Kuklinkova, Mahmoud K Eldahshoury,
Chinedu Anthony Anene <br> doi:
<https://doi.org/10.1101/2025.04.25.650683>

**For eRNA's link to target mRNA stability and translation efficiency, please cite:** 
<br> eRNAs Modulate mRNA Stability and Translation
Efficiency to Bridge Transcriptional and Post-transcriptional Gene
Regulation <br> Rene Kuklinkova, Natalia Benova, Chinedu Anthony Anene <br> 
doi: <https://doi.org/10.1101/2025.07.06.663389>
