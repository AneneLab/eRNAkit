← [Back to eRNAkit](../README.md)

# TranCi
**Mechanistic prioritisation of functional non-coding variants through enhancer RNA-mediated regulation**

TranCi is a computational framework for prioritising functional non-coding variants by integrating enhancer RNA (eRNA) expression with experimentally derived eRNA–mRNA interaction maps. 
Rather than relying solely on chromatin proximity or sequence-based annotation, TranCi identifies variants that alter eRNA activity and links them to downstream regulatory targets through RNA-mediated regulatory networks.
We refer to this as the post-transcriptional side of non-coding variant prioritisation.

The framework supports both cohort-level analysis using matched genomic and transcriptomic datasets and exploratory analysis of individual genomes using reference interaction annotations.


## Installation
TranCi is included as part of the **eRNAkit** R package and does not require a separate installation.

Please see the main installation guide:
**[Installation](docs/installation)**


## Functions

### `TranCi()`
The primary function for cohort-based non-coding variant prioritisation.
It integrates:

- Variant calls (VCF)
- eRNA expression (expression matrix)
- mRNA expression (expression matrix)
- eRNAkitDB (eRNAkitDB.rds)

to identify non-coding variants with evidence of functional eRNA-mediated regulation. 
Variants are prioritised according to integrated molecular evidence and classified into three confidence tiers:

- **Empirical** – supported by both expression changes and downstream regulatory evidence.
- **Predictive** – supported by interaction and functional annotations.
- **Weak** – limited supporting evidence.


### `TranCi_explore()`
This provides an exploratory single variant mode for analysing individual genomes or variant sets when matched RNA-seq data are unavailable.
Variants are mapped onto the eRNAkit reference database to identify known or predicted eRNA-mediated regulatory interactions, making the function suitable for:

- personal genomes
- clinical sequencing datasets
- candidate variant interpretation
- exploratory analyses


## Documentation
Complete documentation for every function is available directly within R.

```r
?TranCi
?TranCi_explore
```

or

```r
help(package = "eRNAkit")
```

## Citation

TranCi is described in:
Benova, N., Kuklinkova, R., Haigh, J.L., Boyne, J.R. and Anene, C.A., 2025. 
**Prioritising Functional Noncoding Variants via eRNA Post-transcriptional Interaction Maps in Human Samples.** bioRxiv, pp.2025-07.

TranCi is built on the integrated resources contained within **eRNAkitDB**, including harmonised human eRNA annotations, transcript-resolved eRNA models, experimentally derived RNA–RNA interaction maps, RNA stability annotations, translation efficiency measurements, chromatin interaction data and additional regulatory annotations. 
Together, these resources provide the molecular context required to connect non-coding variants with downstream regulatory consequences.
See **[core](docs/core.md)** for more.
