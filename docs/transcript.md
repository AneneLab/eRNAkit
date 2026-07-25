# Transcript Models

← [Back to eRNAkit](../README.md)

## Overview

Enhancer RNAs (eRNAs) have traditionally been represented as genomic coordinates rather than defined RNA molecules. This limits the ability to investigate how individual eRNAs are processed, localised, structured or functionally regulated.
The Transcript Models module provides a transcript-resolved annotation of human enhancer RNAs generated through pan-transcriptome reconstruction of stranded RNA-seq datasets across diverse human tissues, cell types and subcellular compartments.
Using experimentally defined enhancer coordinates as a reference framework, RNA-seq reads were restricted to enhancer-associated regions and assembled to generate transcript models. The resulting catalogue resolves enhancer transcription 
into individual RNA molecules with defined boundaries, exon structures and strand information.

The resource provides:

- Transcript-resolved eRNA GTF annotation.
- Exon-exon junction BED files.
- Cross-reference between transcript models and original eRNAkit enhancer coordinates.
- Strand-resolved models for bidirectional enhancer transcription.
- Annotation of transcript structures including multi-exonic eRNA isoforms.

These models enable RNA-centric analysis of eRNA biology, including regulation and function.
We built this due to the Anene's Lab interest in post-trancriptional eRNA functions.


## Standalone GTF Files
 **[GTFs](../gtf/transcript_models)**
| File | Description |
|---|---|
| `eRNATmod_nchr.gtf` | Chromosomes are listed without the "chr" prefix, so 1, 2, 3 etc. |
| `eRNATmod_withchr.gtf` | Chromosomes are listed with the "chr" prefix, so chr1, chr2, chr3 etc.. |

The GTF files is also include as a slot `transcript` inside the **[eRNAkitDB](../eRNAkitDB.rds)**


## GTF Annotation Information
The GTF follows standard GTF formatting, with transcript-specific information stored in column 9 attributes.
Example: chr1  eRNAkit  transcript  51001294  51005245  .  +  .  gene_id “EN4528”; transcript_id “EN4528.1”; eRNAkit “en4528”;


### **Column 9 Attributes**

This is a standard GTF column 9 formatting using `key "value";` attribute pairs.

| Attribute | Description |
|---|---|
| `gene_id` | Identifier assigned to the eRNA gene-level locus. Multiple transcript isoforms derived from the same reconstructed enhancer share the same gene identifier. |
| `transcript_id` | Unique identifier assigned to each reconstructed transcript model. |
| `eRNAkit` | Identifier linking reconstructed transcript models to the original eRNAkit enhancer annotation. Enables integration with eRNAkitDB expression, localisation and regulatory datasets. |
| `gene_source` | Source annotation of the feature. `eRNAkit` indicates the enhancer locus annotation, while `STRINGTIE` indicates reconstructed transcript models. |
| `gene_biotype` | Transcript classification. For eRNAkit annotations this is defined as `eRNA`. |
| `cluster` | Identifier of the merged enhancer cluster where applicable. `NA` indicates no assigned cluster. |
| `n_stranded` | Indicates whether the enhancer locus shows stranded transcription evidence and how many.|
| `embedded_stranded` | Indicates whether transcript reconstruction supports strand-specific embedding within the enhancer region. |
| `exon_number` | Position of each exon within the reconstructed transcript model. |
| `cov` | Transcript coverage estimated during transcript reconstruction. |
| `FPKM` | Fragments per kilobase of transcript per million mapped reads (FPKM) expression estimate. |
| `TPM` | Transcripts per million (TPM) expression estimate. |
---

## Transcript Reconstruction Strategy

Transcript models were generated using:

1. Curated enhancer coordinates from experimentally defined human enhancer resources.
2. Stranded RNA-seq datasets representing diverse tissues, cell types and subcellular fractions (n=121).
3. Read selection restricted to enhancer-associated genomic regions (enhancer cordinates are fixed).
4. Pooled transcript assembly to maximise recovery of low-abundance eRNA molecules.
5. Filtering to remove gene-associated transcription and unlikely transcript models.

The final catalogue provides:

- 18,740 reconstructed eRNA loci.
- 19,484 transcript models.
- 1,538 reconstructed splice junctions.
- Multi-exonic transcript structures representing processed eRNA species.


## Biological Applications

The transcript models enable:

| Application | Description |
|---|---|
| Isoform-aware quantification | Quantify individual eRNA transcript models rather than enhancer-level activity. |
| Splicing analysis | Investigate exon connectivity and regulated splice-junction usage. |
| RNA localisation studies | Analyse transcript-specific behaviour across cellular compartments. |
| RNA interaction studies | Map RNA-binding proteins and RNA-RNA interactions to defined eRNA molecules. |
| Experimental design | Design transcript-specific targeting strategies including primers, siRNAs and cloning constructs. |


## Integration with eRNAkitDB

Transcript models are integrated into **eRNAkitDB** through the `eRNAkit` identifier.

This enables direct connection between transcript-resolved models and:

- eRNA expression profiles.
- Subcellular localisation datasets.
- RNA-RNA interaction maps.
- RNA stability measurements.
- Translation-associated datasets.


## Citation <br>
For the transcript-resolved enhancer RNA catalogue, please cite: <br>
Benova, N., Kuklinkova, R., Ibenye, E., Boyne, J.R. and Anene, C.A., 2026. **Reconstructing the human enhancer RNA transcriptome.** bioRxiv, pp.2026-03.
