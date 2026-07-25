# Transcript Models

← [Back to eRNAkit](../README.md)

## Overview

Enhancer RNAs (eRNAs) have traditionally been represented as genomic coordinates rather than defined RNA molecules. This limits the ability to investigate how individual eRNAs are processed, localised, structured or functionally regulated.

The Transcript Models module provides a transcript-resolved annotation of human enhancer RNAs generated through pan-transcriptome reconstruction of stranded RNA-seq datasets across diverse human tissues, cell types and subcellular compartments.

Using experimentally defined enhancer coordinates as a reference framework, RNA-seq reads were restricted to enhancer-associated regions and assembled to generate transcript models. The resulting catalogue resolves enhancer transcription into individual RNA molecules with defined boundaries, exon structures and strand information.

The resource provides:

- Transcript-resolved eRNA GTF annotation.
- Exon-exon junction BED files.
- Cross-reference between transcript models and original eRNAkit enhancer coordinates.
- Strand-resolved models for bidirectional enhancer transcription.
- Annotation of transcript structures including multi-exonic eRNA isoforms.

These models enable RNA-centric analysis of eRNA expression, processing, localisation and molecular regulation.

---

## Resource Files

| File | Description |
|---|---|
| `eRNAkit_transcripts.gtf` | Transcript-resolved eRNA annotation containing transcript, exon and locus information. |
| `eRNAkit_junctions.bed` | BED-format exon-exon junction annotation for splice junction quantification. |
| `eRNAkit_transcripts.fa` | FASTA sequences corresponding to reconstructed transcript models (when available). |

---

## GTF Annotation

The GTF file follows standard GTF formatting, with transcript-specific information stored in column 9 attributes.

Example: chr1  eRNAkit  transcript  51001294  51005245  .  +  .  gene_id “EN4528”; transcript_id “EN4528.1”; eRNAkit “en4528”;


### Column 9 Attributes

| Attribute | Description |
|---|---|
| `gene_id` | Identifier assigned to the reconstructed eRNA locus. Multiple transcript isoforms from the same enhancer share the same gene-level identifier. |
| `transcript_id` | Unique identifier for each reconstructed transcript model. |
| `eRNAkit` | Link between reconstructed transcript models and the original eRNAkit enhancer coordinate annotation. Allows integration with eRNAkitDB expression, localisation and regulatory datasets. |
| `gene_type` | Transcript classification assigned by the assembly pipeline. |
| `transcript_type` | Transcript classification assigned to the reconstructed model. |
| `exon_number` | Position of each exon within the transcript model. |
| `reference_id` | Identifier linking transcript models to external or originating reference annotations where applicable. |

---

## Transcript Reconstruction Strategy

Transcript models were generated using:

1. Curated enhancer coordinates from experimentally defined human enhancer resources.
2. Stranded RNA-seq datasets representing diverse tissues, cell types and subcellular fractions.
3. Read selection restricted to enhancer-associated genomic regions.
4. Pooled transcript assembly to maximise recovery of low-abundance eRNA molecules.
5. Filtering to remove gene-associated transcription and unlikely transcript models.

The final catalogue provides:

- 18,740 reconstructed eRNA loci.
- 19,484 transcript models.
- 1,538 reconstructed splice junctions.
- Multi-exonic transcript structures representing processed eRNA species.

---

## Biological Applications

The transcript models enable:

| Application | Description |
|---|---|
| Isoform-aware quantification | Quantify individual eRNA transcript models rather than enhancer-level activity. |
| Splicing analysis | Investigate exon connectivity and regulated splice-junction usage. |
| RNA localisation studies | Analyse transcript-specific behaviour across cellular compartments. |
| RNA interaction studies | Map RNA-binding proteins and RNA-RNA interactions to defined eRNA molecules. |
| Experimental design | Design transcript-specific targeting strategies including primers, siRNAs and cloning constructs. |

---

## Integration with eRNAkitDB

Transcript models are integrated into **eRNAkitDB** through the `eRNAkit` identifier.

This enables direct connection between transcript-resolved models and:

- eRNA expression profiles.
- Subcellular localisation datasets.
- RNA-RNA interaction maps.
- RNA stability measurements.
- Translation-associated datasets.

---

## Citation

For the transcript-resolved enhancer RNA catalogue, please cite:

**Reconstructing the human enhancer RNA transcriptome**

Natalia Benova, Rene Kuklinkova, Emanuela Ibenye, James R. Boyne, Chinedu A. Anene

DOI: [add DOI]
