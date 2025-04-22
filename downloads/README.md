# eRNA_v1 Annotation Files – Download Description

This directory contains processed files from the **eRNAkit** project, including genomic coordinates, sequence data, and annotation files. Below is a summary of each file and how it was generated.

---

## 📁 File Descriptions

| File Name                 | Description |
|--------------------------|-------------|
| `eRNA_v1.gtf` **core**   | GTF-formatted file with full eRNA annotations. Can be used with HTSeq-count and IGV. |
| `eRNA_v1_chr.gtf`        | Chr prefix (e.g. chr1) version of the `eRNA_v1.gtf` file. |
| `eRNA_v1_seq.fa`         | FASTA file of eRNA sequences based on the coordinates in `eRNA_v1.bed`. Useful for creating primers, siRNAs etc. |
| `eRNA_v1_seq.tab`        | Tab-delimited version of `eRNA_v1_seq.fa` with eRNA names and sequences. Useful for downstream sequence analysis in R or Python. |
| `eRNA_v1.bed`            | BED-4 formatted eRNA genomic coordinates. useful for interval operations and works well with BEDtools |
| `eRNA_v1_chr.bed`        | Chr prefix (e.g. chr1) version of the `eRNA_v1.bed` file. |
| `eRNA_v1_ext5000.bed`    | Genomics coordinates of eRNA extended by ±5000bp from their ends. Useful for find Proximal target genes |
| `eRNA_v1_chr_ext5000.bed`| Chr prefix (e.g. chr1) version of the `eRNA_v1_ext5000.bed` file.|
| `eRNA_v1_w100bp.bed`     | Genomic coordinates of eRNA windowed by 100bp. Reproduce by `eRNAkit::make_window()`. |
| `eRNA_v1_chr_w100bp.bed` | Chr prefix (e.g. chr1) version of the `eRNA_v1_w100bp.bed` file. |
| `R2R_v1_exact.bed`       | Exact intervals of high-confidence eRNA–mRNA interactions, based on KARR-Seq, RIC-Seq and PARIS data. |

---

## 🧬 Notes on Generation

- All files are gziped to enable distribution. MD5 is included for checking what to expect after download.
- FASTA and `.tab` sequence files were derived using `bedtools getfasta -fi genome.fa -bed eRNA_v1.bed -fo eRNA_v1_seq.fa -name` or option `-tab`
- Sequence naming conventions match `eRNA_v1.gtf` for consistency.
- Mitochondrial Chromosome is included in the annotation files above. However, they are masked in the db "see below" due to nuclear to cytoplasmic "i.e. counts = 0".
- Sequence for `R2R_v1_exact.bed` can be extracted using `bedtools getfasta -fi genome.fa -bed R2R_v1_exact.bed -fo R2R_v1_exact.tab -name -tab` (~1.6TB)

---

## 📁 eRNAkit db `emi.rds and emi_ext.rds`
These database files and included with the installation of eRNAkit R package.
You can also obtain the db files without installation by going to the `inst/extdata` from the root directory of the repo.


## 📌 Recommended Tools for Viewing

- BED/GTF: [IGV](https://software.broadinstitute.org/software/igv/), UCSC Genome Browser
- FASTA: `less`, `seqkit`, or Bioconductor in R
- `.tab`: Any spreadsheet editor or R/Python for parsing

---

## 🧾 License

These files are part of the eRNAkit and distributed under the GPL-3 license.

---
