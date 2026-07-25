# Publications

← [Back to eRNAkit](../README.md)

**Benova, N., Kuklinkova, R., Eldahshoury, M.K. and Anene, C.A., 2025. eRNAkit: Expanding the Functional Atlas of human Enhancer RNAs Beyond the Nucleus. bioRxiv, pp.2025-04.**
Enhancer RNAs (eRNAs) are a class of non-coding RNAs transcribed from active enhancers that regulate various aspects of transcription. 
Although traditionally viewed as nuclear-localised and unstable transcripts, large transcriptomic studies reveal they vary widely in localisation and biochemical properties, 
including detectable accumulation in cytoplasmic compartments. These observations suggest potential non-nuclear functions for eRNAs, yet existing databases remain focused on nuclear, 
cis-regulatory roles, limiting systematic exploration of their broader regulatory repertoire. Here we present eRNAkit, a comprehensive and accessible resource designed to address this 
gap by integrating subcellular localisation, RNA-RNA interaction and expression data for annotated eRNAs. Leveraging fractionation-based RNA-Seq datasets, eRNAkit profiles eRNA distribution 
across nuclear and cytoplasmic compartments. It incorporates gene expression data spanning major human organs and primary cell types, enabling tissue-specific analysis of eRNA function. 
Crucially, eRNAkit includes experimentally derived RNA-RNA interaction data from RIC-Seq, PARIS, and KARR-Seq, supporting exploration of trans-acting and cytoplasmic roles for eRNAs grounded 
in physical interaction evidence. eRNAkit expands current eRNA resources beyond the enhancer-promoter paradigm, offering a robust platform for dissecting non-canonical functions of eRNAs and 
advancing our understanding of their full regulatory potential in human biology. 
**Related module:** [Core](docs/core.md)


**Kuklinkova, R., Benova, N. and Anene, C.A., 2026. eRNAs modulate mRNA stability and translation efficiency to bridge transcriptional and post-transcriptional gene regulation. RNA, 32(5), pp.704-723.**
Enhancer RNAs (eRNAs) are best known for their role in transcriptional regulation, where they facilitate enhancer–promoter communication and chromatin remodeling. Yet growing evidence suggests that their 
function may extend beyond the nucleus. Here, we systematically characterize the decay kinetics of eRNAs across human cell types using time-resolved transcriptomics and kinetic modeling. While most eRNAs 
undergo canonical exponential decay, a subset displays nonlinear dynamics, suggesting context-dependent degradation mechanisms. Perturbation of core decay regulators, including components of the m6A and 
CCR4–NOT pathways, reveals that eRNA stability is modulated by a patchwork of pathways governing mRNA turnover. Integrating transcriptome-wide ribosome profiling, RNA-seq, and half-life data, we identify 
eRNAs associated with changes in mRNA stability and translation efficiency of their target protein-coding transcripts. Functional validation of one such eRNA, en4528, shows that it regulates CDKN2C and 
FAF1 mRNA independently of transcription and impacts cell migration. These findings redefine the regulatory scope of eRNAs, positioning them as active participants in post-transcriptional gene control and 
cellular behavior. 
**Related module:** [Core](docs/core.md) and [TranCi](docs/tranci.md)


**Benova, N., Kuklinkova, R., Ibenye, E., Boyne, J.R. and Anene, C.A., 2026. Reconstructing the human enhancer RNA transcriptome. bioRxiv, pp.2026-03.**
Transcript-resolved models of RNA enable functional interrogation of RNA biology by linking processing, structure, localisation, and regulatory interactions to specific RNA molecules. Across coding and 
noncoding transcriptomes, such models have been essential for defining RNA-level mechanisms relevant to physiology and disease. Enhancer RNAs (eRNAs), however, remain largely characterised without 
transcript-level definitions, and no widely adopted transcript-resolved reference exists, limiting investigation of how individual eRNAs are processed, localised, and participate in transcriptional 
regulation or their emerging post-transcriptional functions. Here, we reconstruct a transcript-resolved catalogue of human eRNAs by pan-transcriptome assembly across diverse tissues, cell types and 
compartments, defining 36,536 transcripts, including a subset with multi-exonic structure.
We show that eRNA splice junctions are reproducible features that exhibit cell-type specificity, subcellular localisation bias, and sensitivity to spliceosome perturbation. In perturbation experiments, 
eRNA splice junction usage responded to SF3B1 mutation, nuclear–cytoplasmic partitioning, and pharmacological inhibition of RNA export, demonstrating regulation across multiple layers of RNA biology. 
In head and neck squamous cell carcinoma, a subset of these junctions showed altered usage between tumour and matched normal tissue, indicating that processing varies in disease contexts. Across three 
validation contexts, nearly one-fifth of reconstructed junctions were detectable, with some showing regulated usage, supporting biological reproducibility.
**Related module:** [Transcript models](docs/transcript.md)


**Benova, N., Kuklinkova, R., Haigh, J.L., Boyne, J.R. and Anene, C.A., 2025. Prioritising Functional Noncoding Variants via eRNA Post-transcriptional Interaction Maps in Human Samples. bioRxiv, pp.2025-07.**
Noncoding variants and mutations outnumber their coding counterparts but remain challenging to interpret functionally. We present TranCi, a method that prioritises human genetic variations by integrating enhancer RNA (eRNA) 
expression with eRNA-mRNA interactome maps. By linking variant-associated changes in eRNA to downstream gene regulation, TranCi captures functional effects missed by sequence-based or chromatin-centric approaches. 
In esophageal squamous cell cancer, TranCi identifies noncoding mutations with roles in disease initiation and progression. A personalised mode enables analysis at single-patient resolution, uncovering potential individual-specific 
regulatory variants. TranCi thus provides a mechanistic framework for interpreting noncoding variations and uniquely identifies their downstream targets, where standard methods often fall short. 
**Related module:** [TranCi](docs/tranci.md)


**Kuklinkova, R., Benova, N., Kohli, J., Boyne, J.R., Roberts, W. and Anene, C.A., 2026. Stress-responsive enhancer RNAs couple chromatin reprogramming to post-transcriptional control of senescence. bioRxiv, pp.2026-03.**
Cellular senescence is accompanied by extensive epigenomic reprogramming leading to changes in enhancer RNA levels, yet how enhancer activity is translated into functional RNA-level regulation remains unclear. Here we investigate 
how enhancer reprogramming during senescence impacts functional RNA-level regulation by eRNAs.
By integrating time-resolved transcriptomic analyses across multiple primary human cell types, we identify a set of recurrently dysregulated senescence-associated enhancer RNAs (SAeRs). We focus on one of these transcripts, EN526, 
which is reproducibly repressed during senescence while its locus remains broadly stable across cell states. EN526 eRNA exhibits cytoplasmic localisation and extensive eRNA-mRNA interactions, and cytoplasmic depletion of EN526 
recapitulates its senescence-associated loss and alters the stability and translation of the cell-cycle regulator CDKN2C. EN526 perturbation further mediates stress responses, cellular survival, and extracellular remodelling 
associated with the senescence phenotype.
Together, these findings show that SAeRs changes accompanying enhancer reprogramming in senescence are not merely passive events but can act as functional intermediates linking enhancer dynamics to post-transcriptional regulatory 
networks that phenocopy key senescence-associated cellular features. Extending this model, genetic associations at the EN526 locus further connect this regulatory axis to age-related traits and circulating protein phenotypes, 
supporting its broader relevance to human ageing and disease.
**Related module:** [Stress.eDB](docs/stress.md)







