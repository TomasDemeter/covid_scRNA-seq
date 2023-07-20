# A molecular single-cell lung atlas of lethal COVID-19

This paper by Melms et al.  presents a comprehensive analysis of the cellular landscape and interactions in the lungs of 19 individuals who died of COVID-19 and 7 control individuals. The authors used single-nucleus RNA sequencing (snRNA-seq) to profile about 116,000 nuclei from snap-frozen lung tissue samples collected within hours of death. They identified substantial alterations in cellular composition, transcriptional cell states, and cell-to-cell interactions in COVID-19 lungs, revealing the molecular mechanisms of tissue damage and impaired regeneration. They also inferred protein activity and ligand-receptor interactions to identify potential therapeutic targets for severe COVID-19.

The main findings of the paper are:

- The lungs from individuals with COVID-19 were highly inflamed, with dense infiltration of aberrantly activated monocyte-derived macrophages and alveolar macrophages, but had impaired T cell responses.
- Alveolar type 2 cells adopted an inflammation-associated transient progenitor cell state and failed to undergo full transition into alveolar type 1 cells, resulting in impaired lung regeneration.
- The authors identified expansion of pathological fibroblasts that contributed to rapidly ensuing pulmonary fibrosis in COVID-19.
- Monocyte/macrophage-derived interleukin-1β and epithelial cell-derived interleukin-6 were unique features of SARS-CoV-2 infection compared to other viral and bacterial causes of pneumonia.
- The authors also observed ectopic tuft-like cells in the lung parenchyma of COVID-19 patients, which may have a role in airway inflammation and tissue regeneration.

The paper provides a valuable resource for understanding the pathophysiology of lethal COVID-19 and developing novel therapeutic strategies.

## scvi and scanpy

To reanalyse the data obtained from this paper, I have used mainly two tools: scvi and scanpy.

scvi is a probabilistic framework for single-cell omics analysis that leverages deep generative models. scvi can perform various tasks such as dimensionality reduction, batch correction, differential expression analysis, imputation, denoising, and integration of multi-modal data. scvi is implemented as a Python library that is compatible with scanpy.

scanpy is a scalable toolkit for analysing single-cell gene expression data. scanpy can perform preprocessing, visualization, clustering, trajectory inference, differential expression testing, gene set enrichment analysis, and simulation of gene regulatory networks. scanpy is also implemented as a Python library that integrates with scvi.

Using these tools, I have performed the following steps:

- I downloaded the raw count matrices from the Gene Expression Omnibus (GEO) database  and loaded them into scanpy.
- I filtered out low-quality nuclei based on the number of genes, UMIs, mitochondrial reads, and doublets.
- I integrated the individual samples using scvi's `scvi.data.setup_anndata` and `scvi.model.SCVI` functions to create a latent representation that accounts for batch effects.
- I performed dimensionality reduction using scanpy's `sc.tl.umap` function to obtain a UMAP embedding of the integrated data.
- I performed clustering using scanpy's `sc.tl.leiden` function to identify cell types based on their gene expression profiles.
- I performed differential expression analysis using scvi's `scvi.model.SCVI.differential_expression` function to compare gene expression between clusters or conditions.
- I performed gene set enrichment analysis using scanpy's `sc.tl.rank_genes_groups` function to identify enriched pathways or signatures in each cluster or condition.
- I performed Gene Ontology (GO) enrichment analysis and KEGG pathway analysis using the GO_Biological_Process_2023 and KEGG_2021_Human libraries.
- I performed ETV5 gene expression comparison between controls and COVID19 samples and determined statistical significance using Mann–Whitney U test.
- I performed gene signatures scoring using scanpy’s `sc.tl.score_genes` function. To calculate a score for each cell based on the average expression of a set of genes, normalized by the average expression of a set of control genes.

: Melms JC et al. A molecular single-cell lung atlas of lethal COVID-19. Nature. 2021 Jul;595(7865):114-119. doi: 10.1038/s41586-021-03569-1. Epub 2021 Apr 29. PMID: 33915569; PMCID: PMC8084440.

: Lopez R et al. A joint model of unpaired data from scRNA-seq and spatial transcriptomics for imputing missing gene expression measurements. Bioinformatics. 2020 Aug 1;36(Suppl_1):i308-i317. doi: 10.1093/bioinformatics/btaa460. PMID: 32657414; PMCID: PMC7355239.

: Wolf FA et al. Scanpy: large-scale single-cell gene expression data analysis. Genome Biol. 2018 Apr 6;19(1):15. doi: 10.1186/s13059-017-1382-0. PMID: 29626926; PMCID: PMC5889549.

: Melms JC et al. A molecular single-cell lung atlas of lethal COVID-19. Gene Expression Omnibus (GEO). 2021 . Available from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157103