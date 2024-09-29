## Brain Single-Cell Drug Treatment Annotation Repository

This repository contains data for annotating single-cell brain queries that were exposed to various drug treatments at **0H**, **24H**, and **72H** time points. The data is specifically designed for brain cells undergoing different drug treatments and aims to classify cell types accurately using modern computational tools.

### Reference Datasets:
- **hiPSC Data**: The hiPSC data used as a reference is sourced from the study published in **Nature Neuroscience** (PubMed ID: [37451260](https://pubmed.ncbi.nlm.nih.gov/37451260/)). The dataset can be accessed and read using the **Cell Ranger** file format with **scanpy’s `read_10x`** function.
- **Full Fetal Reference**: The fetal reference data used in this repository is derived from the study published in **Science** ([DOI: 10.1126/science.adf1226](https://www.science.org/doi/10.1126/science.adf1226)).

### Packages Used:
- **Python**:
  - `scanpy` (for reading and preprocessing single-cell RNA data)
  - `celltypist` (for cell-type annotation using pre-trained models)
  - `anndata` (for handling single-cell data formats)
- **R**:
  - `SingleR` (for cell-type annotation in R)

Make sure all packages are installed and updated to the latest version to ensure compatibility with the datasets provided.
