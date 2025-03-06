## Brain Single-Cell Analysis Repository

### Method Benchmark: scGPT vs. CellTypist

This repository includes a method benchmark comparing **scGPT** and **CellTypist** for the task of brain cell. The benchmark compares both tools in terms of cell-type prediction performance.

You can find the code and detailed results in the provided scripts for both `scGPT` and `CellTypist`. Refer to the `classification/method_benchmark` files for more information.

This dataset contains a benchmark comparision of methods os cell type annotation alongside with single cell data analysis.

### Reference Datasets:
- **hiPSC Data**: The hiPSC data used as a reference is sourced from the study published in **Nature Neuroscience** (PubMed ID: [37451260](https://pubmed.ncbi.nlm.nih.gov/37451260/)). The dataset can be accessed and read using the **Cell Ranger** file format with **scanpy’s `read_10x`** function.
- **Full Fetal Reference**: The fetal reference data used in this repository is derived from the study published in **Science** ([DOI: 10.1126/science.adf1226](https://www.science.org/doi/10.1126/science.adf1226)).


### Packages Used:
- **Python**:
  - `scanpy` (for reading and preprocessing single-cell RNA data)
  - `celltypist` (for cell-type annotation using pre-trained models)
  - `anndata` (for handling single-cell data formats)
  - `scGPT` (for cell-type annotation)
- **R**:
  - `SingleR` (for cell-type annotation in R)

