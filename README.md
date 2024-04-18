# RBET
Reference-based Batch Effect Testing (RBET) is a statistical framework aimed to guide the selection of batch effect correction (BEC) methods. It consists of two parts: (1) selecting reference genes, and (2) testing batch effect on reference genes. Reference genes can be selected either from candidate housekeeping genes or directly from data, which correspond to literature-based RBET and data-based RBET.

## Requirements
R version: >= 4.0 <br />
R packages:<br />
Seurat == 5.0.3 <br />
uwot == 0.1.16 <br />
RVAideMemoire == 0.9-83-7 <br />
dplyr == 1.1.4 <br />
Rcpp == 1.0.12 <br />

## Installation
Install directly from github (install the latest release):<br />
`devtools::install_github("zlyx26/RBET")`

## Tutorial
`Tutorial.Rmd` provides instructions on how to use RBET. The example dataset in the tutorial can be downloaded from <https://drive.google.com/drive/folders/18rdL-L8nHL3MsojmuSIUiiI2NDPldDD9?usp=sharing>.

## Citation
Xiaoyue Hu#, He Li#, Ming Chen, Junbin Qian*  and Hangjin Jiang*. RBET guides case-specific choice of integration methods
for scRNA-seq datasets from different sources. (Submitted)
