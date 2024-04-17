# RBET
RBET (Reference-based batch effect testing) is designed to test the differences between two or more batches in single-cell data. RBET is implemented in a two-step strategy, which are (1) reference genes selection, and (2) batch effect testing. Reference genes can be selected either from candidate housekeeping genes or directly from data, which correspond to literature-based RBET and data-based RBET. 

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