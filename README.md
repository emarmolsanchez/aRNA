# aRNA
This repository is inteded to host code scripts and files needed to reproduce the analyses and plots that form part of the paper 'Ancient RNA expression profiles of the extinct wolly mammoth'.

## Figures and Tables
R scripts and files able to reproduce main and supplementary figures that form part of the paper are available in the `Figure_RCode` and `Figure_Tables` folders. 

To generate all figures by running the corresponding R code, please download all files and R scripts to a local folder and modify the code to include the `input` and `output` paths where files are located and where figure plots will be saved once created. 

Please be aware to install all the needed R packages detailed at the beginning of each corresponding R script before running it.

## Meta analyses
Metagenomics and metatranscriptomics interactive `.tml` files are available in the `Meta` folder.

## Ancient RNA pipeline
Before processing aRNA sequencing data, users need to generate a series of annotation files from their species of interest that are required to perform alingment and gene quantification. 

In order to do that, users can find `bash` scripts in the `Annotation` folder. The only file needed is a `gtf` and/or a `gff` file detailing the available gene annotation for a given species of interest. In this study, we have implemented the aRNA pipeline for both the Asian elephant ([mEleMax1]) and the woolly mammoth ([Sandoval-Velasco et al. 2024]) genome assemblies.









[mEleMax1]:https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_024166365.1/
[Sandoval-Velasco et al. 2024]: https://www.cell.com/cell/fulltext/S0092-8674(24)00642-1?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867424006421%3Fshowall%3Dtrue







