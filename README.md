# aRNA
This repository is inteded to host code scripts and files needed to reproduce the analyses and plots that form part of the paper 'Ancient RNA expression profiles of the extinct wolly mammoth'.

&nbsp;
&nbsp;
&nbsp;

## Figures and Tables
R scripts and files able to reproduce main and supplementary figures that form part of the paper are available in the `Figure_RCode` and `Figure_Tables` folders. 

To generate all figures by running the corresponding R code, please download all files and R scripts to a local folder and modify the code to include the `input` and `output` paths where files are located and where figure plots will be saved once created. 

Please be aware to install all the needed R packages detailed at the beginning of each corresponding R script before running it.

&nbsp;

## Meta analyses
Metagenomics and metatranscriptomics interactive `.tml` files are available in the `Meta` folder.

&nbsp;

## Ancient RNA pipeline
Before processing aRNA sequencing data, users need to generate a series of annotation files from their species of interest that are required to perform alingment and gene quantification. 

Users can find `bash` scripts in the `Annotation` folder. The only file needed is a `gtf` and/or a `gff` file detailing the available gene annotation for a given species of interest. In this study, we have implemented the aRNA pipeline for both the Asian elephant ([mEleMax1]) and the woolly mammoth ([Sandoval-Velasco et al. 2024]) genome assemblies.

The gene anotation files in both `gff` and `gtf` format for the [mEleMax1] assembly are needed to run the `Script_generate_Annotation_Asian-Elephant.sh` script.

The gene annotation file in `gtf` format for the woolly mammoth assembly ([Sandoval-Velasco et al. 2024]) needed to run the `Script_generate_Annotation_Woolly-Mammoth.sh` script has been generated as part of the works of this study.

These files are the starting point from which a series of custom annotation files are created to be later used during alingment and quantification analyses. Due to their size, they have been deposited in a Figshare online repository. Pllease download them by cliking in the following [link].

Jointly with the `gtf` and `gff` annotation files, you will be able to download already processed bowtie2 indexes for both the Asian elephant and the woolly mammoth assemblies, as well as their full genomic sequence data in `.fasta` format and aRNA sequencing data of one of thw wolly mammoths used in our study in `.fastq` format.













[mEleMax1]:https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_024166365.1/
[Sandoval-Velasco et al. 2024]: https://www.cell.com/cell/fulltext/S0092-8674(24)00642-1?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867424006421%3Fshowall%3Dtrue
[link]: https://figshare.com/account/projects/229998/articles/29590415







