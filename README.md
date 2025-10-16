# aRNA
This repository is inteded to host code scripts and files needed to reproduce the analyses and plots that form part of the paper **"Ancient RNA expression profiles of the extinct wolly mammoth"**.

&nbsp;

## Contact

emilio.marmol.sanchez@gmail.com

&nbsp;
&nbsp;

----------------------------------------------------------------------------------------------------------------------------------------------------------------------

&nbsp;
&nbsp;
&nbsp;

## Figures and Tables
R scripts and files able to reproduce main and supplementary figures that form part of the paper are available in the `Figure_RCode` and `Figure_Tables` folders. 

To generate all figures by running the corresponding R code, please download all files and R scripts to a local folder and modify the code to include the `input` and `output` paths where files are located and where figure plots will be saved once created. 

Please be aware to install all the needed R packages detailed at the beginning of each corresponding R script before running it.

&nbsp;

## Meta analyses
Metagenomics and metatranscriptomics interactive `.html` files are available in the `Meta` folder.

&nbsp;

## Ancient RNA pipeline
Before processing aRNA sequencing data, users need to generate a series of annotation files from their species of interest that are required to perform alingment and gene quantification. 

Users can find `bash` scripts in the `Annotation` folder. The only file needed is a `gtf` and/or a `gff` file detailing the available gene annotation for a given species of interest. Optionally, in case mitochondrial genes are not included in the annotation, users shoudl supply it in `bed` format. In this study, we have implemented the aRNA pipeline for both the Asian elephant ([mEleMax1]) and the woolly mammoth ([Sandoval-Velasco et al. 2024]) genome assemblies, jointly with the [DQ188829.2] mammoth mitochondrial assembly.

The gene anotation files in both `gff` and `gtf` format for the [mEleMax1] assembly are needed to run the `Script_generate_Annotation_Asian-Elephant.sh` script.

The gene annotation file in `gtf` format for the woolly mammoth assembly ([Sandoval-Velasco et al. 2024]) needed to run the `Script_generate_Annotation_Woolly-Mammoth.sh` script has been generated as part of the works of this study.

These files are the starting point from which a series of custom annotation files are created to be later used during alingment and quantification analyses. Due to their size, they have been deposited in a Figshare online repository. Pllease download them by cliking in the following [link]. Once downloaded, please uncompress the files that are gziped with `gunzip`.

Jointly with the `gtf` and `gff` annotation files, users can download already processed bowtie2 indexes for both the Asian elephant and the woolly mammoth assemblies, as well as their full genomic sequence data in `.fasta` format, mitochondrial and microRNA gene annotation in `bed` format, and aRNA sequencing data of one of the wolly mammoths used in our study in `.fastq` format.

Prior to running the `bash` scripts included in the pipeline (see `Annotation` folder), users should make sure that the following software is installed and running in their computers:

+ [gffread]
+ [ripgrep]
+ [perl]
+ [bedtools]
+ [seqkit]
+ [bowtie2]
+ [cutadapt]
+ [fastx-toolkit]
+ [sambamba]
+ [UMI-tools]

&nbsp;
&nbsp;

The scripts use the command `gsed` (GNU sed) typically installed in MacOS systems. If users are working on a Linux or other related environments, you can likely use sed instead. Please be aware to change all `gsed` calls to sed in all scripts allocated in the `Annotation` folder if you do not use GNU sed in your system.

Once all the needed software is properly installed, users can download the necessary `gtf` and `gff`, as well as the `bed` annotation file for mitochondrial genes and the `fasta` file of either Asian elephant or woolly mammoth genome assemblies. Please click in this [link] to download all needed files and store them in the same folder locally on your computer. Following this, please perform the following steps:

1. Run the `Script_generate_Annotation_Woolly-Mammoth.sh` and/or `Script_generate_Annotation_Asian-elephant.sh` scripts.
2. Run the `Script_preprocessing.sh` script.
3. Run the `Script_alignment-quantification_Woolly-Mammoth.sh` and/or `Script_alignment-quantification_Asian-elephant.sh` scripts.

&nbsp;
&nbsp;

After all the pipeline has run, users will find a number of files in their folder of choice containing gene quantification, coverage and overlaps with exon-exon and exon-intron junctions.














[mEleMax1]: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_024166365.1/
[Sandoval-Velasco et al. 2024]: https://www.cell.com/cell/fulltext/S0092-8674(24)00642-1?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867424006421%3Fshowall%3Dtrue
[DQ188829.2]: https://www.ncbi.nlm.nih.gov/nuccore/DQ188829.2/
[link]: https://figshare.com/articles/dataset/Index_and_annotation_files/29590415
[gffread]: https://github.com/gpertea/gffread
[ripgrep]: https://github.com/BurntSushi/ripgrep
[perl]: https://www.perl.org/
[bedtools]: https://github.com/arq5x/bedtools2
[seqkit]: https://github.com/shenwei356/seqkit
[bowtie2]: https://github.com/BenLangmead/bowtie2
[cutadapt]: https://github.com/marcelm/cutadapt/
[fastx-toolkit]: https://github.com/Debian/fastx-toolkit
[sambamba]: https://github.com/biod/sambamba
[UMI-tools]: https://github.com/CGATOxford/UMI-tools







