#! /bin/bash/



echo ""
echo "Script to preprocess ancient RNA sequencing data"
echo "            Emilio Mármol-Sánchez. 2025"

#Trimming adapters
echo ""
echo "> > > > > > Trimming adapters and removing sequences <18nt"
cutadapt -j 8 -a TGGAATTCTCGGGTGCCAAGG -e 0.1 -m 26 --trimmed-only --trim-n -o Mammoth1_Trimmed.fastq Mammoth1.fastq.gz 1> Mammoth1_trimming.log
echo ""
echo "Compressing data"
bgzip Mammoth1_Trimmed.fastq

echo "DONE..."



#Collapse (PCR deduplication)
echo ""
echo "> > > > > > PCR deduplication"
gzcat Mammoth1_Trimmed.fastq.gz | fastx_collapser -o Mammoth1_Trimmed_Collapsed.fasta

echo "DONE..."



#Remove UMIs and put in header
echo ""
echo "> > > > > > Reshaping UMIs"
awk '$1 ~ /^>/ {print $1}' Mammoth1_Trimmed_Collapsed.fasta | gsed 's/\-/_x/g' | gsed 's/>/>seq/g' > header1
awk '$1 ~ /^[A-Z]/ {print substr($1,1,4)}' Mammoth1_Trimmed_Collapsed.fasta > umi1 
awk '$1 ~ /^[A-Z]/ {print substr($1,length($1)-3,length($1))}' Mammoth1_Trimmed_Collapsed.fasta > umi2

paste umi1 umi2 | gsed 's/\t//g' > umis
paste header1 umis | gsed 's/\t/_/g' > header
seqtk trimfq -b 4 -e 4 Mammoth1_Trimmed_Collapsed.fasta | awk '{print $1}' > Mammoth1_Trimmed_Collapsed_noUMIs.fasta

rm Mammoth1_Trimmed_Collapsed.fasta

awk '$1 ~ /^[A-Z]/ {print $1}' Mammoth1_Trimmed_Collapsed_noUMIs.fasta > seqs
paste header seqs | awk '{print $1 ORS $2}' > Mammoth1_Trimmed_Collapsed.fasta

rm Mammoth1_Trimmed_Collapsed_noUMIs.fasta
rm header*
rm umi*
rm seqs

echo "DONE..."


