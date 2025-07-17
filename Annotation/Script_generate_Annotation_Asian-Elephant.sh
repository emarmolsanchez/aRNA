#! bin/bash/



echo ""
echo "Script to generate Gene annotation files from GFF"
echo "            Emilio Mármol-Sánchez. 2025"
echo ""

#Generate Transcriptome
echo ""
echo "> > > > > > Generating Transcriptome from GFF file"
gffread -g AsianElephant_MitoMammoth.fa --force-exons --gene2exon --t-adopt -Z -j AsianElephant_MitoMammoth_intronCoords.txt -w AsianElephant_MitoMammoth_Transcriptome.fa -W --tlf -v AsianElephant_MitoMammoth.gff3
rg ">" AsianElephant_MitoMammoth_Transcriptome.fa > AsianElephant_MitoMammoth_Transcriptome.txt

#Put single header in Transcriptome
awk '{print $1}' AsianElephant_MitoMammoth_Transcriptome.fa > AsianElephant_MitoMammoth_Transcriptome_SingleHeader.fa

#Linearize Transcriptome
perl -pe '/^>/ ? print "\n" : chomp' AsianElephant_MitoMammoth_Transcriptome_SingleHeader.fa | tail -n +2 > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear.fa

#Calculate transcript lengths
rg ">" AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear.fa | gsed 's/>//g' > header
rg -v ">" AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear.fa | awk '{print length($1)}' > lengths
paste header lengths | sort -k1,1 > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths.txt

rm header
rm lengths

echo "DONE..."

#Generate BED files by gene type (nuclear genes)
echo ""
echo "> > > > > > Generating BED files for nuclear genes"
#Proteincoding
cat AsianElephant_MitoMammoth.gff3 | awk '$3=="gene" {print}' | rg protein_coding | gsed 's/;/\t/g' | awk -F "\t" '{print $1,$4,$5,$12,"protein_coding",$7}' | gsed 's/gene=//g' | gsed 's/ /\t/g' | sort -f -k1,1 -k2,2 -V > AsianElephant_MitoMammoth_ProteinCoding_genes.bed
cat AsianElephant_MitoMammoth_ProteinCoding_genes.bed | awk '{print $4}' | sort > AsianElephant_MitoMammoth_ProteinCoding_genes_IDs.txt
cat AsianElephant_MitoMammoth.gtf | awk '$3=="transcript" {print}' | rg -w -f AsianElephant_MitoMammoth_ProteinCoding_genes_IDs.txt | rg -v "TRNAC|TRNAS|TRNAF" | awk '{print $10,$12}' | gsed 's/"//g' | gsed 's/;//g' | gsed 's/ /\t/g' | sort -k2,2 > AsianElephant_MitoMammoth_ProteinCoding_genes-transcript_Equivalence.txt
cat AsianElephant_MitoMammoth_ProteinCoding_genes-transcript_Equivalence.txt | awk '{print $2}' | sort | uniq > AsianElephant_MitoMammoth_ProteinCoding_genes-transcript_IDs.txt

#rRNAs
cat AsianElephant_MitoMammoth.gff3 | awk '$3=="gene" {print}' | rg rRNA | rg -v protein_coding | rg -v pseudogene | gsed 's/;/\t/g' | awk -F "\t" '{print $1,$4,$5,$12,"rRNA",$7}' | gsed 's/gene=//g' | gsed 's/ /\t/g' | sort -f -k1,1 -k2,2 -V > AsianElephant_MitoMammoth_rRNA_genes.bed
cat AsianElephant_MitoMammoth_rRNA_genes.bed | awk '{print $4}' | sort > AsianElephant_MitoMammoth_rRNA_genes_IDs.txt
cat AsianElephant_MitoMammoth.gtf | awk '$3=="transcript" {print}' | rg -w -f AsianElephant_MitoMammoth_rRNA_genes_IDs.txt | awk '{print $10,$12}' | gsed 's/"//g' | gsed 's/;//g' | gsed 's/ /\t/g'  | sort -k2,2 > AsianElephant_MitoMammoth_rRNA_genes-transcript_Equivalence.txt
cat AsianElephant_MitoMammoth_rRNA_genes-transcript_Equivalence.txt | awk '{print $2}' | sort | uniq > AsianElephant_MitoMammoth_rRNA_genes-transcript_IDs.txt

#snRNAs
cat AsianElephant_MitoMammoth.gff3 | awk '$3=="gene" {print}' | rg snRNA | rg -v protein_coding | rg -v pseudogene | gsed 's/;/\t/g' | awk -F "\t" '{print $1,$4,$5,$12,"snRNA",$7}' | gsed 's/gene=//g' | gsed 's/ /\t/g' | sort -f -k1,1 -k2,2 -V > AsianElephant_MitoMammoth_snRNA_genes.bed
cat AsianElephant_MitoMammoth_snRNA_genes.bed | awk '{print $4}' | sort > AsianElephant_MitoMammoth_snRNA_genes_IDs.txt
cat AsianElephant_MitoMammoth.gtf | awk '$3=="transcript" {print}' | rg -w -f AsianElephant_MitoMammoth_snRNA_genes_IDs.txt | awk '{print $10,$12}' | gsed 's/"//g' | gsed 's/;//g' | gsed 's/ /\t/g' | sort -k2,2 > AsianElephant_MitoMammoth_snRNA_genes-transcript_Equivalence.txt
cat AsianElephant_MitoMammoth_snRNA_genes-transcript_Equivalence.txt | awk '{print $2}' | sort | uniq > AsianElephant_MitoMammoth_snRNA_genes-transcript_IDs.txt

#snoRNAs
cat AsianElephant_MitoMammoth.gff3 | awk '$3=="gene" {print}' | rg snoRNA | rg -v protein_coding | rg -v pseudogene | gsed 's/;/\t/g' | awk -F "\t" '{print $1,$4,$5,$12,"snoRNA",$7}' | gsed 's/gene=//g' | gsed 's/ /\t/g' | sort -f -k1,1 -k2,2 -V > AsianElephant_MitoMammoth_snoRNA_genes.bed
cat AsianElephant_MitoMammoth_snoRNA_genes.bed | awk '{print $4}' | sort > AsianElephant_MitoMammoth_snoRNA_genes_IDs.txt
cat AsianElephant_MitoMammoth.gtf | awk '$3=="transcript" {print}' | rg -w -f AsianElephant_MitoMammoth_snoRNA_genes_IDs.txt | awk '{print $10,$12}' | gsed 's/"//g' | gsed 's/;//g' | gsed 's/ /\t/g' | sort -k2,2 > AsianElephant_MitoMammoth_snoRNA_genes-transcript_Equivalence.txt
cat AsianElephant_MitoMammoth_snoRNA_genes-transcript_Equivalence.txt | awk '{print $2}' | sort | uniq > AsianElephant_MitoMammoth_snoRNA_genes-transcript_IDs.txt

#lncRNAs
cat AsianElephant_MitoMammoth.gff3 | awk '$3=="gene" {print}' | rg lncRNA | rg -v protein_coding | rg -v pseudogene | gsed 's/;/\t/g' | awk -F "\t" '{print $1,$4,$5,$12,"lncRNA",$7}' | gsed 's/gene=//g' | gsed 's/ /\t/g' | sort -f -k1,1 -k2,2 -V > AsianElephant_MitoMammoth_lncRNA_genes.bed
cat AsianElephant_MitoMammoth_lncRNA_genes.bed | awk '{print $4}' | sort > AsianElephant_MitoMammoth_lncRNA_genes_IDs.txt
cat AsianElephant_MitoMammoth.gtf | awk '$3=="transcript" {print}' | rg -w -f AsianElephant_MitoMammoth_lncRNA_genes_IDs.txt | awk '{print $10,$12}' | gsed 's/"//g' | gsed 's/;//g' | gsed 's/ /\t/g' | sort -k2,2 > AsianElephant_MitoMammoth_lncRNA_genes-transcript_Equivalence.txt
cat AsianElephant_MitoMammoth_lncRNA_genes-transcript_Equivalence.txt | awk '{print $2}' | sort | uniq > AsianElephant_MitoMammoth_lncRNA_genes-transcript_IDs.txt

#tRNAs
cat AsianElephant_MitoMammoth.gff3 | awk '$3=="gene" {print}' | rg tRNA | rg -v protein_coding | rg -v pseudogene | gsed 's/;/\t/g' | awk -F "\t" '{print $1,$4,$5,$12,"tRNA",$7}' | gsed 's/gene=//g' | gsed 's/ /\t/g' | sort -f -k1,1 -k2,2 -V > AsianElephant_MitoMammoth_tRNA_genes.bed
cat AsianElephant_MitoMammoth_tRNA_genes.bed | awk '{print $4}' | sort > AsianElephant_MitoMammoth_tRNA_genes_IDs.txt
cat AsianElephant_MitoMammoth.gtf | awk '$3=="transcript" {print}' | rg -w -f AsianElephant_MitoMammoth_tRNA_genes_IDs.txt | awk '{print $10,$12}' | gsed 's/"//g' | gsed 's/;//g' | gsed 's/ /\t/g' | sort -k2,2 > AsianElephant_MitoMammoth_tRNA_genes-transcript_Equivalence.txt
cat AsianElephant_MitoMammoth_tRNA_genes-transcript_Equivalence.txt | awk '{print $2}' | sort | uniq > AsianElephant_MitoMammoth_tRNA_genes-transcript_IDs.txt



#Merge all gene types (protein coding, rRNA, snRNA, snoRNA, lncRNA, tRNA and miRNA)
cat AsianElephant_MitoMammoth_ProteinCoding_genes.bed AsianElephant_MitoMammoth_rRNA_genes.bed AsianElephant_MitoMammoth_snRNA_genes.bed AsianElephant_MitoMammoth_snoRNA_genes.bed AsianElephant_MitoMammoth_lncRNA_genes.bed AsianElephant_MitoMammoth_tRNA_genes.bed  AsianElephant_MitoMammoth_miRNA_genes.bed >  AsianElephant_MitoMammoth_PC-NC-miRNA_genes.bed

echo "DONE..."



#Generate BED files by gene type (mitochondrial genes)
echo ""
echo "> > > > > > Generating BED files for mitochondrial genes"
cat DQ188829.2_MitoMammoth.bed | rg protein_coding > DQ188829.2_MitoMammoth_ProteinCoding.bed
awk '{print $1,1,$3-$2+1, $4,$5,$6}' DQ188829.2_MitoMammoth_ProteinCoding.bed | gsed 's/ /\t/g' > DQ188829.2_MitoMammoth_ProteinCoding_Transcriptome.bed
cat DQ188829.2_MitoMammoth.bed | rg rRNA > DQ188829.2_MitoMammoth_rRNA.bed
awk '{print $1,1,$3-$2+1, $4,$5,$6}' DQ188829.2_MitoMammoth_rRNA.bed | gsed 's/ /\t/g' > DQ188829.2_MitoMammoth_rRNA_Transcriptome.bed 
cat DQ188829.2_MitoMammoth.bed | rg tRNA > DQ188829.2_MitoMammoth_tRNA.bed
awk '{print $1,1,$3-$2+1, $4,$5,$6}' DQ188829.2_MitoMammoth_tRNA.bed | gsed 's/ /\t/g' > DQ188829.2_MitoMammoth_tRNA_Transcriptome.bed

echo "DONE..."



#Generate Linear Fastas and BED by gene type with proper header ID and longest per gene (if applicable)
echo ""
echo "> > > > > > Generating Fasta and BED files (longest per gene)"
#Mitochondrial genes
bedtools getfasta -fi AsianElephant_MitoMammoth.fa -bed DQ188829.2_MitoMammoth_rRNA.bed -s -name | gsed 's/:/\t/g' | awk '{print $1}' > DQ188829.2_MitoMammoth_rRNA.fa
bedtools getfasta -fi AsianElephant_MitoMammoth.fa -bed DQ188829.2_MitoMammoth_tRNA.bed -s -name | gsed 's/:/\t/g' | awk '{print $1}' > DQ188829.2_MitoMammoth_tRNA.fa
bedtools getfasta -fi AsianElephant_MitoMammoth.fa -bed DQ188829.2_MitoMammoth_ProteinCoding.bed -s -name | gsed 's/:/\t/g' | awk '{print $1}' > DQ188829.2_MitoMammoth_ProteinCoding.fa

#Nuclear genes
#rRNA
rg -w -f AsianElephant_MitoMammoth_rRNA_genes-transcript_IDs.txt AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths.txt > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_rRNA.txt
rg -w -A 1 -f AsianElephant_MitoMammoth_rRNA_genes-transcript_IDs.txt AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear.fa | grep -v "\--" | seqkit sort --line-width 0 > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_rRNA.fa
rg ">" AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_rRNA.fa > header
rg -v ">" AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_rRNA.fa > seqs
paste AsianElephant_MitoMammoth_rRNA_genes-transcript_Equivalence.txt header | awk '{print ">"$1"|"$2}' > headerf
paste headerf seqs | awk '{print $1 ORS $2}' > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_rRNA_FINAL.fa
paste headerf AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_rRNA.txt | awk '{print $1,1,$3}' | gsed 's/ /\t/g' | gsed 's/>//g'  > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_rRNA.bed

#snRNA
rg -w -f AsianElephant_MitoMammoth_snRNA_genes-transcript_IDs.txt AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths.txt > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_snRNA.txt
rg -w -A 1 -f AsianElephant_MitoMammoth_snRNA_genes-transcript_IDs.txt AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear.fa | grep -v "\--" | seqkit sort --line-width 0 > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_snRNA.fa
rg ">" AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_snRNA.fa > header
rg -v ">" AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_snRNA.fa > seqs
paste AsianElephant_MitoMammoth_snRNA_genes-transcript_Equivalence.txt header | awk '{print ">"$1"|"$2}' > headerf
paste headerf seqs | awk '{print $1 ORS $2}' > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_snRNA_FINAL.fa
paste headerf AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_snRNA.txt | awk '{print $1,1,$3}' | gsed 's/ /\t/g' | gsed 's/>//g' > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_snRNA.bed

#snoRNA
rg -w -f AsianElephant_MitoMammoth_snoRNA_genes-transcript_IDs.txt AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths.txt > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_snoRNA.txt
rg -w -A 1 -f AsianElephant_MitoMammoth_snoRNA_genes-transcript_IDs.txt AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear.fa | grep -v "\--" | seqkit sort --line-width 0 > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_snoRNA.fa
rg ">" AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_snoRNA.fa > header
rg -v ">" AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_snoRNA.fa > seqs
paste AsianElephant_MitoMammoth_snoRNA_genes-transcript_Equivalence.txt header | awk '{print ">"$1"|"$2}' > headerf
paste headerf seqs | awk '{print $1 ORS $2}' > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_snoRNA_FINAL.fa
paste headerf AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_snoRNA.txt | awk '{print $1,1,$3}' | gsed 's/ /\t/g' | gsed 's/>//g' > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_snoRNA.bed

#tRNA
rg -w -f AsianElephant_MitoMammoth_tRNA_genes-transcript_IDs.txt AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths.txt > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_tRNA.txt
rg -w -A 1 -f AsianElephant_MitoMammoth_tRNA_genes-transcript_IDs.txt AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear.fa | grep -v "\--" | seqkit sort --line-width 0 > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_tRNA.fa
rg ">" AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_tRNA.fa > header
rg -v ">" AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_tRNA.fa > seqs
paste AsianElephant_MitoMammoth_tRNA_genes-transcript_Equivalence.txt header | awk '{print ">"$1"|"$2}' > headerf
paste headerf seqs | awk '{print $1 ORS $2}' > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_tRNA_FINAL.fa
paste headerf AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_tRNA.txt | awk '{print $1,1,$3}' | gsed 's/ /\t/g' | gsed 's/>//g' > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_tRNA.bed

#miRNA
awk '{print $1,$2-30,$3+30,$4,$5,$6}' AsianElephant_MitoMammoth_miRNA_genes.bed | gsed 's/_pre/_pri/g' | awk '{if($3<0) {print $1,1,$3,$4,$5,$6} else {print $1,$2,$3,$4,$5,$6}}' | gsed 's/ /\t/g' > AsianElephant_MitoMammoth_miRNA_genes_pri_+-30nt.bed
bedtools getfasta -fi AsianElephant_MitoMammoth.fa -bed AsianElephant_MitoMammoth_miRNA_genes_pri_+-30nt.bed -s -name | gsed 's/:/\t/g' | awk '{print $1}' > AsianElephant_MitoMammoth_miRNA_genes_pri_+-30nt.fa
cat AsianElephant_MitoMammoth_miRNA_genes_pri_+-30nt.bed | awk '{print $4,1,$3-$2+1}' | gsed 's/ /\t/g' | sort -k1,1 > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_miRNA.bed
awk '{print $4}' AsianElephant_MitoMammoth_miRNA_genes_pri_+-30nt.bed > AsianElephant_MitoMammoth_miRNA_genes_pri_+-30nt_IDs.txt

#lncRNA
rg -w -f AsianElephant_MitoMammoth_lncRNA_genes-transcript_IDs.txt AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths.txt > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA.txt
paste AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA.txt AsianElephant_MitoMammoth_lncRNA_genes-transcript_Equivalence.txt | awk '{print $3,$4,$2}' | gsed 's/ /\t/g' > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL.txt
cat AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL.txt | awk '{if ($3 > max[$1]) {max[$1] = $3; id[$1] = $0}} END {for (i in id) print id[i]}' | sort -k2,2 > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene.txt
awk '{print $2}' AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene.txt > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene_IDs.txt

seqkit grep -f AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene_IDs.txt AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear.fa | seqkit sort --line-width 0 > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene.fa
rg ">" AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene.fa > header
rg -v ">" AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene.fa > seqs
paste AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene.txt header | awk '{print ">"$1"|"$2}' > headerf
paste headerf seqs | awk '{print $1 ORS $2}' > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene_FINAL.fa

cat AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene.txt | awk '{print $1"|"$2,1,$3}' | gsed 's/ /\t/g' > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene.bed


#ProteinCoding
grep -w -f AsianElephant_MitoMammoth_ProteinCoding_genes-transcript_IDs.txt AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths.txt > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding.txt
paste AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding.txt AsianElephant_MitoMammoth_ProteinCoding_genes-transcript_Equivalence.txt | awk '{print $3,$4,$2}' | gsed 's/ /\t/g' > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL.txt
cat AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL.txt | awk '{if ($3 > max[$1]) {max[$1] = $3; id[$1] = $0}} END {for (i in id) print id[i]}' | sort -k2,2 > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene.txt
awk '{print $2}' AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene.txt > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene_IDs.txt

seqkit grep -f AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene_IDs.txt AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear.fa | seqkit sort --line-width 0 > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene.fa
rg ">" AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene.fa > header
rg -v ">" AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene.fa > seqs
paste AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene.txt header | awk '{print ">"$1"|"$2}' > headerf
paste headerf seqs | awk '{print $1 ORS $2}' > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene_FINAL.fa

cat AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene.txt | awk '{print $1"|"$2,1,$3}' | gsed 's/ /\t/g' > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene.bed



rm header
rm headerf
rm seqs



#Merge all gene types longest per gene (protein coding, rRNA, snRNA, snoRNA, lncRNA, tRNA and miRNA)
cat DQ188829.2_MitoMammoth_rRNA.fa DQ188829.2_MitoMammoth_tRNA.fa DQ188829.2_MitoMammoth_ProteinCoding.fa AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_rRNA_FINAL.fa AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_snRNA_FINAL.fa AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_snoRNA_FINAL.fa AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_tRNA_FINAL.fa AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene_FINAL.fa AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene_FINAL.fa AsianElephant_MitoMammoth_miRNA_genes_pri_+-30nt.fa > AsianElephant_MitoMammoth_PC-NC-miRNA_genes_Transcriptome.fa

cat AsianElephant_MitoMammoth_rRNA_genes-transcript_IDs.txt AsianElephant_MitoMammoth_tRNA_genes-transcript_IDs.txt AsianElephant_MitoMammoth_snRNA_genes-transcript_IDs.txt AsianElephant_MitoMammoth_snoRNA_genes-transcript_IDs.txt AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene_IDs.txt AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene_IDs.txt AsianElephant_MitoMammoth_miRNA_genes_pri_+-30nt_IDs.txt > AsianElephant_MitoMammoth_PC-NC-miRNA_genes_transcript_IDs.txt

echo "DONE..."



#Build Bowtie2 index for Transcriptome
echo ""
echo "> > > > > > Builidng Bowtie2 index for Transcriptome"
bowtie2-build --threads 8 AsianElephant_MitoMammoth_PC-NC-miRNA_genes_Transcriptome.fa AsianElephant_MitoMammoth_PC-NC-miRNA_genes_Transcriptome

echo "DONE..."

#Get Intronic longest per gene regions
echo ""
echo "> > > > > > Calculating intronic regions"
cat AsianElephant_MitoMammoth_intronCoords.txt | awk '{print $1,$2,$3,$5,0,$4}' | gsed 's/ /\t/g' | sort -f -k1,1 -k2,2 -V > AsianElephant_MitoMammoth_intronCoords.bed
rg -w -f AsianElephant_MitoMammoth_PC-NC-miRNA_genes_transcript_IDs.txt AsianElephant_MitoMammoth_intronCoords.bed > AsianElephant_MitoMammoth_intronCoords_Transcriptome_Longest-per_Gene_Genomic.bed

echo "DONE..."

#Get Exonic longest per gene regions
echo ""
echo "> > > > > > Calculating exonic regions"
rg -F -f AsianElephant_MitoMammoth_PC-NC-miRNA_genes_transcript_IDs.txt AsianElephant_MitoMammoth.gff3 | awk '$3=="exon" {exon_num = "NA"; if (match($0, /exon_number "[^"]+"/)) {exon_num = substr($0, RSTART+13, RLENGTH-13);} print $1, $4-1, $5, $10 "|" $12 "_exon_" exon_num, 0, $7}' | gsed 's/"//g' | gsed 's/;//g' | gsed 's/ /\t/g' > temp
cat temp AsianElephant_MitoMammoth_miRNA_genes_pri_+-30nt.bed DQ188829.2_MitoMammoth_rRNA.bed DQ188829.2_MitoMammoth_tRNA.bed DQ188829.2_MitoMammoth_ProteinCoding.bed > AsianElephant_MitoMammoth_ExonCoords_Transcriptome_Longest-per_Gene_Genomic.bed

rm temp

echo "DONE..."



#Get Exon-Exon junction transcriptome BED
echo ""
echo "> > > > > > Generating exon-exon junction BED files"
rg -w -f AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene_IDs.txt AsianElephant_MitoMammoth_Transcriptome.txt | rg , > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene_IDs_ExonStructure.txt
rg -w -f AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene_IDs.txt AsianElephant_MitoMammoth_Transcriptome.txt | rg , > AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene_IDs_ExonStructure.txt

# Input and output files
input_file="AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene_IDs_ExonStructure.txt"           # Input file with the strings
output_file="AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene_Exon-Exon_junctions.bed"      # Output file

# Empty the output file
> "$output_file"

# Process each line
while IFS= read -r line; do
    # First field
    first_field=$(echo "$line" | cut -d' ' -f1)
    first_field="${first_field#>}"  # remove '>' only at start

    # Extract the segs field
    segs_string=$(echo "$line" | sed -n 's/.*segs:\([^ ]*\).*/\1/p')

    # Skip if no segs found
    if [[ -z "$segs_string" ]]; then
        continue
    fi

    # Split segs into an array
    IFS=',' read -ra segs <<< "$segs_string"

    # Iterate and print junctions
    for ((i=0; i<${#segs[@]}-1; i++)); do
        current=${segs[i]}
        next=${segs[i+1]}
        current_end=$(echo "$current" | cut -d'-' -f2)
        next_start=$(echo "$next" | cut -d'-' -f1)
        echo -e "${first_field}\t${current_end}\t${next_start}" >> "$output_file"
    done
done < "$input_file"


# Input and output files
input_file="AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene_IDs_ExonStructure.txt"           # Input file with the strings
output_file="AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene_Exon-Exon_junctions.bed"      # Output file

# Empty the output file
> "$output_file"

# Process each line
while IFS= read -r line; do
    # First field
    first_field=$(echo "$line" | cut -d' ' -f1)
    first_field="${first_field#>}"  # remove '>' only at start

    # Extract the segs field
    segs_string=$(echo "$line" | sed -n 's/.*segs:\([^ ]*\).*/\1/p')

    # Skip if no segs found
    if [[ -z "$segs_string" ]]; then
        continue
    fi

    # Split segs into an array
    IFS=',' read -ra segs <<< "$segs_string"

    # Iterate and print junctions
    for ((i=0; i<${#segs[@]}-1; i++)); do
        current=${segs[i]}
        next=${segs[i+1]}
        current_end=$(echo "$current" | cut -d'-' -f2)
        next_start=$(echo "$next" | cut -d'-' -f1)
        echo -e "${first_field}\t${current_end}\t${next_start}" >> "$output_file"
    done
done < "$input_file"

#Merge Exon-Exon junctions and fix the ID
cat AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene_Exon-Exon_junctions.bed AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene_Exon-Exon_junctions.bed | sort -f -k1,1 -k2,2 -V > AsianElephant_MitoMammoth_Transcriptome_PC-NC_Longest-per_Gene_Exon-Exon_junctions.bed

sort -k2,2 AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene.txt > temp
sort -k2,2 AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene.txt > temp2
cat temp temp2 | sort -f -k2,2 > temp3
join -1 1 -2 2 -t $'\t' AsianElephant_MitoMammoth_Transcriptome_PC-NC_Longest-per_Gene_Exon-Exon_junctions.bed temp3 | awk '{print $4"|"$1,$2,$3}' | gsed 's/ /\t/g' > AsianElephant_MitoMammoth_Transcriptome_PC-NC_Longest-per_Gene_Exon-Exon_junctions_FINAL.bed

rm temp*

echo "DONE..."



#Get Exon-Intron junctions genome-wide
echo ""
echo "> > > > > > Generating exon-intron junction BED files"
# Input and output files
input_file="AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene_IDs_ExonStructure.txt"           # Input file with the strings
output_file="AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene_Exon-Intron_junctions.bed"      # Output file

# Empty the output file first
> "$output_file"

while IFS= read -r line; do
    # Extract gene ID and remove leading '>'
    gene_id=$(echo "$line" | cut -d' ' -f1)
    gene_id="${gene_id#>}"

    # Extract scaffold name after 'loc:' (can be any string like HiC_scaffold_38 or chr1)
    scaffold=$(echo "$line" | sed -n 's/.*loc:\([^|]*\).*/\1/p')

    # Extract exon intervals
    exon_field=$(echo "$line" | grep -o 'exons:[^ ]*')
    exon_string=${exon_field#exons:}

    # Skip if no exons field
    if [[ -z "$exon_string" ]]; then
        continue
    fi

    # Split exons into array
    IFS=',' read -ra exons <<< "$exon_string"
    exon_count=${#exons[@]}

    # Loop through all but the last exon
    for ((i=0; i<exon_count-1; i++)); do
        exon_end=$(echo "${exons[i]}" | cut -d'-' -f2)
        exon_next_start=$((exon_end + 1))
        # Print scaffold, exon_end, exon_next_start, and gene_id to the output file
        echo -e "${scaffold}\t${exon_end}\t${exon_next_start}\t${gene_id}" >> "$output_file"
    done
done < "$input_file"


# Input and output files
input_file="AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene_IDs_ExonStructure.txt"           # Input file with the strings
output_file="AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene_Exon-Intron_junctions.bed"      # Output file

# Empty the output file first
> "$output_file"

while IFS= read -r line; do
    # Extract gene ID and remove leading '>'
    gene_id=$(echo "$line" | cut -d' ' -f1)
    gene_id="${gene_id#>}"

    # Extract scaffold name after 'loc:' (can be any string like HiC_scaffold_38 or chr1)
    scaffold=$(echo "$line" | sed -n 's/.*loc:\([^|]*\).*/\1/p')

    # Extract exon intervals
    exon_field=$(echo "$line" | grep -o 'exons:[^ ]*')
    exon_string=${exon_field#exons:}

    # Skip if no exons field
    if [[ -z "$exon_string" ]]; then
        continue
    fi

    # Split exons into array
    IFS=',' read -ra exons <<< "$exon_string"
    exon_count=${#exons[@]}

    # Loop through all but the last exon
    for ((i=0; i<exon_count-1; i++)); do
        exon_end=$(echo "${exons[i]}" | cut -d'-' -f2)
        exon_next_start=$((exon_end + 1))
        # Print scaffold, exon_end, exon_next_start, and gene_id to the output file
        echo -e "${scaffold}\t${exon_end}\t${exon_next_start}\t${gene_id}" >> "$output_file"
    done
done < "$input_file"

cat AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene_Exon-Intron_junctions.bed AsianElephant_MitoMammoth_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene_Exon-Intron_junctions.bed | sort -f -k1,1 -k2,2 -V > AsianElephant_MitoMammoth_Transcriptome_PC-NC_Longest-per_Gene_Genomic_Exon-Intron_junctions.bed

echo "DONE..."

