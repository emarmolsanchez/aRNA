#! /bin/bash/

#Alignment to genome and sort
echo ""
echo "> > > > > > Mapping to the Asian elephant genome"
bowtie2 --end-to-end --sensitive --no-unal -p 8 -x Reference-assisted_3D_woolly_mammoth_assembly -f Mammoth1_Trimmed_Collapsed.fasta -S Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly.sam 2>Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_alignment.log

echo "DONE..."

echo ""
echo "> > > > > > Convert to BAM and sort"
sambamba view -f bam -o Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly.bam -S Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly.sam
sambamba sort -o Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_sorted.bam Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly.bam

echo "DONE..."


#UMI deduplication
echo ""
echo "> > > > > > UMI deduplication"
umi_tools dedup --output-stats=Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Stats.txt --stdin=Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_sorted.bam --log=Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP.log > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP.bam

sambamba view -f sam -o Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP.sam Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP.bam

echo "DONE..."

#Generate BED from alignment file
echo ""
echo "> > > > > > Generating BED from alingment files"
cat Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP.sam | grep -v "@" | awk '{FS=OFS="\t"}{print $3, $4, length($10)+$4-1, $1,$10}' > temp
cat Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP.sam | grep -v "@" | awk '{FS=OFS="\t"}{print $2}' | gsed 's/0/+/g' | gsed 's/16/-/g' > temp2
paste temp temp2 > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP.bed

rm *.sam
rm temp*

echo "DONE..."




#Determine Genic, Intergenic & Mitochondrial sequences
echo ""
echo "> > > > > > Determining genic, intergenic and mitochondrial sequences"
#Genic
bedtools intersect -a Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP.bed -b Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes.bed -s -wao | rg -v "\t\-1" | rg -v "DQ188829.2" | awk '{print $4}' | sort | uniq > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic_ReadIDs.txt
rg -f Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic_ReadIDs.txt Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP.bed > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic.bed

#Intergenic
bedtools intersect -a Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP.bed -b Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes.bed -s -wao | rg "\t\-1" | rg -v "DQ188829.2" | awk '{print $4}' | sort | uniq > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Intergenic_ReadIDs.txt
rg -f Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Intergenic_ReadIDs.txt Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP.bed > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Intergenic.bed

#Mitochondrial
cat Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP.bed | rg "DQ188829" | awk '{print $4}' | sort | uniq > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Mitochondrial_ReadIDs.txt
rg -f Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Mitochondrial_ReadIDs.txt Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP.bed > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Mitochondrial.bed

echo "DONE..."



#Quantification & Coverage of Mitochondrial genes
echo ""
echo "> > > > > > Quantification & Coverage of Mitochondrial genes"
#Protein coding
bedtools coverage -a DQ188829.2_MitoMammoth_ProteinCoding.bed -b Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Mitochondrial.bed -wao | awk '{FS=OFS="\t"}{print $4,"protein_coding",$7,$10}' > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Mitochondrial_ProteinCoding_Count-Coverage.txt

#tRNA
bedtools coverage -a DQ188829.2_MitoMammoth_tRNA.bed -b Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Mitochondrial.bed -wao | awk '{FS=OFS="\t"}{print $4,"tRNA",$7,$10}' > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Mitochondrial_tRNA_Count-Coverage.txt

#rRNA
bedtools coverage -a DQ188829.2_MitoMammoth_rRNA.bed -b Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Mitochondrial.bed -wao | awk '{FS=OFS="\t"}{print $4,"rRNA",$7,$10}' > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Mitochondrial_rRNA_Count-Coverage.txt

echo "DONE..."



#Determine genomic Exonic & Intronic
echo ""
echo "> > > > > > Determining exonic and intronic sequences genome-wide"
#Exonic
bedtools intersect -a Reference-assisted_3D_woolly_mammoth_assembly_ExonCoords_Transcriptome_Longest-per_Gene_Genomic.bed -b Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic.bed -s -wao | rg -v "\t\-1" | rg -v "DQ188829" | awk '{print $10}' | sort | uniq > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic_Exonic_ReadIDs.txt
rg -f Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic_Exonic_ReadIDs.txt Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic.bed > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic_Exonic.bed

#Intronic
bedtools intersect -a Reference-assisted_3D_woolly_mammoth_assembly_intronCoords_Transcriptome_Longest-per_Gene_Genomic.bed -b Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic.bed -s -wao | grep -v "\t\-1" | grep -v "DQ188829.2" | awk '{print $10}' | sort | uniq > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic_Intronic_ReadIDs.txt
rg -f Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic_Intronic_ReadIDs.txt Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic.bed > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic_Intronic.bed

echo "DONE..."



#Alignment to Transcriptome and sort
echo ""
echo "> > > > > > Mapping to the Asian elephant transcriptome"
bowtie2 --end-to-end --sensitive --no-unal --norc -p 8 -x Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome -f Mammoth1_Trimmed_Collapsed.fasta -S Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome.sam 2>Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_alignment.log

echo "DONE..."

echo ""
echo "> > > > > > Convert to BAM and sort"
sambamba view -f bam -o Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome.bam -S Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome.sam
sambamba sort -o Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_sorted.bam Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome.bam

echo "DONE..."


#UMI deduplication
echo ""
echo "> > > > > > UMI deduplication"
umi_tools dedup --output-stats=Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_Stats.txt --stdin=Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_sorted.bam --log=Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP.log > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP.bam

sambamba view -f sam -o Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP.sam Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP.bam

echo "DONE..."



#Generate BED from alignment file
echo ""
echo "> > > > > > Generating BED from alingment files"
cat Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP.sam | grep -v "@" | awk '{FS=OFS="\t"}{print $3, $4, length($10)+$4-1, $1,$10}' > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP.bed

rm *.sam

echo "DONE..."



#Clean Intergenic, Intronic & Mitochondrial sequences from transcriptome alignment
echo ""
echo "> > > > > > Cleaning intergenic, intronic & mitochondrial sequences from transcriptome alignment"
awk 'BEGIN { if (getline == 0) print "none"; else { print; while (getline) print } close(FILENAME) }' Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Mitochondrial_ReadIDs.txt > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Mitochondrial_ReadIDs2.txt
rg -v -f Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Mitochondrial_ReadIDs2.txt Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP.bed > temp

awk 'BEGIN { if (getline == 0) print "none"; else { print; while (getline) print } close(FILENAME) }'  Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Intergenic_ReadIDs.txt >  Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Intergenic_ReadIDs2.txt
rg -v -f Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Intergenic_ReadIDs2.txt temp > temp2

awk 'BEGIN { if (getline == 0) print "none"; else { print; while (getline) print } close(FILENAME) }'  Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic_Intronic_ReadIDs.txt > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic_Intronic_ReadIDs2.txt
rg -v -f Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic_Intronic_ReadIDs2.txt temp2 > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN.bed

rm temp*
rm *_ReadIDs2.txt

echo "DONE..."



#Generate BED from clean data
echo ""
echo "> > > > > > Generating BED from clean files"
#Protein coding
rg -f Reference-assisted_3D_woolly_mammoth_assembly_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene_IDs.txt Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN.bed > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_ProteinCoding.bed

#rRNA
rg -f Reference-assisted_3D_woolly_mammoth_assembly_rRNA_genes-transcript_IDs.txt Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN.bed > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_rRNA.bed

#tRNA
rg -f Reference-assisted_3D_woolly_mammoth_assembly_tRNA_genes-transcript_IDs.txt Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN.bed > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_tRNA.bed

#snRNA
rg -f Reference-assisted_3D_woolly_mammoth_assembly_snRNA_genes-transcript_IDs.txt Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN.bed > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_snRNA.bed

#snoRNA
rg -f Reference-assisted_3D_woolly_mammoth_assembly_snoRNA_genes-transcript_IDs.txt Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN.bed > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_snoRNA.bed

#lncRNA
rg -f Reference-assisted_3D_woolly_mammoth_assembly_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene_IDs.txt Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN.bed > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_lncRNA.bed

#miRNA
rg "Ema" Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN.bed > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_miRNA.bed
bedtools intersect -a Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_miRNA.bed -b Reference-assisted_3D_woolly_mammoth_assembly_Transcriptome_SingleHeader_Linear_Lengths_miRNA.bed -wao | awk '{if($3>$8-34 && $3>$8-34-23) {print}}' > temp_3p
bedtools intersect -a Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_miRNA.bed -b Reference-assisted_3D_woolly_mammoth_assembly_Transcriptome_SingleHeader_Linear_Lengths_miRNA.bed -wao | awk '{if($2>31-3 && $2<54+3) {print}}' > temp_5p
cat temp_5p temp_3p > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_miRNA_Corrected.bed

rm temp*

echo "DONE..."



#Quantification & Coverage of nuclear genes
echo ""
echo "> > > > > > Quantification & Coverage of Nuclear genes"
#miRNA
cat Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_miRNA_Corrected.bed | awk '{print $1}' | sort | uniq -c | gsed 's/^ *//' | awk '{print $2,$1}' | gsed 's/ /\t/g' > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_miRNA_Corrected_Count.txt

#Proteiin coding
bedtools coverage -a Reference-assisted_3D_woolly_mammoth_assembly_Transcriptome_SingleHeader_Linear_Lengths_ProteinCoding_FINAL_Longest-per_Gene.bed -b Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_ProteinCoding.bed -wao | awk '{FS=OFS="\t"}{print $1,"protein_coding",$4,$7}' | sort -k1,1 > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_ProteinCoding_Count-Coverage.txt

#lncRNA
bedtools coverage -a Reference-assisted_3D_woolly_mammoth_assembly_Transcriptome_SingleHeader_Linear_Lengths_lncRNA_FINAL_Longest-per_Gene.bed -b Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_lncRNA.bed -wao | awk '{FS=OFS="\t"}{print $1,"lncRNA",$4,$7}' | sort -k1,1 > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_lncRNA_Count-Coverage.txt

#rRNA
bedtools coverage -a Reference-assisted_3D_woolly_mammoth_assembly_Transcriptome_SingleHeader_Linear_Lengths_rRNA.bed -b Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_rRNA.bed -wao | awk '{FS=OFS="\t"}{print $1,"rRNA",$4,$7}' | sort -k1,1 > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_rRNA_Count-Coverage.txt

#tRNA
bedtools coverage -a Reference-assisted_3D_woolly_mammoth_assembly_Transcriptome_SingleHeader_Linear_Lengths_tRNA.bed -b Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_tRNA.bed -wao | awk '{FS=OFS="\t"}{print $1,"tRNA",$4,$7}' | sort -k1,1 > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_tRNA_Count-Coverage.txt

#snRNA
bedtools coverage -a Reference-assisted_3D_woolly_mammoth_assembly_Transcriptome_SingleHeader_Linear_Lengths_snRNA.bed -b Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_snRNA.bed -wao | awk '{FS=OFS="\t"}{print $1,"snRNA",$4,$7}' | sort -k1,1 > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_snRNA_Count-Coverage.txt

#snoRNA
bedtools coverage -a Reference-assisted_3D_woolly_mammoth_assembly_Transcriptome_SingleHeader_Linear_Lengths_snoRNA.bed -b Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_snoRNA.bed -wao | awk '{FS=OFS="\t"}{print $1,"snoRNA",$4,$7}' | sort -k1,1 > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_snoRNA_Count-Coverage.txt

echo "DONE..."



#Determine Exon-Exon & Exon-Intron overlapping sequences
echo ""
echo "> > > > > > Determining Exon-Exon & Exon-Intron overlaps"
#Exon-Exon
bedtools intersect -a Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN.bed -b Reference-assisted_3D_woolly_mammoth_assembly_Transcriptome_PC-NC_Longest-per_Gene_Exon-Exon_junctions_FINAL.bed -wao | rg -v "\-1" | awk '{if($7-$2>=4 && $3-$8>=4) {print $1,$2,$3,$4,$5}}' | gsed 's/ /\t/g' > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_PC-NC-miRNA_genes_Transcriptome_UMIDEDUP_CLEAN_Exon-Exon_junction_Overlaps.bed

#Exon-Intron
bedtools intersect -a Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic.bed -b Reference-assisted_3D_woolly_mammoth_assembly_Transcriptome_PC-NC_Longest-per_Gene_Genomic_Exon-Intron_junctions.bed -wao | rg -v "\-1" | awk '{if($8-$2>=4 && $3-$9>=4) {print $1,$2,$3,$4,$5}}' | gsed 's/ /\t/g' > Mammoth1_Trimmed_Collapsed_vs_Reference-assisted_3D_woolly_mammoth_assembly_UMIDEDUP_Nuclear_Genic_Exon-Intron_Overlaps.bed

echo "DONE..."



