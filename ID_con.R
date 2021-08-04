#Data1 #GRCh37.p13 non coding 
 #awk conversion

#awk '{                                
 #   for (i = 1; i <= NF; i++) {
  #      if ($i ~ /gene_id|gene_name/) {
   #         printf "%s ", $(i+1)
    #    }
    #}
    #print ""
#}' Homo_sapiens.GRCh37.70.gtf | sed -e 's/"//g' -e 's/;//g' -e 's/ /\t/' | sort -k1,1 | uniq > Homo_sapiens.GRCh37.70.txt
#for RNA seq samples
#Data1 #GRCh37.p13 non coding 
setwd("~/GenomeDK/faststorage/reference/UCSC/human")
annot = read.delim('gencode.v19.long_noncoding_RNAs_sorted_filtered_modified.gtf', header=F)
colnames(annot) <- c("Gene_ID","gene_type","gene_name","transcript_type","transcript_status")

setwd("~/GenomeDK/faststorage/F20_new/Sym/F20FTSEUHT1552_HUMjmvN/Clean")
S1=as.matrix(read.table("1_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1))
S2=as.matrix(read.table("2_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1))
S3=as.matrix(read.table("3_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1))
S4=as.matrix(read.table("4_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1))
S5=as.matrix(read.table("5_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1))
S6=as.matrix(read.table("6_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1))
S7=as.matrix(read.table("7_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1))
S8=as.matrix(read.table("8_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1))
S9=as.matrix(read.table("9_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1))
S10=as.matrix(read.table("10_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1))
S11=as.matrix(read.table("11_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1))
S12=as.matrix(read.table("12_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1))

col=cbind(S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12)
colnames(col)=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12")
col <- as.data.frame(cbind(rownames(col), col)) # rowname is a metadata and rest all the rows are data so to convert metadata to data one should run this command 
rownames(col) <- NULL
colnames(col) <- c("Gene_id","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12")

merged.file = merge(col, annot, by.x='Gene_id', by.y='Gene_ID')
write.table(merged.file, file = "All_annot.txt", row.names = F, sep = "\t", quote = F)

#in awk
#cat All_annot.txt | awk '{counter1[$16] += $2} {counter2[$16] += $3} {counter3[$16] += $4} {counter4[$16] += $5} {counter5[$16] += $6} {counter6[$16] += $7} {counter7[$16] += $8} {counter8[$16] += $9} {counter9[$16] += $10} {counter10[$16] += $11} {counter11[$16] += $12} {counter12[$16] += $13} END {for(i in counter1){print i"\t"counter1[i]"\t"counter2[i]"\t"counter3[i]"\t"counter4[i]"\t"counter5[i]"\t"counter6[i]"\t"counter7[i]"\t"counter8[i]"\t"counter9[i]"\t"counter10[i]"\t"counter11[i]"\t"counter12[i]}}' >> Count_annotated_new_RNA_species_count.txt 

#Data2
#for htacseq samples #GRCh37.p13 non coding 
setwd("~/GenomeDK/faststorage/reference/UCSC/human")
annot = read.delim('gencode.v19.long_noncoding_RNAs_sorted_filtered_modified.gtf', header=F)
colnames(annot) <- c("Gene_ID","gene_type","gene_name","transcript_type","transcript_status")

setwd("~/GenomeDK/faststorage/chromatin_retension/All/07_bamfilter2")
S1=as.matrix(read.table("200811_I29_V300060997_L3_a0_1_1.fq.gz.bbduk.fq.gz_Aligned.sortedByCoord.out.bam_new.bam_filter.bam.sam_htseq_count_union_nc.txt", header = FALSE, row.names = 1))
S2=as.matrix(read.table("200811_I29_V300060997_L3_a0_2_1.fq.gz.bbduk.fq.gz_Aligned.sortedByCoord.out.bam_new.bam_filter.bam.sam_htseq_count_union_nc.txt", header = FALSE, row.names = 1))
S3=as.matrix(read.table("200811_I29_V300060997_L3_a0_3_1.fq.gz.bbduk.fq.gz_Aligned.sortedByCoord.out.bam_new.bam_filter.bam.sam_htseq_count_union_nc.txt", header = FALSE, row.names = 1))
S4=as.matrix(read.table("200811_I29_V300060997_L3_B0_1_1.fq.gz.bbduk.fq.gz_Aligned.sortedByCoord.out.bam_new.bam_filter.bam.sam_htseq_count_union_nc.txt", header = FALSE, row.names = 1))
S5=as.matrix(read.table("200811_I29_V300060997_L3_B0_2_1.fq.gz.bbduk.fq.gz_Aligned.sortedByCoord.out.bam_new.bam_filter.bam.sam_htseq_count_union_nc.txt", header = FALSE, row.names = 1))
S6=as.matrix(read.table("200811_I29_V300060997_L3_B0_3_1.fq.gz.bbduk.fq.gz_Aligned.sortedByCoord.out.bam_new.bam_filter.bam.sam_htseq_count_union_nc.txt", header = FALSE, row.names = 1))
S7=as.matrix(read.table("200811_I29_V300060997_L3_m0_1_1.fq.gz.bbduk.fq.gz_Aligned.sortedByCoord.out.bam_new.bam_filter.bam.sam_htseq_count_union_nc.txt", header = FALSE, row.names = 1))
S8=as.matrix(read.table("200811_I29_V300060997_L3_m0_3_1.fq.gz.bbduk.fq.gz_Aligned.sortedByCoord.out.bam_new.bam_filter.bam.sam_htseq_count_union_nc.txt", header = FALSE, row.names = 1))
S9=as.matrix(read.table("200811_I29_V300060997_L4_a120_1_1.fq.gz.bbduk.fq.gz_Aligned.sortedByCoord.out.bam_new.bam_filter.bam.sam_htseq_count_union_nc.txt", header = FALSE, row.names = 1))
S10=as.matrix(read.table("200811_I29_V300060997_L4_a120_2_1.fq.gz.bbduk.fq.gz_Aligned.sortedByCoord.out.bam_new.bam_filter.bam.sam_htseq_count_union_nc.txt", header = FALSE, row.names = 1))
S11=as.matrix(read.table("200811_I29_V300060997_L4_a120_3_1.fq.gz.bbduk.fq.gz_Aligned.sortedByCoord.out.bam_new.bam_filter.bam.sam_htseq_count_union_nc.txt", header = FALSE, row.names = 1))
S12=as.matrix(read.table("200811_I29_V300060997_L4_B120_1_1.fq.gz.bbduk.fq.gz_Aligned.sortedByCoord.out.bam_new.bam_filter.bam.sam_htseq_count_union_nc.txt", header = FALSE, row.names = 1))
S13=as.matrix(read.table("200811_I29_V300060997_L4_B120_2_1.fq.gz.bbduk.fq.gz_Aligned.sortedByCoord.out.bam_new.bam_filter.bam.sam_htseq_count_union_nc.txt", header = FALSE, row.names = 1))
S14=as.matrix(read.table("200811_I29_V300060997_L4_B120_3_1.fq.gz.bbduk.fq.gz_Aligned.sortedByCoord.out.bam_new.bam_filter.bam.sam_htseq_count_union_nc.txt", header = FALSE, row.names = 1))
S15=as.matrix(read.table("200811_I29_V300060997_L4_m120_1_1.fq.gz.bbduk.fq.gz_Aligned.sortedByCoord.out.bam_new.bam_filter.bam.sam_htseq_count_union_nc.txt", header = FALSE, row.names = 1))
S16=as.matrix(read.table("200811_I29_V300060997_L4_m120_2_1.fq.gz.bbduk.fq.gz_Aligned.sortedByCoord.out.bam_new.bam_filter.bam.sam_htseq_count_union_nc.txt", header = FALSE, row.names = 1))
S17=as.matrix(read.table("200811_I29_V300060997_L4_m120_3_1.fq.gz.bbduk.fq.gz_Aligned.sortedByCoord.out.bam_new.bam_filter.bam.sam_htseq_count_union_nc.txt", header = FALSE, row.names = 1))

col=cbind(S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16,S17)
colnames(col)=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15","S16","S17")
col <- as.data.frame(cbind(rownames(col), col)) # rowname is a metadata and rest all the rows are data so to convert metadata to data one should run this command 
rownames(col) <- NULL
colnames(col) <- c("Gene_id","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15","S16","S17")

merged.file = merge(col, annot, by.x='Gene_id', by.y='Gene_ID')
write.table(merged.file, file = "All_annot.txt", row.names = F, sep = "\t", quote = F)

#in awk
#cat All_annot.txt | awk '{counter1[$21] += $2} {counter2[$21] += $3} {counter3[$21] += $4} {counter4[$21] += $5} {counter5[$21] += $6} {counter6[$21] += $7} {counter7[$21] += $8} {counter8[$21] += $9} {counter9[$21] += $10} {counter10[$21] += $11} {counter11[$21] += $12} {counter12[$21] += $13} {counter13[$21] += $14} {counter14[$21] += $15} {counter15[$21] += $16} {counter16[$21] += $17} {counter17[$21] += $18} END {for(i in counter1){print i"\t"counter1[i]"\t"counter2[i]"\t"counter3[i]"\t"counter4[i]"\t"counter5[i]"\t"counter6[i]"\t"counter7[i]"\t"counter8[i]"\t"counter9[i]"\t"counter10[i]"\t"counter11[i]"\t"counter12[i]"\t"counter13[i]"\t"counter14[i]"\t"counter15[i]"\t"counter16[i]"\t"counter17[i]}}' >> Count_annotated_new_RNA_species_count.txt



#Data3 #GRCh37.p13 non coding 
setwd("~/GenomeDK/faststorage/reference/UCSC/human")
annot = read.delim('gencode.v19.long_noncoding_RNAs_sorted_filtered_modified.gtf', header=F)
colnames(annot) <- c("Gene_ID","gene_type","gene_name","transcript_type","transcript_status")

setwd("~/GenomeDK/faststorage/F19FTSEUET0176_HUMzyhE/Clean/All_new")
S1=as.matrix(read.table("V300013149_L4_DKRT190705OligoTB49-67_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1)) #0
S2=as.matrix(read.table("V300013149_L4_DKRT190705OligoTB50-68_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1)) #30
S3=as.matrix(read.table("V300013149_L4_DKRT190705OligoTB51-69_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1)) #60
S4=as.matrix(read.table("V300027153_L3_DKRT191031TruseqTotal1-49_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1)) #0  
S5=as.matrix(read.table("V300027153_L3_DKRT191031TruseqTotal2-50_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1)) #30
S6=as.matrix(read.table("V300027153_L3_DKRT191031TruseqTotal3-51_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1)) #31
S7=as.matrix(read.table("V300027153_L3_DKRT191031TruseqTotal4-52_1.fq.gz_Aligned.sortedByCoord.out.txt", header = FALSE, row.names = 1)) #60

col=cbind(S1,S2,S3,S4,S5,S6,S7)
colnames(col)=c("S1","S2","S3","S4","S5","S6","S7")
col <- as.data.frame(cbind(rownames(col), col)) # rowname is a metadata and rest all the rows are data so to convert metadata to data one should run this command 
rownames(col) <- NULL
colnames(col) <- c("Gene_id","S1","S2","S3","S4","S5","S6","S7")

merged.file = merge(col, annot, by.x='Gene_id', by.y='Gene_ID')
write.table(merged.file, file = "All_annot.txt", row.names = F, sep = "\t", quote = F)

#in awk
#cat All_annot.txt | awk '{counter1[$21] += $2} {counter2[$21] += $3} {counter3[$21] += $4} {counter4[$21] += $5} {counter5[$21] += $6} {counter6[$21] += $7} {counter7[$21] += $8} {counter8[$21] += $9} {counter9[$21] += $10} {counter10[$21] += $11} {counter11[$21] += $12} {counter12[$21] += $13} {counter13[$21] += $14} {counter14[$21] += $15} {counter15[$21] += $16} {counter16[$21] += $17} {counter17[$21] += $18} END {for(i in counter1){print i"\t"counter1[i]"\t"counter2[i]"\t"counter3[i]"\t"counter4[i]"\t"counter5[i]"\t"counter6[i]"\t"counter7[i]"\t"counter8[i]"\t"counter9[i]"\t"counter10[i]"\t"counter11[i]"\t"counter12[i]"\t"counter13[i]"\t"counter14[i]"\t"counter15[i]"\t"counter16[i]"\t"counter17[i]}}' >> Count_annotated_new_RNA_species_count.txt

####Data1 - coding #GRCh37 coding 
setwd("~/GenomeDK/faststorage/reference/UCSC/human")
annot = read.delim('gencode.v36lift37.annotation.filtered_modified.txt', header=F)
colnames(annot) <- c("Gene_ID","gene_type","gene_name")

setwd("~/GenomeDK/faststorage/F20_new/Sym/F20FTSEUHT1552_HUMjmvN/Clean/coding_count_gencode")
S1=as.matrix(read.table("1_1.fq.gz_Aligned.sortedByCoord.out.coding_gencode.txt_new.txt", header = FALSE, row.names = 1))
S2=as.matrix(read.table("2_1.fq.gz_Aligned.sortedByCoord.out.coding_gencode.txt_new.txt", header = FALSE, row.names = 1))
S3=as.matrix(read.table("3_1.fq.gz_Aligned.sortedByCoord.out.coding_gencode.txt_new.txt", header = FALSE, row.names = 1))
S4=as.matrix(read.table("4_1.fq.gz_Aligned.sortedByCoord.out.coding_gencode.txt_new.txt", header = FALSE, row.names = 1))
S5=as.matrix(read.table("5_1.fq.gz_Aligned.sortedByCoord.out.coding_gencode.txt_new.txt", header = FALSE, row.names = 1))
S6=as.matrix(read.table("6_1.fq.gz_Aligned.sortedByCoord.out.coding_gencode.txt_new.txt", header = FALSE, row.names = 1))
S7=as.matrix(read.table("7_1.fq.gz_Aligned.sortedByCoord.out.coding_gencode.txt_new.txt", header = FALSE, row.names = 1))
S8=as.matrix(read.table("8_1.fq.gz_Aligned.sortedByCoord.out.coding_gencode.txt_new.txt", header = FALSE, row.names = 1))
S9=as.matrix(read.table("9_1.fq.gz_Aligned.sortedByCoord.out.coding_gencode.txt_new.txt", header = FALSE, row.names = 1))
S10=as.matrix(read.table("10_1.fq.gz_Aligned.sortedByCoord.out.coding_gencode.txt_new.txt", header = FALSE, row.names = 1))
S11=as.matrix(read.table("11_1.fq.gz_Aligned.sortedByCoord.out.coding_gencode.txt_new.txt", header = FALSE, row.names = 1))
S12=as.matrix(read.table("12_1.fq.gz_Aligned.sortedByCoord.out.coding_gencode.txt_new.txt", header = FALSE, row.names = 1))

col=cbind(S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12)
colnames(col)=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12")
col <- as.data.frame(cbind(rownames(col), col)) # rowname is a metadata and rest all the rows are data so to convert metadata to data one should run this command 
rownames(col) <- NULL
colnames(col) <- c("Gene_id","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12")

merged.file = merge(col, annot, by.x='Gene_id', by.y='Gene_ID')
write.table(merged.file, file = "All_annot.txt", row.names = F, sep = "\t", quote = F)

#in awk
#cat All_annot.txt | awk '{counter1[$14] += $2} {counter2[$16] += $3} {counter3[$16] += $4} {counter4[$16] += $5} {counter5[$16] += $6} {counter6[$16] += $7} {counter7[$16] += $8} {counter8[$16] += $9} {counter9[$16] += $10} {counter10[$16] += $11} {counter11[$16] += $12} {counter12[$16] += $13} END {for(i in counter1){print i"\t"counter1[i]"\t"counter2[i]"\t"counter3[i]"\t"counter4[i]"\t"counter5[i]"\t"counter6[i]"\t"counter7[i]"\t"counter8[i]"\t"counter9[i]"\t"counter10[i]"\t"counter11[i]"\t"counter12[i]}}' >> Count_annotated_new_RNA_species_count.txt 

####Data1 - coding #GRCh37 noncoding 
setwd("~/GenomeDK/faststorage/reference/UCSC/human/GRCh37")
annot = read.delim('gencode.v36lift37.long_noncoding_RNAs_filtered_modified.gtf', header=F)
colnames(annot) <- c("Gene_ID","gene_type","gene_name")

setwd("~/GenomeDK/faststorage/F20_new/Sym/F20FTSEUHT1552_HUMjmvN/Clean/03_Noncoding_count_GRCh37")
S1=as.matrix(read.table("1_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.NC.txt", header = FALSE, row.names = 1))
S2=as.matrix(read.table("2_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.NC.txt", header = FALSE, row.names = 1))
S3=as.matrix(read.table("3_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.NC.txt", header = FALSE, row.names = 1))
S4=as.matrix(read.table("4_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.NC.txt", header = FALSE, row.names = 1))
S5=as.matrix(read.table("5_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.NC.txt", header = FALSE, row.names = 1))
S6=as.matrix(read.table("6_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.NC.txt", header = FALSE, row.names = 1))
S7=as.matrix(read.table("7_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.NC.txt", header = FALSE, row.names = 1))
S8=as.matrix(read.table("8_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.NC.txt", header = FALSE, row.names = 1))
S9=as.matrix(read.table("9_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.NC.txt", header = FALSE, row.names = 1))
S10=as.matrix(read.table("10_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.NC.txt", header = FALSE, row.names = 1))
S11=as.matrix(read.table("11_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.NC.txt", header = FALSE, row.names = 1))
S12=as.matrix(read.table("12_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.NC.txt", header = FALSE, row.names = 1))

col=cbind(S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12)
colnames(col)=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12")
col <- as.data.frame(cbind(rownames(col), col)) # rowname is a metadata and rest all the rows are data so to convert metadata to data one should run this command 
rownames(col) <- NULL
colnames(col) <- c("Gene_id","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12")

merged.file = merge(col, annot, by.x='Gene_id', by.y='Gene_ID')
write.table(merged.file, file = "All_annot.txt", row.names = F, sep = "\t", quote = F)

#in awk
#cat All_annot.txt | awk '{counter1[$14] += $2} {counter2[$16] += $3} {counter3[$16] += $4} {counter4[$16] += $5} {counter5[$16] += $6} {counter6[$16] += $7} {counter7[$16] += $8} {counter8[$16] += $9} {counter9[$16] += $10} {counter10[$16] += $11} {counter11[$16] += $12} {counter12[$16] += $13} END {for(i in counter1){print i"\t"counter1[i]"\t"counter2[i]"\t"counter3[i]"\t"counter4[i]"\t"counter5[i]"\t"counter6[i]"\t"counter7[i]"\t"counter8[i]"\t"counter9[i]"\t"counter10[i]"\t"counter11[i]"\t"counter12[i]}}' >> Count_annotated_new_RNA_species_count.txt 

#Data1 #GRCh37.p13 coding 
setwd("~/GenomeDK/faststorage/reference/UCSC/human/GRCh37.p13")
annot = read.delim('gencode.v19.chr_patch_hapl_scaff.annotation_filtered_modified.gtf', header=F)
colnames(annot) <- c("Gene_ID","gene_type","gene_name","transcript_type","transcript_status")

setwd("~/GenomeDK/faststorage/F20_new/Sym/F20FTSEUHT1552_HUMjmvN/Clean/01_C_count_GRCh37.p13")
S1=as.matrix(read.table("1_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.p13C.txt", header = FALSE, row.names = 1))
S2=as.matrix(read.table("2_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.p13C.txt", header = FALSE, row.names = 1))
S3=as.matrix(read.table("3_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.p13C.txt", header = FALSE, row.names = 1))
S4=as.matrix(read.table("4_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.p13C.txt", header = FALSE, row.names = 1))
S5=as.matrix(read.table("5_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.p13C.txt", header = FALSE, row.names = 1))
S6=as.matrix(read.table("6_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.p13C.txt", header = FALSE, row.names = 1))
S7=as.matrix(read.table("7_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.p13C.txt", header = FALSE, row.names = 1))
S8=as.matrix(read.table("8_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.p13C.txt", header = FALSE, row.names = 1))
S9=as.matrix(read.table("9_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.p13C.txt", header = FALSE, row.names = 1))
S10=as.matrix(read.table("10_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.p13C.txt", header = FALSE, row.names = 1))
S11=as.matrix(read.table("11_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.p13C.txt", header = FALSE, row.names = 1))
S12=as.matrix(read.table("12_1.fq.gz_Aligned.sortedByCoord.out.GRCh37.p13C.txt", header = FALSE, row.names = 1))

col=cbind(S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12)
colnames(col)=c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12")
col <- as.data.frame(cbind(rownames(col), col)) # rowname is a metadata and rest all the rows are data so to convert metadata to data one should run this command 
rownames(col) <- NULL
colnames(col) <- c("Gene_id","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12")

merged.file = merge(col, annot, by.x='Gene_id', by.y='Gene_ID')
write.table(merged.file, file = "All_annot.txt", row.names = F, sep = "\t", quote = F)

#in awk
#cat All_annot.txt | awk '{counter1[$16] += $2} {counter2[$16] += $3} {counter3[$16] += $4} {counter4[$16] += $5} {counter5[$16] += $6} {counter6[$16] += $7} {counter7[$16] += $8} {counter8[$16] += $9} {counter9[$16] += $10} {counter10[$16] += $11} {counter11[$16] += $12} {counter12[$16] += $13} END {for(i in counter1){print i"\t"counter1[i]"\t"counter2[i]"\t"counter3[i]"\t"counter4[i]"\t"counter5[i]"\t"counter6[i]"\t"counter7[i]"\t"counter8[i]"\t"counter9[i]"\t"counter10[i]"\t"counter11[i]"\t"counter12[i]}}' >> Count_annotated_new_RNA_species_count.txt 


#Data1 #GRCh37.p13 all
setwd("~/GenomeDK/faststorage/reference/UCSC/mm9/ens")
annot = read.delim('ens_genename_ID_name.gtf', header=T)

annot1 = read.delim('3UTR_ens_new.bed', header=F)
colnames(annot1) <- c("Chr","Start","End","ens_id")
merged.file = merge(annot1, annot, by.x='ens_id', by.y='mm9.ensemblToGeneName.name')
write.table(merged.file, file = "3UTR_ens_annot.bed", row.names = F, sep = "\t", quote = F)

annot2 = read.delim('5UTR_ens_new.bed', header=F)
colnames(annot2) <- c("Chr","Start","End","ens_id")
merged.file1 = merge(annot2, annot, by.x='ens_id', by.y='mm9.ensemblToGeneName.name')
write.table(merged.file1, file = "5UTR_ens_annot.bed", row.names = F, sep = "\t", quote = F)

annot3 = read.delim('Intron_ens_new.bed', header=F)
colnames(annot3) <- c("Chr","Start","End","ens_id")
merged.file2 = merge(annot3, annot, by.x='ens_id', by.y='mm9.ensemblToGeneName.name')
write.table(merged.file2, file = "Intron_ens_annot.bed", row.names = F, sep = "\t", quote = F)

annot4 = read.delim('CodingExon_ens_new.bed', header=F)
colnames(annot4) <- c("Chr","Start","End","ens_id")
merged.file3 = merge(annot4, annot, by.x='ens_id', by.y='mm9.ensemblToGeneName.name')
write.table(merged.file3, file = "CodingExon_ens_annot.bed", row.names = F, sep = "\t", quote = F)

annot5 = read.delim('Exon_ens_new.bed', header=F)
colnames(annot5) <- c("Chr","Start","End","ens_id")
merged.file4 = merge(annot5, annot, by.x='ens_id', by.y='mm9.ensemblToGeneName.name')
write.table(merged.file4, file = "Exon_ens_annot.bed", row.names = F, sep = "\t", quote = F)

#Overall normalization technique
#background corrction
#there are two methods to normalize background correction however there are many other background correction option
MA <- normalizeWithinArrays(RG) #1  
RGb <- backgroundCorrect(RG, method="subtract") #2 is the same as above
MA <- normalizeWithinArrays(RGb)#2
RG <- backgroundCorrect(RG, method="normexp", offset=50) #3 for the purpose of DE

#withinarrary normalization
MA <- normalizeWithinArrays(RG) #1
MA <- normalizeWithinArrays(RG, method="loess") #2
MA <- normalizeWithinArrays(RG, method="robustspline") #3
MA <- normalizeWithinArrays(RG, weights=NULL) #4

#between array normalization 
MA.pAq <- normalizeBetweenArrays(MA.p, method="Aquantile") #1 Aquantile
plotDensities(MA.pAq)

MA.q <- normalizeBetweenArrays(RG.b, method="quantile") #2 Quantile
plotDensities(MA.q, col="black")

MA.q <- normalizeBetweenArrays(RG.b, method="vsn") #3 variance stabilizing normalization 
plotDensities(MA.q, col="black")

# All of the step in combination
MA.p <-normalizeWithinArrays(RG.b)
plotDensities(MA.p)

#Spike-in analysis using RUV
setwd("~/GenomeDK/faststorage/reference/UCSC/human")
annot = read.delim('gencode.v19.long_noncoding_RNAs_sorted_filtered_modified.gtf', header=F)
colnames(annot) <- c("Gene_ID","gene_type","gene_name","transcript_type","transcript_status")

#for file in /home/anksi/faststorage/F20_new/Sym/F20FTSEUHT1552_HUMjmvN/Clean/ERCC/*ERCC.SU.txt; do
 # cat $file | sort | uniq >> ${file/.txt/.SU.txt}
#done

#for file in /home/anksi/faststorage/F20_new/Sym/F20FTSEUHT1552_HUMjmvN/Clean/*/*ERCC.coding_gencode.txt; do
# cp $file /home/anksi/faststorage/F20_new/Sym/F20FTSEUHT1552_HUMjmvN/Clean/ERCC_new
#done

#for file in /home/anksi/faststorage/F20_new/Sym/F20FTSEUHT1552_HUMjmvN/Clean/ERCC/*ERCC.SU.txt; do
#  grep -v ^__ $file >> ${file/.ERCC.SU.txt/.ERCC.SU.new.txt}
#done

#for file in /home/anksi/faststorage/F20_new/Sym/F20FTSEUHT1552_HUMjmvN/Clean/ERCC_new/*.txt; do
#  grep -v ^__ $file >> ${file/.txt/.new.txt}
#done

setwd("~/GenomeDK/faststorage/F20_new/Sym/F20FTSEUHT1h552_HUMjmvN/Clean/ERCC_new")
S1=as.matrix(read.table("1_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S2=as.matrix(read.table("2_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S3=as.matrix(read.table("3_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S4=as.matrix(read.table("4_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S5=as.matrix(read.table("5_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S6=as.matrix(read.table("6_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S7=as.matrix(read.table("7_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S8=as.matrix(read.table("8_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S9=as.matrix(read.table("9_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S10=as.matrix(read.table("10_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S11=as.matrix(read.table("11_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S12=as.matrix(read.table("12_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))

##
col=cbind(S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12)
#for nuclear/nucleolar NAT10 vs WT
col=cbind(S1,S2,S3,S4,S5,S6)
col=cbind(S7,S8,S9,S10,S11,S12)

#for nuclear WT vs nucleolar WT /for nuclear NAT10 vs nucleolar NAT10
col=cbind(S1,S2,S3,S7,S8,S9)
col=cbind(S4,S5,S6,S10,S11,S12)

#colnames(col)=c("WT_O1","NAT10_O1","WT_i1","NAT10_i1","WT_O2","NAT10_O2","WT_i2","NAT10_i2","WT_O3","NAT10_O3","WT_i3","NAT10_i3")

colnames(col)=c("WT_O1","WT_O2","WT_O3","NAT10_O1","NAT10_O2","NAT10_O3","WT_i1","WT_i2","WT_i3","NAT10_i1","NAT10_i2","NAT10_i3")

#for nuclear/nucleolar NAT10 vs WT
colnames(col)=c("WT_O1","WT_O2","WT_O3","NAT10_O1","NAT10_O2","NAT10_O3")
colnames(col)=c("WT_i1","WT_i2","WT_i3","NAT10_i1","NAT10_i2","NAT10_i3")

#for nuclear WT vs nucleolar WT /for nuclear NAT10 vs nucleolar NAT10
colnames(col)=c("WT_O1","WT_O2","WT_O3","WT_i1","WT_i2","WT_i3")
colnames(col)=c("NAT10_O1","NAT10_O2","NAT10_O3","NAT10_i1","NAT10_i2","NAT10_i3")

#colnames(col)=c("WT1","NAT101","WT2","NAT102","WT3","NAT103","WT4","NAT104","WT5","NAT105","WT6","NAT106")
#col <- as.data.frame(cbind(rownames(col), col)) # rowname is a metadata and rest all the rows are data so to convert metadata to data one should run this command 
#rownames(col) <- NULL
#colnames(col) <- c("Gene_id","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S12")
filter <- apply(col, 1, function(x) length(x[x>5])>=2)
filter2 <- apply(col, 1, function(x) length(x[x>1])>=2)
filtered <- col[filter,]
filtered2 <- col[filter2,]
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]

#x <- as.factor(rep(c("WT", "NAT10"), each=3))
#x <- as.factor(c("WT_O","NAT10_O","WT_i1","NAT10_i1","WT_O2","NAT10_O2","WT_i2","NAT10_i2","WT_O3","NAT10_O3","WT_i3","NAT10_i3"))

x <- as.factor(rep(c("WTO","NAT10O","WTi","NAT10i"), each=6))
x <- as.factor(rep(c("WTO","NAT10O","WTi","NAT10i"), each=3))
#for nuclear/nucleolar NAT10 vs WT
x <- as.factor(rep(c("WTO","NAT10O"), each=3))
x <- as.factor(rep(c("WTi","NAT10i"), each=3))

#for nuclear WT vs nucleolar WT /for nuclear NAT10 vs nucleolar NAT10
x <- as.factor(rep(c("WTO","WTi"), each=3))
x <- as.factor(rep(c("NAT10O","NAT10i"), each=3))

BiocManager::install("EDASeq")
BiocManager::install("parallel")
library(parallel)
library(BiocGenerics)
library(Biobase)
library(EDASeq)
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered)))
set

library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

set0 <- betweenLaneNormalization(set, which="upper")
plotRLE(set0, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set0, col=colors[x], cex=1.2)

BiocManager::install("RUVSeq")
library(RUVSeq)
set1 <- RUVg(set0, spikes, k=1)
pData(set1)
plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set1, col=colors[x], cex=1.2)

head(normCounts(set))

head(normCounts(set0))
set0_normCount <- normCounts(set0)
set0_normCount_new <- cbind(rownames(set0_normCount), set0_normCount)
rownames(set0_normCount_new) <- NULL
colnames(set0_normCount_new) <- c("SYMBOL", "WT_O1","WT_O2","WT_O3","NAT10_O1","NAT10_O2","NAT10_O3","WT_i1","WT_i2","WT_i3","NAT10_i1","NAT10_i2","NAT10_i3")
write.table(set0_normCount_new, file = "ERCC.V19.C_nor1.txt", row.names = F, sep = "\t", quote = F)

head(normCounts(set1))
set1_normCount <- normCounts(set1)
set1_normCount_new <- cbind(rownames(set1_normCount), set1_normCount)
rownames(set1_normCount_new) <- NULL
colnames(set1_normCount_new) <- c("SYMBOL", "WT_O1","WT_O2","WT_O3","NAT10_O1","NAT10_O2","NAT10_O3","WT_i1","WT_i2","WT_i3","NAT10_i1","NAT10_i2","NAT10_i3")
write.table(set1_normCount_new, file = "ERCC.V19.C_nor2.txt", row.names = F, sep = "\t", quote = F)

head(counts(set1))

#DE edgeR
design <- model.matrix(~x + W_1, data=pData(set1))
y <- DGEList(counts=counts(set1), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)
summary(dt<-decideTestsDGE(lrt, p.value = 0.05))
deg <- topTags(lrt, n = Inf, p = 0.05)$table

up <- deg[deg$logFC > 0,]
down <- deg[deg$logFC < 0,]
nrow(lrt)

#while comparing WT and NAT10 we found 8 genes to be regulated and 12 genes to be downregulated

# while comparing WTi vs WTo, none of the genes are found to be sig up/down regulated 
Wt_i_o <- cbind(rownames(lrt$table),lrt$table)
colnames(Wt_i_o) <- c("SYMBOL","logFC","logCPM","LR","PValue")                 
write.table(Wt_i_o, file = "WT_ovsI_entireTable_new.txt", row.names = F, sep = "\t", quote = F)

NAT10_i_o <- cbind(rownames(lrt$table),lrt$table)
colnames(NAT10_i_o) <- c("SYMBOL","logFC","logCPM","LR","PValue")                 
write.table(NAT10_i_o, file = "NAT10_ovsI_entireTable_new_new.txt", row.names = F, sep = "\t", quote = F)

up <- cbind(rownames(up), up)
rownames(up) <- NULL
colnames(up) <- c("SYMBOL", "logFC", "logCPM", "LR","PValue","FDR")
write.table(up, file = "up_nucleolous.txt", row.names = F, sep = "\t", quote = F)

down <- cbind(rownames(down), down)
rownames(down) <- NULL
colnames(up) <- c("SYMBOL", "logFC", "logCPM", "LR","PValue","FDR")
write.table(down, file = "down_nucleolous.txt", row.names = F, sep = "\t", quote = F)

up <- cbind(rownames(up), up)
rownames(up) <- NULL
colnames(up) <- c("SYMBOL", "logFC", "logCPM", "LR","PValue","FDR")
write.table(up, file = "up_nucleus.txt", row.names = F, sep = "\t", quote = F)

down <- cbind(rownames(down), down)
rownames(down) <- NULL
colnames(up) <- c("SYMBOL", "logFC", "logCPM", "LR","PValue","FDR")
write.table(down, file = "down_nucleus.txt", row.names = F, sep = "\t", quote = F)


#DE deseq2
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts(set1),
                              colData = pData(set1),
                              design = ~ W_1 + x)
dds <- DESeq(dds)
res <- results(dds)
res

dds <- DESeq(dds, test="LRT", reduced=as.formula("~ W_1"))
res <- results(dds)
res

#DE using replicate samples 
differences <- makeGroups(x)
differences
set2 <- RUVs(set0, genes, k=1, differences)
pData(set2)

#It is true that glmLRT() does not allow you to specify multiple contrast vectors to get an F-like test.
#thats why I have done it separately for nucleous and nucleolous
#DE using residual
design <- model.matrix(~x, data=pData(set0))
y1 <- DGEList(counts=counts(set0), group=x)
y1 <- calcNormFactors(y, method="upperquartile")
y1 <- estimateGLMCommonDisp(y1, design)
y1 <- estimateGLMTagwiseDisp(y1, design)
fit <- glmFit(y1, design)
res <- residuals(fit, type="deviance")

set3 <- RUVr(set0, genes, k=1, res)
pData(set3)

#this one is redone so as to compare WT with NAT10
setwd("~/GenomeDK/faststorage/F20_new/Sym/F20FTSEUHT1552_HUMjmvN/Clean/ERCC.V19.Cod")
S1=as.matrix(read.table("1_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #WO
S2=as.matrix(read.table("5_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #WO
S3=as.matrix(read.table("9_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #WO
S4=as.matrix(read.table("3_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #Wi
S5=as.matrix(read.table("7_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #Wi
S6=as.matrix(read.table("11_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #Wi
S7=as.matrix(read.table("2_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #NO
S8=as.matrix(read.table("6_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #NO
S9=as.matrix(read.table("10_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #NO
S10=as.matrix(read.table("4_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #Ni
S11=as.matrix(read.table("8_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #Ni
S12=as.matrix(read.table("12_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #Ni
colnames(col)=c("WT_O1","WT_O2","WT_O3","WT_i1","WT_i2","WT_i3","NAT10_O1","NAT10_O2","NAT10_O3","NAT10_i1","NAT10_i2","NAT10_i3")
x <- as.factor(rep(c("WTOi","NAT10Oi"), each=6))
#rough
S1=as.matrix(read.table("1_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #WO
S2=as.matrix(read.table("5_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #WO
S3=as.matrix(read.table("9_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #WO
S4=as.matrix(read.table("2_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #NO
S5=as.matrix(read.table("6_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #NO
S6=as.matrix(read.table("10_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #NO
S7=as.matrix(read.table("3_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #Wi
S8=as.matrix(read.table("7_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #Wi
S9=as.matrix(read.table("11_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #Wi
S10=as.matrix(read.table("4_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #Ni
S11=as.matrix(read.table("8_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #Ni
S12=as.matrix(read.table("12_1.fq.gz_Aligned.sortedByCoord.out.V19.C.ERCC.new.txt", header = FALSE, row.names = 1)) #Ni



setwd("~/GenomeDK/faststorage/F20_new/Sym/F20FTSEUHT1552_HUMjmvN/Clean/ERCC.V19.NC")
S1=as.matrix(read.table("1_1.fq.gz_Aligned.sortedByCoord.outV19.NC.new.txt", header = FALSE, row.names = 1)) #WO
S2=as.matrix(read.table("5_1.fq.gz_Aligned.sortedByCoord.outV19.NC.new.txt", header = FALSE, row.names = 1)) #WO
S3=as.matrix(read.table("9_1.fq.gz_Aligned.sortedByCoord.outV19.NC.new.txt", header = FALSE, row.names = 1)) #WO
S4=as.matrix(read.table("2_1.fq.gz_Aligned.sortedByCoord.outV19.NC.new.txt", header = FALSE, row.names = 1)) #NO
S5=as.matrix(read.table("6_1.fq.gz_Aligned.sortedByCoord.outV19.NC.new.txt", header = FALSE, row.names = 1)) #NO
S6=as.matrix(read.table("10_1.fq.gz_Aligned.sortedByCoord.outV19.NC.new.txt", header = FALSE, row.names = 1)) #NO
S7=as.matrix(read.table("3_1.fq.gz_Aligned.sortedByCoord.outV19.NC.new.txt", header = FALSE, row.names = 1)) #Wi
S8=as.matrix(read.table("7_1.fq.gz_Aligned.sortedByCoord.outV19.NC.new.txt", header = FALSE, row.names = 1)) #Wi
S9=as.matrix(read.table("11_1.fq.gz_Aligned.sortedByCoord.outV19.NC.new.txt", header = FALSE, row.names = 1)) #Wi
S10=as.matrix(read.table("4_1.fq.gz_Aligned.sortedByCoord.outV19.NC.new.txt", header = FALSE, row.names = 1)) #Ni
S11=as.matrix(read.table("8_1.fq.gz_Aligned.sortedByCoord.outV19.NC.new.txt", header = FALSE, row.names = 1)) #Ni
S12=as.matrix(read.table("12_1.fq.gz_Aligned.sortedByCoord.outV19.NC.new.txt", header = FALSE, row.names = 1)) #Ni

setwd("~/GenomeDK/faststorage/F20_new/Sym/F20FTSEUHT1552_HUMjmvN/Clean/ERCC.GRCh37.Cod")
S1=as.matrix(read.table("1_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S2=as.matrix(read.table("5_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S3=as.matrix(read.table("9_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S4=as.matrix(read.table("2_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S5=as.matrix(read.table("6_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S6=as.matrix(read.table("10_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S7=as.matrix(read.table("3_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S8=as.matrix(read.table("7_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S9=as.matrix(read.table("11_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S10=as.matrix(read.table("4_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S11=as.matrix(read.table("8_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))
S12=as.matrix(read.table("12_1.fq.gz_Aligned.sortedByCoord.out.ERCC.coding_gencode.new.txt", header = FALSE, row.names = 1))

setwd("~/GenomeDK/faststorage/F20_new/Sym/F20FTSEUHT1552_HUMjmvN/Clean/ERCC.GRCh37.NC")
S1=as.matrix(read.table("1_1.fq.gz_Aligned.sortedByCoord.out.ERCC.GRCh37.NC.new.txt", header = FALSE, row.names = 1))
S2=as.matrix(read.table("5_1.fq.gz_Aligned.sortedByCoord.out.ERCC.GRCh37.NC.new.txt", header = FALSE, row.names = 1))
S3=as.matrix(read.table("9_1.fq.gz_Aligned.sortedByCoord.out.ERCC.GRCh37.NC.new.txt", header = FALSE, row.names = 1))
S4=as.matrix(read.table("2_1.fq.gz_Aligned.sortedByCoord.out.ERCC.GRCh37.NC.new.txt", header = FALSE, row.names = 1))
S5=as.matrix(read.table("6_1.fq.gz_Aligned.sortedByCoord.out.ERCC.GRCh37.NC.new.txt", header = FALSE, row.names = 1))
S6=as.matrix(read.table("10_1.fq.gz_Aligned.sortedByCoord.out.ERCC.GRCh37.NC.new.txt", header = FALSE, row.names = 1))
S7=as.matrix(read.table("3_1.fq.gz_Aligned.sortedByCoord.out.ERCC.GRCh37.NC.new.txt", header = FALSE, row.names = 1))
S8=as.matrix(read.table("7_1.fq.gz_Aligned.sortedByCoord.out.ERCC.GRCh37.NC.new.txt", header = FALSE, row.names = 1))
S9=as.matrix(read.table("11_1.fq.gz_Aligned.sortedByCoord.out.ERCC.GRCh37.NC.new.txt", header = FALSE, row.names = 1))
S10=as.matrix(read.table("4_1.fq.gz_Aligned.sortedByCoord.out.ERCC.GRCh37.NC.new.txt", header = FALSE, row.names = 1))
S11=as.matrix(read.table("8_1.fq.gz_Aligned.sortedByCoord.out.ERCC.GRCh37.NC.new.txt", header = FALSE, row.names = 1))
S12=as.matrix(read.table("12_1.fq.gz_Aligned.sortedByCoord.out.ERCC.GRCh37.NC.new.txt", header = FALSE, row.names = 1))


##########################
############## corelation analysis ##############

install.packages("ggpubr")
library("ggpubr")
my_data1 <- read.delim("try_sort_ERCC.V19.C_nor2_annot_sorted_unique3.txt", header = FALSE)
colnames(my_data1) <- c("SYMBOL", "enrichment")
my_data2 <- read.delim("ERCC.V19.C_nor2_annot_sorted_unique_try_sorted3.txt", header = FALSE)
#my_data2 <- cbind(rownames(my_data2), my_data2)
colnames(my_data2) <- c("SYM","WT_O1","WT_O2","WT_O3","NAT10_O1","NAT10_O2","NAT10_O3","WT_i1","WT_i2","WT_i3","NAT10_i1","NAT10_i2","gene_id", "V1", "Gene", "V2","status")

rownames(my_data1) <- my_data1$SYMBOL
rownames(my_data2) <- my_data2$SYM

my_data = merge(my_data1, my_data2, by.x='SYMBOL', by.y='SYM')
rownames(my_data) <- my_data$SYMBOL

ggscatter(my_data, x = "enrichment", y = "NAT10_i1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "enrichment", ylab = "WT_01")

ggscatter(my_data, x = "NAT10_O1", y = "NAT10_i1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "NAT10_O1", ylab = "NAT10_i1")

shapiro.test(my_data$enrichment)
shapiro.test(my_data$WT_O3)

ggqqplot(my_data$enrichment, ylab = "enrichment")
ggqqplot(my_data$NAT10_i1, ylab = "WT_01")

res <- cor.test(my_data$enrichment, my_data$WT_O1, method = "pearson")
res

res <- cor.test(my_data$NAT10_O1, my_data$NAT10_i1, method = "pearson")
res

#merging various replicates 
my_data2_new <- my_data2[1:13]

results=data.frame(apply(array(as.matrix(my_data2_new[,-1]), c(nrow(my_data2_new),3, ncol(my_data2_new)/3)),3, rowMeans))
results=cbind(my_data2_new$SYM, results)
head(results)
colnames(results) <- c("Sym", "WT_O", "NAT10_O", "WT_i", "NAT10_i")
rownames(results) <- results$Sym

#qqplot to see whether the data is normally distributed or not
# if the data is normally distributed then perform parametric test (pearson correlation) else perform non-parametric test (spearman correlation) 

#does the data follow normal distribution?
shapiro.test(results$WT_O)
#W = 0.055676, p-value < 2.2e-16
#null hypothesis is rejected and hence the data is not normally distributed 
shapiro.test(results$WT_i)
#W = 0.022886, p-value < 2.2e-16
#null hypothesis is rejected and hence the data is not normally distributed 
shapiro.test(results$NAT10_O)
#W = 0.10793, p-value < 2.2e-16
#null hypothesis is rejected and hence the data is not normally distributed 
shapiro.test(results$NAT10_i)
#W = 0.15093, p-value < 2.2e-16
#null hypothesis is rejected and hence the data is not normally distributed 

#Visual inspection of if the data is normally distributed?
par(mfrow=c(2,2))
g1 <- ggqqplot(results$WT_O, ylab = "WT_O")
g2 <-ggqqplot(results$WT_i, ylab = "WT_i")
g3 <-ggqqplot(results$NAT10_O, ylab = "NAT10_O")
g4 <- ggqqplot(results$NAT10_i, ylab = "NAT10_i")

ggarrange(g1, g2, g3, g4, ncol=2, nrow =2)

# summary: none of the data is normally distributed

#Is the covarion linear or not?
#install.packages("UsingR")
library(UsingR) # for dividing the page into four sections #did not work
par(mfrow=c(2,2)) # did not work
par(mfrow=c(2,2)) # did not work
p1 <- ggscatter(results, x = "WT_O", y = "WT_i", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "WT_O", ylab = "WT_i")

p2 <- ggscatter(results, x = "NAT10_O", y = "NAT10_i", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "NAT10_O", ylab = "NAT10_i")

p3 <- ggscatter(results, x = "WT_O", y = "NAT10_O", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "WT_O", ylab = "NAT10_O")

p4 <- ggscatter(results, x = "WT_i", y = "NAT10_i", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "WT_i", ylab = "NAT10_i")

grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2) # does not work
ggarrange(p1, p2, p3, p4, ncol=2, nrow =2)
#Yes, from the above scatter plot it seems that the covarion is linear

#correlation analysis
res1 <- cor.test(results$NAT10_O, results$NAT10_i, method = "pearson")
res2 <- cor.test(results$NAT10_O, results$NAT10_i, method = "kendall")

res3_1 <- cor.test(results$NAT10_O, results$NAT10_i, method = "spearman", exact = FALSE)
res3_2 <- cor.test(results$WT_O, results$WT_i, method = "spearman", exact = FALSE)
res3_3 <- cor.test(results$WT_O, results$NAT10_O, method = "spearman", exact = FALSE)
res3_4 <- cor.test(results$WT_i, results$NAT10_i, method = "spearman", exact = FALSE)

res3_1
res3_2
res3_3
res3_4

result_new <- results[2:5]
corelMat <- cor(result_new, method = c("spearman"))
head(corelMat)
#library(corrplot)
par(mfrow=c(1,1))
corrplot(corelMat)
################################################

#enrichment analysis for coding 
s1=as.matrix(read.table("ac4c0_peaks.narrowPeak_sorted.narrowPeak_corel_V19_NC_new.SU.count.final.txt", header = FALSE, row.names = 1))
s2=as.matrix(read.table("ac4c120_set1_new_peaks.narrowPeak_sorted.narrowPeak_corel_V19_cod_new.SU.count.final.txt", header = FALSE, row.names = 1))
s3=as.matrix(read.table("ac4c120_set2_new_peaks.narrowPeak_sorted.narrowPeak_corel_V19_cod_new.SU.count.final.txt", header = FALSE, row.names = 1))
s4=as.matrix(read.table("ac4c120_set3_new_peaks.narrowPeak_sorted.narrowPeak_corel_V19_cod_new.SU.count.final.txt", header = FALSE, row.names = 1))
s5=as.matrix(read.table("m6a0_set1_new_peaks.narrowPeak_sorted.narrowPeak_corel_V19_cod_new.SU.count.final.txt", header = FALSE, row.names = 1))
s6=as.matrix(read.table("m6a0_set2_new_peaks.narrowPeak_sorted.narrowPeak_corel_V19_cod_new.SU.count.final.txt", header = FALSE, row.names = 1))
s7=as.matrix(read.table("m6a120_set1_new_peaks.narrowPeak_sorted.narrowPeak_corel_V19_cod_new.SU.count.final.txt", header = FALSE, row.names = 1))
s8=as.matrix(read.table("m6a120_set2_new_peaks.narrowPeak_sorted.narrowPeak_corel_V19_cod_new.SU.count.final.txt", header = FALSE, row.names = 1))
s9=as.matrix(read.table("m6a120_set3_new_peaks.narrowPeak_sorted.narrowPeak_corel_V19_cod_new.SU.count.final.txt", header = FALSE, row.names = 1))

#enrichment analysis for non coding 
s1=as.matrix(read.table("ac4c0_peaks.narrowPeak_sorted.narrowPeak_corel_V19_NC_new.SU.count.txt", header = FALSE, row.names = 1))
s2=as.matrix(read.table("ac4c120_set1_new_new_peaks.narrowPeak_sorted.narrowPeak_corel_V19_NC_new.SU.count.txt", header = FALSE, row.names = 1))
s3=as.matrix(read.table("ac4c120_set2_new_peaks.narrowPeak_sorted.narrowPeak_corel_V19_NC_new.SU.count.txt", header = FALSE, row.names = 1))
s4=as.matrix(read.table("ac4c120_set3_new_peaks.narrowPeak_sorted.narrowPeak_corel_V19_NC_new.SU.count.txt", header = FALSE, row.names = 1))
s5=as.matrix(read.table("m6a0_set1_new_peaks.narrowPeak_sorted.narrowPeak_corel_V19_NC_new.SU.count.txt", header = FALSE, row.names = 1))
s6=as.matrix(read.table("m6a0_set2_new_peaks.narrowPeak_sorted.narrowPeak_corel_V19_NC_new.SU.count.txt", header = FALSE, row.names = 1))
s7=as.matrix(read.table("m6a120_set1_new_peaks.narrowPeak_sorted.narrowPeak_corel_V19_NC_new.SU.count.txt", header = FALSE, row.names = 1))
s8=as.matrix(read.table("m6a120_set2_new_peaks.narrowPeak_sorted.narrowPeak_corel_V19_NC_new.SU.count.txt", header = FALSE, row.names = 1))
s9=as.matrix(read.table("m6a120_set3_new_peaks.narrowPeak_sorted.narrowPeak_corel_V19_NC_new.SU.count.txt", header = FALSE, row.names = 1))

#enrichment analysis after fetchinh the read count of the enriched peaks
#setwd("~/GenomeDK/faststorage/chromatin_retension/All/07_bamfilter2")
s1_enrich_count=as.matrix(read.table("counts_m6a120.cnt.ERCC.V19cod_SU_uniq.bed", header = FALSE, row.names = 1))
s2_enrich_count=as.matrix(read.table("counts_ac4c120.cnt.ERCC.V19cod_SU_uniq.bed", header = FALSE, row.names = 1))
#s3=as.matrix(read.table("counts_m6a0.cnt.ERCC.V19cod_SU_uniq.bed", header = FALSE, row.names = 1))
s3_enrich_count=as.matrix(read.table("counts_ac4c0.cnt.ERCC.V19cod_SU_uniq.bed", header = FALSE, row.names = 1))

setwd("~/GenomeDK/faststorage/F19FTSEUET0176_HUMzyhE/Clean/All_new")
#RNA seq data at 0,30 and 60 time points 
my_data3 <- read.delim("All_annot.txt", header = FALSE)
my_data4 <- read.delim("All_annot_repAdded_duprem.txt", header = FALSE)

s1_enrich_count <- cbind(rownames(s1_enrich_count), s1_enrich_count)
colnames(s1_enrich_count) <- c("Sym", "enrich")
#s1_merged = merge(s1_enrich_count, results, by.x='Sym', by.y='Sym')
s1_merged = merge(s1_enrich_count, my_data3, by.x='Sym', by.y='V1')

s1_merged = merge(s1_enrich_count, my_data4, by.x='Sym', by.y='V1')

s1_merged$enrich <- as.numeric(s1_merged$enrich)
s1_merged$V2 <- as.numeric(s1_merged$V2)
s1_merged$V3 <- as.numeric(s1_merged$V3)
s1_merged$V4 <- as.numeric(s1_merged$V4)

ggscatter(s1_merged, x = "enrich", y = "V4", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "enrich", ylab = "V4")
#to remove duplicate rows
s1_merged_new <- s1_merged[!duplicated(s1_merged$Sym), ]
rownames(s1_merged_new) <- s1_merged_new$Sym
s1_merged_new$enrich <- as.numeric(s1_merged_new$enrich)
s1_merged_new$V2 <- as.numeric(s1_merged_new$V2)
s1_merged_new$V3 <- as.numeric(s1_merged_new$V3)
s1_merged_new$V4 <- as.numeric(s1_merged_new$V4)
s1_merged_new$V5 <- as.numeric(s1_merged_new$V5)
s1_merged_new$V6 <- as.numeric(s1_merged_new$V6)
s1_merged_new$V7 <- as.numeric(s1_merged_new$V7)
s1_merged_new$V8 <- as.numeric(s1_merged_new$V8)

res3_1 <- cor.test(as.numeric(s1_merged_new$enrich), as.numeric(s1_merged_new$V2), method = "spearman", exact = FALSE)

s2_enrich_count <- cbind(rownames(s2_enrich_count), s2_enrich_count)
s3_enrich_count <- cbind(rownames(s3_enrich_count), s3_enrich_count)

s1_new <- cbind(rownames(s1), s1)
colnames(s1_new) <- c("Sym", "enrich")
s1_merged = merge(s1_new, results, by.x='Sym', by.y='Sym')
s1_merged_new = merge(s1_new, results(), by.x='Sym', by.y='Sym')

s2_new <- cbind(rownames(s2), s2)
colnames(s2_new) <- c("Sym", "enrich")
s2_merged = merge(s2_new, results, by.x='Sym', by.y='Sym')

s3_new <- cbind(rownames(s3), s3)
colnames(s3_new) <- c("Sym", "enrich")
s3_merged = merge(s3_new, results, by.x='Sym', by.y='Sym')

s4_new <- cbind(rownames(s4), s4)
colnames(s4_new) <- c("Sym", "enrich")
s4_merged = merge(s4_new, results, by.x='Sym', by.y='Sym')

s5_new <- cbind(rownames(s5), s5)
colnames(s5_new) <- c("Sym", "enrich")
s5_merged = merge(s5_new, results, by.x='Sym', by.y='Sym')

s6_new <- cbind(rownames(s6), s6)
colnames(s6_new) <- c("Sym", "enrich")
s6_merged = merge(s6_new, results, by.x='Sym', by.y='Sym')

s7_new <- cbind(rownames(s7), s7)
colnames(s7_new) <- c("Sym", "enrich")
s7_merged = merge(s7_new, results, by.x='Sym', by.y='Sym')

s8_new <- cbind(rownames(s8), s8)
colnames(s8_new) <- c("Sym", "enrich")
s8_merged = merge(s8_new, results, by.x='Sym', by.y='Sym')

s9_new <- cbind(rownames(s9), s9)
colnames(s9_new) <- c("Sym", "enrich")
s9_merged = merge(s9_new, results, by.x='Sym', by.y='Sym')

#find out wheather enrichment data is normally distributed or not
shapiro.test(s1_merged$enrich) # shows error
class(s1_merged$enrich) # it shows factor
shapiro.test(as.numeric(s1_merged$enrich))
#W = 0.96479, p-value < 2.2e-16
#null hypothesis is rejected and hence the data is not normally distributed
rownames(s1_merged) <- s1_merged$Sym
s1_merged$enrich <- as.numeric(s1_merged$enrich)

p1 <- ggscatter(s1_merged, x = "enrich", y = "WT_O", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "spearman",
                xlab = "enrich", ylab = "WT_O")

p2 <- ggscatter(s1_merged, x = "enrich", y = "NAT10_O", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "spearman",
                xlab = "enrich", ylab = "NAT10_O")

p3 <- ggscatter(s1_merged, x = "enrich", y = "WT_i", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "spearman",
                xlab = "enrich", ylab = "WT_i")

p4 <- ggscatter(s1_merged, x = "enrich", y = "NAT10_i", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "spearman",
                xlab = "enrich", ylab = "NAT10_i")

ggarrange(p1, p2, p3, p4, ncol=2, nrow =2)

#col_new=cbind(s1,s2,s3,s4,s5,s6,s7,s8,s9)
#colnames(col)=c("ac4c0_O1","ac4c120_O1","ac4c120_O2","ac4c120_O3","m6a0_O1","m6a0_O2","m6a120_O1","m6a120_O2","m6a120_O3")

shapiro.test(as.numeric(s2_merged$enrich))
#W = 0.96479, p-value < 2.2e-16
#null hypothesis is rejected and hence the data is not normally distributed
rownames(s2_merged) <- s2_merged$Sym
s2_merged$enrich <- as.numeric(s2_merged$enrich)

p1_s2 <- ggscatter(s2_merged, x = "enrich", y = "WT_O", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "spearman",
                xlab = "enrich", ylab = "WT_O")

p2_s2 <- ggscatter(s2_merged, x = "enrich", y = "NAT10_O", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "spearman",
                xlab = "enrich", ylab = "NAT10_O")

p3_s2 <- ggscatter(s2_merged, x = "enrich", y = "WT_i", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "spearman",
                xlab = "enrich", ylab = "WT_i")

p4_s2 <- ggscatter(s2_merged, x = "enrich", y = "NAT10_i", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "spearman",
                xlab = "enrich", ylab = "NAT10_i")

ggarrange(p1_s2, p2_s2, p3_s2, p4_s2, ncol=2, nrow =2)


shapiro.test(as.numeric(s3_merged$enrich))
#W = 0.96479, p-value < 2.2e-16
#null hypothesis is rejected and hence the data is not normally distributed
rownames(s3_merged) <- s3_merged$Sym
s3_merged$enrich <- as.numeric(s3_merged$enrich)

p1_s3 <- ggscatter(s3_merged, x = "enrich", y = "WT_O", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "spearman",
                   xlab = "enrich", ylab = "WT_O")

p2_s3 <- ggscatter(s3_merged, x = "enrich", y = "NAT10_O", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "spearman",
                   xlab = "enrich", ylab = "NAT10_O")

p3_s3 <- ggscatter(s3_merged, x = "enrich", y = "WT_i", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "spearman",
                   xlab = "enrich", ylab = "WT_i")

p4_s3 <- ggscatter(s3_merged, x = "enrich", y = "NAT10_i", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "spearman",
                   xlab = "enrich", ylab = "NAT10_i")

ggarrange(p1_s3, p2_s3, p3_s3, p4_s3, ncol=2, nrow =2)

#Summary
#NOvsNi: correlation analysis seems to be correlating very well,
#EnrichmentvsExpression: for ac4c and NAt10 samples none of the combination seems to be correlating 

###############################
#m6a correlation analysis
m1=as.matrix(read.table("counts_bedtools_KI_c5_m6a.cnt_2172_2184_2201_corel_su_new_new.bed", header = FALSE, row.names = 1))
m1 <- as.data.frame(m1)
colnames(m1) <- c("V2", "V3")

#qqplot to see whether the data is normally distributed or not
# if the data is normally distributed then perform parametric test (pearson correlation) else perform non-parametric test (spearman correlation) 

#does the data follow normal distribution?
shapiro.test(m1$V2)
shapiro.test(m1$V3)

par(mfrow=c(2,2))
#g1 <- ggqqplot(results$WT_O, ylab = "WT_O")

ggscatter(m1, x = "V2", y = "V3", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "spearman",
                   xlab = "V2", ylab = "V3")
############################################
#NOP58 correlation with the enrichment count
#s1=as.matrix(read.table("ac4c0_120_NOP58A_B_tab_annot_new.txt", header = FALSE, row.names = 1))

s1=as.matrix(read.table("ac4c120_NOP58A_B_tab_annot_tab_new.txt", header = FALSE, row.names = 1))
#s1=as.matrix(read.table("ac4c0_NOP58A_B_tab.txt", header = FALSE, row.names = 1))
s1_new <- as.data.frame(s1)

#s1_new$V2 <- as.numeric(s1_new$V2)
s1_new$V3 <- as.numeric(s1_new$V3)
s1_new$V4 <- as.numeric(s1_new$V4)
s1_new$V5 <- as.numeric(s1_new$V5)
s1_new$V6 <- as.numeric(s1_new$V6)
#s1_new$V7 <- as.numeric(s1_new$V7)
#s1_new$V8 <- as.numeric(s1_new$V8)
#s1_new$V9 <- as.numeric(s1_new$V9)

#does the data follow normal distribution?

shapiro.test(s1_new$V3)
shapiro.test(s1_new$V4)
shapiro.test(s1_new$V5)
shapiro.test(s1_new$V6)

p1_s1 <- ggqqplot(s1_new$V3, ylab = "NOP58A")
p2_s1 <- ggqqplot(s1_new$V4, ylab = "NOP58B")
p3_s1 <- ggqqplot(s1_new$V6, ylab = "ac4c_Enrich120")
#p4_s1 <-ggqqplot(s1_new$V9, ylab = "ac4c_Enrich120")

ggarrange(p1_s1, p2_s1, p3_s1, ncol=2, nrow =2)

par(mfrow=c(2,2))
#g1 <- ggqqplot(results$WT_O, ylab = "WT_O")

p1_s1 <- ggscatter(s1_new, x = "V3", y = "V6", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "NOP58A", ylab = "m6a_enrich120")

p2_s1 <- ggscatter(s1_new, x = "V4", y = "V6", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "NOP58B", ylab = "m6a_enrich120")

p3_s1 <-ggscatter(s1_new, x = "V3", y = "V6", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "spearman",
                   xlab = "NOP58A", ylab = "m6a_enrich120")
p4_s1 <-ggscatter(s1_new, x = "V4", y = "V6", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                  xlab = "NOP58B", ylab = "m6a_enrich120")

p5_s1 <-ggscatter(s1_new, x = "V3", y = "V6", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "kendall",
                   xlab = "NOP58A", ylab = "m6a_enrich120")
p6_s1 <-ggscatter(s1_new, x = "V4", y = "V6", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "kendall",
                  xlab = "NOP58B", ylab = "m6a_enrich120")

res2 <- cor.test(s1_new$V3, s1_new$V6,  method="pearson")
res3 <- cor.test(s1_new$V4, s1_new$V6,  method="pearson")

res4 <- cor.test(s1_new$V3, s1_new$V6,  method="spearman")
res5 <- cor.test(s1_new$V4, s1_new$V6,  method="spearman")

res6 <- cor.test(s1_new$V3, s1_new$V6,  method="kendall")
res7 <- cor.test(s1_new$V4, s1_new$V6,  method="kendall")

res2
res3
res4
res5
res6
res7

#p3_s1 <-ggscatter(s1_new, x = "V7", y = "V9", 
 #         add = "reg.line", conf.int = TRUE, 
  #        cor.coef = TRUE, cor.method = "pearson",
   #       xlab = "NOP58A", ylab = "ac4c_enrich0")
#p4_s1 <-ggscatter(s1_new, x = "V8", y = "V9", 
    #      add = "reg.line", conf.int = TRUE, 
     #     cor.coef = TRUE, cor.method = "kendall",
      #    xlab = "NOP58B", ylab = "ac4c_enrich120")
ggarrange(p1_s1, p2_s1, ncol=1, nrow =2)
ggarrange(p1_s1, p2_s1, p3_s1, p4_s1, p5_s1, p6_s1, ncol=2, nrow =3)

#NOP58 correaltion with ac4c (entire set)
#NOP58 correaltion m6a (entire set)
#NOP58 correaltion with ac4c and m6a (here only those sets of m6a has been considered which is found to be overlapping with ac4c)

#DE with different technique without considering ERCC
#Read all of the htseq files woth the help of as.matrix
A1=as.matrix(read.table("RNA_6Aligned.out.sorted.htseq.txt", header = FALSE, row.names = 1))
A2=as.matrix(read.table("RNA_7Aligned.out.sorted.htseq.txt", header = FALSE, row.names = 1))
A3=as.matrix(read.table("RNA_9Aligned.out.sorted.htseq.txt", header = FALSE, row.names = 1))
A4=as.matrix(read.table("RNA_10Aligned.out.sorted.htseq.txt", header = FALSE, row.names = 1))
A5=as.matrix(read.table("101Aligned.out.sorted.htseq.txt", header = FALSE, row.names = 1))
A6=as.matrix(read.table("121Aligned.out.sorted.htseq.txt", header = FALSE, row.names = 1))
# Combine all of the individula matrices into one matric with the help of cbind 
col=cbind(A1,A2,A3,A4,A5,A6)
#Provide column name to all of the columns 
colnames(col)=c(1.1,2.1,3.1,4.1,5.1,6.1)
#Group the samples on the basis of technical replicates and leave biological replicates apart  
group <- factor(c(rep("K1Exp", 1),rep.int("K1Stat", 1),rep.int("DG44Exp", 1),rep.int("DG44Stat", 1),rep.int("SExp", 1),rep.int("SStat", 1)))
# Tells about the dimention of the matrix 
dim(col) #[1] 26635    12
# Displays top 10 genes over all the samples 
head(col) #can not display columns as it is huge 
colnames(col) # Displays the sample size 
#Filtering the reads. Filtering can be performed in mant different ways. Here, we have filtered out dataset on the basis of the genes found atleast in more than two samples(this is basically decided by the number of the biological replicates i.e. if all of the sample sets are having more than two biological replicates then the same gene is supposed to be present in at least two samples) and by CPM>0.5(which can be decided with the help of the size of the library or number of replicated found to be present across all samples).  
#Here in my case all of the CHOS samples are having biological replicates of 4 and CHOK1, DG44 samples are not having any biological replicate which makes the analysis a bit challenging because by considering the least number of the biological samples across all of my samples i.e. 1 would result in overlooking the one which are present in more amount.
myCPM <- cpm(col) # cpm function of the edgeR counts the library in millions 
head(myCPM) # this gives a logical matix in the form  of 'TRUE' and 'FALSE'
########## Until here things are same for both htseq and edgeR

keep <- rowSums(myCPM) >= 1 # We would ike to keep the genes that are atleast one true in each samples minimum nuners of the biological replicates across all the samples is one. 
head(keep)
counts.keep <- col[keep,] # Subset the rows of countdata to keep the more highly expressed genes
dim(counts.keep)  
#[1] 15513     6
summary(keep) # This command will summerize the gene sets passed filter of 0.5. And total of 16032 genes were having CPM count more than .5 as it was found to be true and 10603 genes failed to pass this filter.  
y <- DGEList(counts.keep) # DGE list object in edgeR is used to store count data
y
y$samples$lib.size
par(mfrow=c(1,1))
par(mar=c(7,7,4,2)+0.1,mgp=c(5,1,0)) # read about the description of this command
barplot(y$samples$lib.size,names=colnames(y),las=2) #barplot of the library size 
abline(h=2.00e+07,col="red") # library size of all of the samples were found to be more than 2.75e+07
title(main = "Barplot of Number of reads", xlab = "Sample number", ylab = "Number of reads")
dev.off()
logcounts <- cpm(y,log=TRUE) # log transformation of the library size has been performed so as to shrink the data so as make it easier for visualization 
boxplot(logcounts, xlab="Sample number ", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue") # and the blue line is the median of all the library sizes 
title("Boxplots of logCPMs (unnormalised)")
plotMDS(y) 
title("MDS plot") # Multi dimentional scaling was performed on the unnormalized log transformes cpm's
sampleinfo <- read.delim("All_cell_line_SampleInfo.txt")
sampleinfo
levels(sampleinfo$CellType)
col.cell <- c("purple","orange","dark green","blue","red","black")[sampleinfo$CellType]
data.frame(sampleinfo$CellType,col.cell)
par(mfrow=c(1,1))
plotMDS(y,col=col.cell)
legend("topright",fill=c("purple","orange","dark green"),legend=levels(sampleinfo$CellType),cex=0.5)
title("Cell type")
levels(sampleinfo$Phase)
col.phase <- c("blue","red","dark green","black")[sampleinfo$Phase]
col.phase
plotMDS(y,col=col.phase)
legend("topright",fill=c("blue","red","dark green","black"),legend=levels(sampleinfo$Phase),cex=0.5)
title("Phase")
#hierarchical clustering with heatmaps. As an alternative to the MDSplot for examining the relationship between the samples

var_genes <- apply(logcounts, 1, var) # Variance of each row has been measured with the help of logcount. And the variation was calculated for each row for logcount matrix 
head(var_genes)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:50] # Top 50 most variable gene sets were selected. 
head(select_var)
highly_variable_lcpm <- logcounts[select_var,] # logcount matix were subsetted into parts
dim(highly_variable_lcpm)
head(highly_variable_lcpm)
mypalette <- brewer.pal(11,"RdYlBu") # selection of the coloue was performed with the help of this step
morecols <- colorRampPalette(mypalette) #setting up the colour for vector 
col.cell <- c("purple","orange")[sampleinfo$CellType]
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 50 most variable genes across samples",ColSideColors=col.cell,scale="row") # heatmap function calculates the distant matrix of eucledian distance from logCPM
ynorm <- calcNormFactors(y) # Normalization was performed with the help of calcNormFactors function in edgeR. TMM normalization was performed to eliminated composition bias across libraries. Where the product of the normalization factor and library size defines the effective library size.
#calcNormFactors calculates the normalization across libraries. So the normalization factor more than one means the library size will be scaled down and normalization factor less than one the library size will be scaled up.  
ynorm$samples 
ynorm$samples$norm.factors
#plotMD(logcounts,column = 1, title(main = "MD plot of sample 1 (unnormalised)")) # MD plot of the sample having least normalization factor 
#abline(h=0,col="grey")
#plotMD(logcounts,column = 67)
#abline(h=0,col="grey")
#par(mfrow=c(1,2))
#plotMD(ynorm,column = 1)
#abline(h=0,col="grey")
#plotMD(ynorm,column = 67)
#abline(h=0,col="grey")
group
design <- model.matrix(~ 0 + group)
design
colnames(design) <- levels(group)
design
par(mfrow=c(1,1))
v <- voom(ynorm,design,plot = FALSE) #This plot can also tell us if there are any genes that look really variable in our data, and if we’ve filtered the low counts adequately.
### with plot TRUE: Error in plot.window(...) : need finite 'ylim' values
### with plot FALSE: Error in approxfun(l, rule = 2) : 
#need at least two non-NA values to interpolate
#voom: transforms the readcounts into logCPM while taking into consideration the mean variance of the data
v
names(v)
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(v$E),col="blue")
fit <- lmFit(v) # lmFit function was used to find out the differentially expresssed genes. It uses voom transformed function along with the design matix that we have already specified in the command.
#fit <- glmQLFit(v)
#lmFit calculates group's mean based on the design matirix as well as gene wise variance.
names(fit)

cont.matrix1 <- makeContrasts(K1Exp.K1Stat = K1Exp-K1Stat, DG44Exp.DG44Stat = DG44Exp-DG44Stat, SExp.SStat = SExp-SStat, levels=design)
cont.matrix1 # with the help of this kind of design matix I am trying to compare different types of cell lines 
fit.cont1 <- contrasts.fit(fit, cont.matrix1) #Then we fit the contrast matrix with the help of fit function 
fit.cont1 <- eBayes(fit.cont1)
dim(fit.cont1)

summa.fit1 <- decideTests(fit.cont1)
summary(summa.fit1)

par(mfrow=c(1,3))
#plotMD(fit.cont1,coef=1,status=summa.fit1[,"CS13L.CS13H"])
volcanoplot(fit.cont1,coef=1,highlight=100,names=rownames(fit.cont1))
#volcanoplot(fit.cont1,coef=2,highlight=100,names=rownames(fit.cont1))
#volcanoplot(fit.cont1,coef=3,highlight=100,names=rownames(fit.cont1))
#Volcano of the case1 (volcano.case1)

par(mfrow=c(1,2))
#vennDiagram(vennDia2)
vennDiagram(summa.fit1, include=c("up"), counts.col=c("red"), circle.col = c("red", "blue", "green3"), main = "Rep1_up")
#vennDiagram(vennDia3)
vennDiagram(summa.fit1, include=c("down"), counts.col=c("blue"), circle.col = c("red", "blue", "green3"), main = "Rep1_down")

write.csv(summa.fit1, "summa.fit1.csv")
write.csv(topTable(fit.cont1,coef="Cont1",sort.by="p", number = "Inf"), "rep1Cont1.csv")
write.csv(topTable(fit.cont1,coef="Cont2",sort.by="p", number = "Inf"), "rep1Cont2.csv")
write.csv(topTable(fit.cont1,coef="Cont3",sort.by="p", number = "Inf"), "rep1Cont3.csv")

toptable1_Cont1 <- topTable(fit.cont1,coef="K1Exp.K1Stat",sort.by="p", number = "Inf")
toptable1_Cont2 <- topTable(fit.cont1,coef="DG44Exp.DG44Stat",sort.by="p", number = "Inf")
toptable1_Cont3 <- topTable(fit.cont1,coef="SExp.SStat",sort.by="p", number = "Inf")
nrow(toptable1_Cont3)
head(toptable1_Cont1)