#WT vs NAT10
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

filter <- apply(col, 1, function(x) length(x[x>5])>=2)
filter2 <- apply(col, 1, function(x) length(x[x>1])>=2)
filtered <- col[filter,]
filtered2 <- col[filter2,]
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]

x <- as.factor(rep(c("WTOi","NAT10Oi"), each=6))

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

up_new <- deg[deg$logFC > 0,]
down_new <- deg[deg$logFC < 0,]
nrow(up_new)
nrow(down_new)
nrow(lrt)

# while comparing WTi vs WTo, none of the genes are found to be sig up/down regulated 
WtiovsNAT10io <- cbind(rownames(lrt$table),lrt$table)
colnames(WtiovsNAT10io) <- c("SYMBOL","logFC","logCPM","LR","PValue")                 
write.table(WtiovsNAT10io, file = "WT_oivsNAT10_oi.txt", row.names = F, sep = "\t", quote = F)

up <- cbind(rownames(up), up)
rownames(up) <- NULL
colnames(up) <- c("SYMBOL", "logFC", "logCPM", "LR","PValue","FDR")
write.table(up, file = "up_modified.txt", row.names = F, sep = "\t", quote = F)

down <- cbind(rownames(down), down)
rownames(down) <- NULL
colnames(down) <- c("SYMBOL", "logFC", "logCPM", "LR","PValue","FDR")
write.table(down, file = "down_modified.txt", row.names = F, sep = "\t", quote = F)

#annotation of protien coding genes
setwd("~/GenomeDK/faststorage/reference/UCSC/human/GRCh37")
annot = read.delim('gencode.v36lift37.annotation.filtered_modified.txt', header=F)
colnames(annot) <- c("Gene_ID","gene_type","gene_name")

merged.file = merge(up, annot, by.x='SYMBOL', by.y='Gene_ID')
write.table(merged.file, file = "up_annot.txt", row.names = F, sep = "\t", quote = F)

merged.file = merge(down, annot, by.x='SYMBOL', by.y='Gene_ID')
write.table(merged.file, file = "down_annot.txt", row.names = F, sep = "\t", quote = F)

#annotation of non-protien coding genes
#Doing with non coding doesnt make much difference because all of the genes reflected here are already present in the coding 
setwd("~/GenomeDK/faststorage/reference/UCSC/human/GRCh37")
annot = read.delim('gencode.v36lift37.long_noncoding_RNAs_filtered_modified.gtf', header=F)
colnames(annot) <- c("Gene_ID","gene_type","gene_name")

merged.file = merge(up, annot, by.x='SYMBOL', by.y='Gene_ID')
write.table(merged.file, file = "up_annotNC.txt", row.names = F, sep = "\t", quote = F)

merged.file = merge(down, annot, by.x='SYMBOL', by.y='Gene_ID')
write.table(merged.file, file = "down_annotNC.txt", row.names = F, sep = "\t", quote = F)
