############################
#move to bamfiles path  #####
#move to bamfiles path  #####
#############################
####bam pathway which will be told by xuelan or xiangyu ######
pathway=./QL/sortedByCoord.out.bam
cd $pathway

############################
#work in R environment  #####
#work in R environment  #####
#############################

project <- c("NFF")##课题负责人
deal_design <- c("15b","control")##课题的case比上control
significant_cutoff <- c(1)##显著水平的cutoff，默认两倍，所填数字是log2之后的
organism <- "mouse"##物种为人就是hsa，是老鼠就是mouse


##################################################################
##################################################################
#######################DO NOT MODIFY ANYTHING#####################
#######################DO NOT MODIFY ANYTHING#####################
#######################DO NOT MODIFY ANYTHING#####################
##################################################################

file_path <- paste(getwd(),"/work_file",sep="")
sample_sampletable.path <- getwd()
dir.create(file_path)
load(file="./GRCm38_anno/ebg_GRCm38.RData")
load(file="./GRCm38_anno/txdb_GRCm38.RData")
deal_case <- deal_design[1]
dela_control <- deal_design[2]
deal_names <- paste(deal_case,dela_control,sep="_VS_")
project_pathway <- paste(".",project,sep="/")
se.save_RData <- paste(project_pathway,"0_RNAseq_se.RData",sep="_")
cluster_sample_heatmap_rld <- paste(project_pathway,"1_cluster_similarity_rld_pheatmap.pdf",sep="_")
cluster_sample_heatmap_vsd <- paste(project_pathway,"1_cluster_similarity_vsd_pheatmap.pdf",sep="_")
cluster_sample_pca <- paste(project_pathway,"1_cluster_similarity_pca.pdf",sep="_")
tpm_csv <- paste(project_pathway,"2_tpm.csv",sep="_")
res_1_file.csv <- paste(project_pathway,"3",deal_names,"result.csv",sep="_")
res_1_filter <- c("deal",deal_design)
colnames_1 <- deal_names
res_1_file_sym_name <- paste(project_pathway,"4",deal_names,"symbol_and_anno.csv",sep="_")
ALL.CSV_FILE <- paste(project_pathway,"5",deal_names,"count_tpm_symbol_and_anno.csv",sep="_")
KEGGupres_1_file.csv <- paste(project_pathway,"6",deal_names,"KEGG_UP.csv",sep="_")
KEGGdownres_1_file.csv <- paste(project_pathway,"7",deal_names,"KEGG_DOWN.csv",sep="_")
KEGGupres_1_pdf <- paste(project_pathway,"8",deal_names,"KEGG_UP.pdf",sep="_")
KEGGdownres_1_pdf <- paste(project_pathway,"9",deal_names,"KEGG_DOWN.pdf",sep="_")
KEGGres_1_all_pdf <- paste(project_pathway,"10",deal_names,"KEGG_all.pdf",sep="_")
KEGGres_1_all_title <- paste(deal_names,"KEGG Pathway Enrichment",sep=" ")
GOres_1_all_UP_csv <- paste(project_pathway,"11",deal_names,"GO_UP.csv",sep="_")
GOres_1_all_DOWN_csv <- paste(project_pathway,"12",deal_names,"GO_DOWN.csv",sep="_")
GOres_1_all_UP_pdf <- paste(project_pathway,"13",deal_names,"GO_UP.pdf",sep="_")
GOres_1_all_DOWN_pdf <- paste(project_pathway,"14",deal_names,"GO_DOWN.pdf",sep="_")

suppressPackageStartupMessages({
	library(Rsamtools)
	library(GenomicFeatures)
	library(GenomicAlignments)
	library(BiocParallel)
	library(pheatmap)
	library(RColorBrewer)
	library(PoiClaClu)
	library(org.Mm.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(pathview)
	library(org.Hs.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(ggplot2)
})
indir <- file.path(sample_sampletable.path)
csvfile <- file.path(indir, "sample_sampletable.csv")
sampleTable <- read.csv(csvfile, row.names = 1)
sampleTable
filenames <- file.path(indir, sampleTable$ids)
file.exists(filenames)
bamfiles <- BamFileList(filenames, yieldSize=2000000)
seqinfo(bamfiles[1])
register(SerialParam())
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
mode="Union",
singleEnd=TRUE,
ignore.strand=TRUE,
fragments=FALSE)
setwd(file_path)
save(se,file=se.save_RData)

colData(se) <- DataFrame(sampleTable)
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ deal)
countdata <- assay(se)
coldata <- colData(se)

rld <- rlog(dds, blind = FALSE)
sampleDists <- dist(t(assay(rld)))
samplePoisDistMatrix <- as.matrix(sampleDists)
rownames(samplePoisDistMatrix) <- paste(rld$sample)
colnames(samplePoisDistMatrix) <- NULL
pdf(file=cluster_sample_heatmap_rld)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix,
clustering_distance_rows = sampleDists,
clustering_distance_cols = sampleDists,
col = colors)
dev.off()

vsd <- vst(dds, blind = FALSE)
sampleDists <- dist(t(assay(vsd)))
samplePoisDistMatrix <- as.matrix(sampleDists)
rownames(samplePoisDistMatrix) <- paste(vsd$sample)
colnames(samplePoisDistMatrix) <- NULL
pdf(file=cluster_sample_heatmap_vsd)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix,
clustering_distance_rows = sampleDists,
clustering_distance_cols = sampleDists,
col = colors)
dev.off()

pdf(file=cluster_sample_pca)
library("ggplot2")
pcaData <- plotPCA(rld, intgroup = c("sample","order"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample)) +
geom_point(size =1) +
geom_text(aes(label=order))+
xlab(paste0("PC1: ", percentVar[1], "% variance")) +
ylab(paste0("PC2: ", percentVar[2], "% variance")) +
coord_fixed()
dev.off()

exons<-as.data.frame(ebg)
exon_length<-data.frame(gene_name=c(0),length=c(0))
for(i in 1:nrow(assay(se))){
  exonsfromgene<-exons[which(exons$group==i),]
  length<-sum(exonsfromgene$width)
  gene_name<-exonsfromgene[1,2]
  exon_length<-rbind(exon_length,c(gene_name,length))
}
exon_length<-exon_length[-1,]
#求出每个样本的总read数（百万）
whole_length<-as.data.frame(round(colSums(assay(dds))/1e6,1))
#计算FPKM矩阵
countM<-as.data.frame(assay(se))
FPKM<-countM
for(i in 1:nrow(FPKM)){
  gene<-rownames(FPKM[i,])
  length<-as.numeric(exon_length[which(exon_length$gene_name==gene),]$length)
  for(j in 1:ncol(FPKM)){
    FPKM[i,j]<-FPKM[i,j]/(whole_length[j,]*length*0.001) 
  }
}
fpkmToTpm <- function(FPKM){
    exp(log(FPKM) - log(sum(FPKM)) + log(1e6))
}
tpm <- fpkmToTpm(FPKM)
count <- assay(dds)
colnames(tpm) <- sampleTable$sample
colnames(count) <- sampleTable$sample
count_and_tpm <- cbind(count,tpm)
write.csv(tpm,file=tpm_csv)
dds <- DESeq(dds)
res_1 <- results(dds, contrast=res_1_filter)
colnames(res_1) = paste(colnames_1,colnames(res_1),sep="_")
write.csv(res_1, file = res_1_file.csv)

if (organism %in% c("mouse")) {
	anno_data=org.Mm.eg.db
	print("organism is mouse")
} else {
	anno_data=org.Hs.eg.db
	print("organism is human")
}

res_1$symbol <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="ENSEMBL",
						column ="SYMBOL",
						multiVals="first")
res_1$entrez <- mapIds(x = anno_data,
                        keys = rownames(res_1),
						keytype ="ENSEMBL",
						column ="ENTREZID",
						multiVals="first")
AA <- res_1$symbol
AA <- as.character(AA)
res_1$GENENAME <- mapIds(x = anno_data,
                        keys = AA,
						keytype ="SYMBOL",
						column ="GENENAME",
						multiVals="first")
write.csv(res_1, file = res_1_file_sym_name)
all_summry <- cbind(count_and_tpm,res_1)
write.csv(all_summry, file = ALL.CSV_FILE)

res_1 <- na.omit(res_1)
y <- res_1[c(1:nrow(res_1)),6]
res_1 <- res_1[with(res_1,y<0.05),]
y <- res_1[c(1:nrow(res_1)),2]
upres_1 <- res_1[with(res_1,y>=significant_cutoff),]
downres_1  <- res_1[with(res_1,y<=-significant_cutoff),]
ee	<-as.matrix(upres_1$entrez)
	dd <- as.vector(ee)
KEGGupres_1 <- enrichKEGG(gene =dd, 
					organism = organism, 
					keyType = "ncbi-geneid",
					 pvalueCutoff = 0.05,
				       pAdjustMethod = "BH", 
				       minGSSize = 10, 
				       maxGSSize = 500,
				       qvalueCutoff = 0.2, 
				       use_internal_data = FALSE)
ee	<-as.matrix(downres_1$entrez)
	dd <- as.vector(ee)
KEGGdownres_1 <- enrichKEGG(gene =dd, 
					organism = organism, 
					keyType = "ncbi-geneid",
					 pvalueCutoff = 0.05,
				       pAdjustMethod = "BH", 
				       minGSSize = 10, 
				       maxGSSize = 500,
				       qvalueCutoff = 0.2, 
				       use_internal_data = FALSE)

write.csv(KEGGupres_1, file = KEGGupres_1_file.csv)
write.csv(KEGGdownres_1, file = KEGGdownres_1_file.csv)
pdf(file = KEGGupres_1_pdf)
dotplot(KEGGupres_1, showCategory=20)
dev.off()
pdf(file = KEGGdownres_1_pdf)
dotplot(KEGGdownres_1, showCategory=20)
dev.off()

pdf(file = KEGGres_1_all_pdf)
KEGGupres_1uporder <- KEGGupres_1[order((KEGGupres_1$Count),decreasing = TRUE),]
KEGGdownres_1downorder <- KEGGdownres_1[order((KEGGdownres_1$Count),decreasing = TRUE),]
down_kegg<-KEGGdownres_1downorder[KEGGdownres_1downorder$p.adjust<0.05,];
down_kegg$group=-1
up_kegg<-KEGGupres_1uporder[KEGGupres_1uporder$p.adjust<0.05,];
up_kegg$group=1
down_kegg1 <- down_kegg[-(21:nrow(down_kegg)),]
up_kegg1 <- up_kegg[-(21:nrow(up_kegg)),]
  dat=rbind(up_kegg1,down_kegg1)
  colnames(dat)
    dat$p.adjust = -log10(dat$p.adjust)
  dat$p.adjust=dat$p.adjust*dat$group 
  dat=dat[order(dat$p.adjust,decreasing = F),]
  library("ggplot2")
  g_kegg<- ggplot(dat, aes(x=reorder(Description,order(p.adjust, decreasing = F)), y=p.adjust, fill=group)) + 
    geom_bar(stat="identity") + 
    scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
    scale_x_discrete(name ="Pathway names") +
    scale_y_continuous(name ="-log10p_adjust") +
    coord_flip() +
    ggtitle(KEGGres_1_all_title)
  print(g_kegg)
dev.off()

ee	<-as.matrix(upres_1$entrez)
	dd <- as.vector(ee)
GOupres_1_all <- enrichGO(gene = dd, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOupres_1_all,GOres_1_all_UP_csv)

ff <- barplot(GOupres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file =GOres_1_all_UP_pdf)
ff1
dev.off()

ee	<-as.matrix(downres_1$entrez)
	dd <- as.vector(ee)
GOdownres_1_all <- enrichGO(gene = dd, 
	           OrgDb = anno_data,
				ont = "all", 
		             pvalueCutoff = 0.01, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05,
                     minGSSize = 10, 
                     maxGSSize = 500, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOdownres_1_all,file=GOres_1_all_DOWN_csv)
ff <- barplot(GOdownres_1_all, split="ONTOLOGY",showCategory=15) + facet_grid(ONTOLOGY~., scale="free")
ff1 <- ff + theme(axis.text.y = element_text(size = 8)) +labs(title = NULL)
pdf(file =GOres_1_all_DOWN_pdf)
ff1
dev.off()



