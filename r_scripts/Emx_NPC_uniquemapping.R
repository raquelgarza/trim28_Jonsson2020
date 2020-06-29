# Settings ----
setwd('/Volumes/Seagate Backup /trim28/09.10.19/')
library('data.table')
library('stringr')
library("RColorBrewer")
library("xlsx")
library("deseqAbstraction")
library('pheatmap')
library('ggplot2')
library('DESeq2')
library('rjson')
library('Hmisc')
# Functions ----
getSignName <- function(x,p,l=0) {
  
  up <- x[!is.na(x$padj) & x$padj < p & x$log2FoldChange > l,]
  down <- x[!is.na(x$padj) & x$padj < p & x$log2FoldChange < -l,]
  return(list(up=rownames(up),down=rownames(down)))
  
}
getAverage <- function(dds) {
  
  baseMeanPerLvl <- sapply( levels(dds$condition), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$condition == lvl] ) )
  baseSDPerLvl <- sapply( levels(dds$condition), function(lvl) apply( counts(dds,normalized=TRUE)[,dds$condition == lvl],1,sd ) )
  colnames(baseSDPerLvl) <- paste("st.dev:",colnames(baseSDPerLvl),sep="")
  return(list(Mean=baseMeanPerLvl,SD=baseSDPerLvl))
  
}
meanPlot_cus <- function(exp,test,c1 = "condition 1",c2 = "condition 2",p=.05,l=0,id=F, ttl="", 
                         repel=TRUE, col1="firebrick3", col2="steelblue4", col3="black", highlights=NA){
  sign <- getSignName(x = test,p = p,l = l)
  u <- sign$up
  d <- sign$down
  
  #color up and down sign..
  colVec <- ifelse(test = (rownames(exp) %in% u),
                   yes = col1,
                   no = ifelse(test = (rownames(exp) %in% d),
                               yes = col2, no =col3))
  colVec[is.na(colVec)] <- "steelblue" ## if NA make sure it's not counted as <p
  #size of points
  cexVec <- ifelse(test = (rownames(exp) %in% u),
                   yes = 0.35,
                   no = ifelse(test = (rownames(exp) %in% d),
                               yes = 0.35, no = 0.3))
  
  exp_log <- as.data.frame(log2(exp[,c(c1, c2)]))
  exp_log$Name <- rownames(exp_log)
  
  exp_log$reg <- factor(ifelse(exp_log$Name %in% u, paste('upregulated in ', c1, ' (', length(u), ')', sep =''),
                               ifelse(exp_log$Name %in% d, paste('downregulated in ', c1,' (', length(d), ')', sep =''), paste('not significant', ' (', (nrow(test) - length(u) - length(d)), ')', sep=''))))
  
  library(ggrepel)
  if(repel == TRUE){
    plt <- ggplot(exp_log, aes(x=get(c2), y=get(c1), label=Name, color=reg)) + geom_point(aes(size=cexVec))+ scale_color_manual(values=c(col2, col3, col1))+ scale_size_continuous(range=c(1,2), guide=FALSE)+ geom_text_repel(data = subset(exp_log, Name %in% u),direction    = "y", nudge_y = 0.4, nudge_x = -0.5)
  }
  else{
    plt <- ggplot(exp_log, aes(x=get(c2), y=get(c1), color=reg)) + geom_point(aes(size=cexVec))+ scale_color_manual(values=c(col2, col3, col1))+ scale_size_continuous(range=c(1,2), guide=FALSE)
  }
  plt <- plt + labs(x=paste("log2(mean ",c2,")",sep=""), 
                    y=paste("log2(mean ",c1,")",sep=""),
                    title=paste(ttl, paste(c1," vs. ",c2,sep=""), sep = ': '),
                    subtitle=paste("p-adj < ",p,", log2(fc) > ",l,sep=""))+theme(
                      plot.title = element_text( size=14, face="bold"),  panel.grid.major = element_line(colour="white"),
                      panel.grid.minor = element_line(colour="white"), panel.background = element_rect(fill = "white",
                                                                                                       colour = "white",
                                                                                                       size = 0.5, linetype = "solid"),
                      axis.line = element_line(size = 0.5, linetype = "solid",
                                               colour = "black"), 
                      legend.title=element_blank()) 
  
  
  if(id==T) {
    
    identify(log2(exp[,1]),log2(exp[,2]),labels = rownames(exp))
    
  }
  
  if(!is.na(highlights)){
    plt <- plt + geom_point(data=exp_log[highlights,], aes(x=get(c2), y=get(c1)), colour="springgreen3", size=5, shape=1, stroke=2)
  }
  return(plt)
  
}

# General ----
TE_classification <- fread('/Volumes/Seagate Backup /annotation/mouse/repeatmasker/mm10_rmsk_TE_classification.tab', data.table = FALSE, header=FALSE, fill=TRUE)
colnames(TE_classification) <- c('TE_id', 'TE_subfamily', 'TE_family', 'TE_class')
transcript_gene <- fread('/Volumes/Seagate Backup /annotation/mouse/gencode/gencode.vM20.annotation.transc.gene.tab', data.table = FALSE, header=FALSE)
colnames(transcript_gene) <- c('transcript_id', 'gene_id', 'gene_name', 'gene_type')

# NPC ----
# Read unique mapping quantification of TEs (repeatmasker)
npc_TE <- fread('/Volumes/Seagate Backup /trim28/09.10.19/1_uniqmapping/2_readcount/TE_count_matrix_2.csv', data.table = F)
npc_TE <- npc_TE[,colnames(npc_TE)[c(1, (which(sapply(strsplit(colnames(npc_TE)[7:ncol(npc_TE)], '/'), `[[`, 3) == 'npc')+6))]]

colnames(npc_TE)[2:ncol(npc_TE)] <- paste(sapply(strsplit(colnames(npc_TE)[2:ncol(npc_TE)], '/'), `[[`, 4), unlist(lapply(strsplit(colnames(npc_TE)[2:ncol(npc_TE)], '/'), `[[`, 5)), sapply(strsplit(colnames(npc_TE)[2:ncol(npc_TE)], '/'), `[[`, 6), sep='_')
colnames(npc_TE)[1] <- 'TE_id'
rownames(npc_TE) <- npc_TE$TE_id
npc_TE <- npc_TE[,-1]

# Create the metadata
npc_coldata <- data.frame(samples=colnames(npc_TE), condition=unlist(lapply(strsplit(colnames(npc_TE), '_'), `[[`, 1)))
rownames(npc_coldata) <- npc_coldata$samples
npc_coldata$samples <- as.character(npc_coldata$samples)

# Differential TE expression analysis testing for differences between conditions (Trim28 KO vs control)
npc_TE_dds <- DESeqDataSetFromMatrix(npc_TE[,rownames(npc_coldata)], npc_coldata, design = ~ condition)
npc_TE_dds <- DESeq(npc_TE_dds)
npc_TE_res <- results(npc_TE_dds)
# Calculate average expression on normalized reads on TEs
npc_TE_exp <- getAverage(npc_TE_dds)
npc_TE_vst <- varianceStabilizingTransformation(npc_TE_dds)

# Plot PCA based on TE expression
npc_TE_pca <- plotPCA(npc_TE_vst) + theme_classic() + ggtitle("PCA - In vitro CRISPR TE expression")
ggsave(npc_TE_pca, file="/Volumes/Seagate Backup /trim28/09.10.19/1_uniqmapping/npc/plots/TE_pca.svg", width=20, height=20, units="cm", dpi=96)

# Mean plot of individual elements
p_TE_meanplot_npc <- meanPlot_cus(npc_TE_exp$Mean, test=npc_TE_res, l=0.5, p=0.05, c1='ko', c2='ctrl',ttl='', repel = F) + labs(title="", subtitle="")
ggsave(p_TE_meanplot_npc, file="/Volumes/Seagate Backup /trim28/09.10.19/1_uniqmapping/npc/plots/TE_meanplot_0.5.svg", width=20, height=20, units="cm", dpi=320)
ggsave(p_TE_meanplot_npc, file="/Volumes/Seagate Backup /trim28/09.10.19/1_uniqmapping/npc/plots/TE_meanplot_0.5_npc.png", dpi=320)

# Read gene expression based on unique mapping 
npc_gene <- fread('/Volumes/Seagate Backup /trim28/09.10.19/1_uniqmapping/2_readcount/gene_count_matrix_2.csv', data.table = F)
npc_gene <- npc_gene[,colnames(npc_gene)[c(1, (which(sapply(strsplit(colnames(npc_gene)[7:ncol(npc_gene)], '/'), `[[`, 3) == 'npc')+6))]]

colnames(npc_gene)[2:ncol(npc_gene)] <- paste(sapply(strsplit(colnames(npc_gene)[2:ncol(npc_gene)], '/'), `[[`, 4), unlist(lapply(strsplit(colnames(npc_gene)[2:ncol(npc_gene)], '/'), `[[`, 5)), sapply(strsplit(colnames(npc_gene)[2:ncol(npc_gene)], '/'), `[[`, 6), sep='_')
colnames(npc_gene)[1] <- 'gene_id'
rownames(npc_gene) <- npc_gene$gene_id
npc_gene <- npc_gene[,-1]

# Differential gene expression analysis testing for differences between conditions (Trim28 KO vs control)
npc_gene_dds <- DESeqDataSetFromMatrix(npc_gene[,rownames(npc_coldata)], npc_coldata, design = ~ condition)
npc_gene_dds <- DESeq(npc_gene_dds)
npc_gene_res <- results(npc_gene_dds)
npc_gene_exp <- getAverage(npc_gene_dds)

# Plot PCA based on protein coding genes
npc_gene_dds_protein <- npc_gene_dds[rownames(npc_gene_dds) %in% unique(subset(transcript_gene, transcript_gene$gene_type == 'protein_coding')$gene_id),]
npc_gene_vst_protein <- varianceStabilizingTransformation(npc_gene_dds_protein)
npc_gene_pca_protein <- plotPCA(npc_gene_vst_protein) + theme_classic() + ggtitle("PCA - In vitro CRISPR gene expression")
ggsave(npc_gene_pca_protein, file="/Volumes/Seagate Backup /trim28/09.10.19/1_uniqmapping/npc/plots/gene_pca_protein.svg", width=20, height=20, units="cm", dpi=96)

# Mean plot of gene expression
p_gene_meanplot_npc <- meanPlot_cus(npc_gene_exp$Mean, test=npc_gene_res, p=0.05, c1='ko', c2='ctrl',ttl='', repel = F, l=0.5)  + labs(title="", subtitle="")
ggsave(p_gene_meanplot_npc, file="/Volumes/Seagate Backup /trim28/09.10.19/1_uniqmapping/npc/plots/gene_meanplot_0.5.svg", width=20, height=20, units="cm", dpi=96)

# Emx ----
# Read unique mapping quantification of TEs (repeatmasker)
emx_TE <- fread('/Volumes/Seagate Backup /trim28/09.10.19/1_uniqmapping/2_readcount/TE_count_matrix_2.csv', data.table = F)
emx_TE <- emx_TE[,colnames(emx_TE)[c(1, (which(sapply(strsplit(colnames(emx_TE)[7:ncol(emx_TE)], '/'), `[[`, 3) == 'emx')+6))]]

colnames(emx_TE)[2:ncol(emx_TE)] <- paste(sapply(strsplit(colnames(emx_TE)[2:ncol(emx_TE)], '/'), `[[`, 4), unlist(lapply(strsplit(colnames(emx_TE)[2:ncol(emx_TE)], '/'), `[[`, 5)), sep='_')
colnames(emx_TE)[1] <- 'TE_id'
rownames(emx_TE) <- emx_TE$TE_id
emx_TE <- emx_TE[,-1]

# Create metadata
emx_coldata <- data.frame(samples=colnames(emx_TE), condition=unlist(lapply(strsplit(colnames(emx_TE), '_'), `[[`, 1)))
rownames(emx_coldata) <- emx_coldata$samples
emx_coldata$samples <- as.character(emx_coldata$samples)
emx_coldata <- subset(emx_coldata, emx_coldata$condition != 'str')

# emx_TE <- merge(emx_TE, TE_classification[,c(1,2)], by.x='row.names', by.y='TE_id')
# emx_TE <- emx_TE[,-1]
# emx_TE <- aggregate(emx_TE[,-ncol(emx_TE)], by=list(emx_TE$TE_subfamily), FUN=sum)
# rownames(emx_TE) <- emx_TE$Group.1
# emx_TE <- emx_TE[,-1]

# Differential TE expression analysis testing for differences between conditions (Trim28 KO vs control)
emx_TE_dds <- DESeqDataSetFromMatrix(emx_TE[,rownames(emx_coldata)], emx_coldata, design = ~ condition)
emx_TE_dds <- DESeq(emx_TE_dds)
emx_TE_res <- results(emx_TE_dds)
emx_TE_exp <- getAverage(emx_TE_dds)
emx_TE_vst <- varianceStabilizingTransformation(emx_TE_dds)

# Plot PCA based on TE expression
emx_TE_pca <- plotPCA(emx_TE_vst) + theme_classic() + ggtitle("PCA - In vitro CRISPR TE expression")
ggsave(emx_TE_pca, file="/Volumes/Seagate Backup /trim28/09.10.19/1_uniqmapping/npc/plots/TE_pca.svg", width=20, height=20, units="cm", dpi=96)

# Mean plot of individual elements 
p_TE_meanplot_emx <- meanPlot_cus(emx_TE_exp$Mean, test=emx_TE_res, l=0.5, p=0.05, c1='ko', c2='ctrl',ttl='', repel = F) + labs(title="", subtitle="")

ggsave(p_TE_meanplot_emx, file="/Volumes/Seagate Backup /trim28/09.10.19/1_uniqmapping/emx/plots/TE_meanplot_0.5.svg", width=20, height=20, units="cm", dpi=320)
ggsave(p_TE_meanplot_emx, file="/Volumes/Seagate Backup /trim28/09.10.19/1_uniqmapping/emx/plots/TE_meanplot_0.5_emx.png", dpi=320)



