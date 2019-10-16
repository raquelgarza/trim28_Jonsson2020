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
# Functions ----
more_10 <- function(row){
  return(length(which(row > 10)))
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
# EMX animals ----
emx <- fread('6_TEtranscripts/invivo_bd/invivo_bd.cntTable', data.table = F)
colnames(emx)[-1] <- paste(unlist(lapply(strsplit(colnames(emx)[-1], '/'), `[[`, 5)), unlist(lapply(strsplit(colnames(emx)[-1], '/'), `[[`, 4)), sep='_')
# Gene count
emx_gene <- subset(emx, startsWith(emx$`gene/TE`, 'ENSMUS'))
colnames(emx_gene)[1] <- 'gene_id'
rownames(emx_gene) <- emx_gene$gene_id
emx_gene <- emx_gene[,-1]
# TE count
emx_TE <- subset(emx, !startsWith(emx$`gene/TE`, 'ENSMUS'))
colnames(emx_TE)[1] <- 'TE_id'
rownames(emx_TE) <- emx_TE$TE_id
emx_TE <- emx_TE[,-1]

emx_ctx <- colnames(emx_gene)[grepl('ctx', colnames(emx_gene))]
emx_coldata <- data.frame(samples=emx_ctx, condition=unlist(lapply(strsplit(emx_ctx, '_'), `[[`, 2)))
rownames(emx_coldata) <- emx_coldata$samples
emx_coldata$samples <- as.character(emx_coldata$samples)

# Gene DEA
emx_genes_dds <- DESeqDataSetFromMatrix(emx_gene[,rownames(emx_coldata)], emx_coldata, design = ~ condition)
emx_genes_dds <- DESeq(emx_genes_dds)
emx_genes_res <- results(emx_genes_dds)
emx_genes_exp <- getAverage(emx_genes_dds)

p_gene_meanplot_emx <- meanPlot_cus(emx_genes_exp$Mean, test=emx_genes_res, p=0.05, c1='ko', c2='ctrl',ttl='Gene DEA in Cre Loxp experiment', repel = FALSE) + labs(title="", subtitle="")
ggsave(p_gene_meanplot_emx, file="6_TEtranscripts/invivo_bd/plots/gene_meanplot.svg", width=20, height=20, units="cm", dpi=96)


emx_gene_norm <- counts(emx_genes_dds, normalized = TRUE)
emx_gene_norm <- merge(emx_gene_norm, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
emx_gene_norm_trim28 <- as.data.frame(t(emx_gene_norm[which(emx_gene_norm$gene_name == 'Trim28'),as.character(emx_coldata$samples)]))
colnames(emx_gene_norm_trim28) <- 'value'
emx_gene_norm_trim28$sample <- rownames(emx_gene_norm_trim28)
emx_gene_norm_trim28$Condition <- ifelse(emx_coldata[as.character(emx_gene_norm_trim28$sample),'condition'] == 'ko', 'Knock out', 'Control')

p_trim28_emx <- ggplot(emx_gene_norm_trim28, aes(x=Condition, y=value, fill=Condition)) + geom_boxplot()+ theme_classic() + labs(y="Normalized read count", title='EMX animals - Trim28 expression')+
  theme(
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15,margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 15,margin = margin(t = 0, r = 20, b = 0, l = 0))) + scale_y_continuous(limits = c(0, 450))

ggsave(p_trim28_emx, file='6_TEtranscripts/invivo_bd/plots/trim28.svg', width=20, height=20, units="cm", dpi=96)

emx_upreg <- as.data.frame(as.matrix(subset(emx_genes_res, emx_genes_res$log2FoldChange > 0 & emx_genes_res$padj < 0.05)))
emx_upreg <- merge(emx_upreg, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
colnames(emx_upreg) <- c('Gene ID', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'P-adj', 'Gene Name')

write.xlsx(emx_upreg[,c(8,1,3,7)], '6_TEtranscripts/invivo_bd/upregulated_genes.xlsx', row.names = FALSE, col.names=TRUE)
write.table(unique(emx_upreg$`Gene Name`), 'GO_analysis/invivo_bd/upregulated/sign_upreg_genes.txt', row.names = FALSE, col.names=FALSE, quote = F)

emx_dwreg <- as.data.frame(as.matrix(subset(emx_genes_res, emx_genes_res$log2FoldChange < 0 & emx_genes_res$padj < 0.05)))
emx_dwreg <- merge(emx_dwreg, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
colnames(emx_dwreg) <- c('Gene ID', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'P-adj', 'Gene Name')

write.table(unique(emx_dwreg$`Gene Name`), 'GO_analysis/invivo_bd/downregulated/sign_downreg_genes.txt', row.names = FALSE, col.names=FALSE, quote=F)

emx_gene_expressed <- emx_gene[which(apply(emx_gene, 1, more_10) > 0),]
write(rownames(emx_gene_expressed), 'GO_analysis/invivo_bd/background_emx.txt')

# TE DEA
emx_TE_dds <- DESeqDataSetFromMatrix(emx_TE[,rownames(emx_coldata)], emx_coldata, design = ~ condition)
emx_TE_dds <- DESeq(emx_TE_dds)
emx_TE_res <- results(emx_TE_dds)
emx_TE_exp <- getAverage(emx_TE_dds)

p_TE_meanplot_emx <- meanPlot_cus(emx_TE_exp$Mean, test=emx_TE_res, p=0.05, c1='ko', c2='ctrl',ttl='TE DEA in Cre Loxp experiment', repel = FALSE, col3='firebrick', col2='black') + labs(title="", subtitle="")
ggsave(p_TE_meanplot_emx, file="6_TEtranscripts/invivo_bd/plots/TE_meanplot.svg", width=20, height=20, units="cm", dpi=96)

emx_TE_norm <- emx_TE
emx_TE_norm[] <- mapply('/', emx_TE_norm, emx_genes_dds$sizeFactor)
emx_TE_norm$TE_subfamily <- as.character(unlist(lapply(strsplit(rownames(emx_TE_norm), ':'), `[[`, 1)))
emx_TE_norm$TE_family <- unlist(lapply(strsplit(rownames(emx_TE_norm), ':'), `[[`, 2))
emx_TE_norm$TE_class <- unlist(lapply(strsplit(rownames(emx_TE_norm), ':'), `[[`, 3))
emx_TE_norm$TE_family <- as.character(ifelse(endsWith((emx_TE_norm$TE_family), '?'), substr(emx_TE_norm$TE_family, 1, (nchar(emx_TE_norm$TE_family)-1)),emx_TE_norm$TE_family))
emx_TE_norm$TE_class <- as.character(ifelse(endsWith((emx_TE_norm$TE_class), '?'), substr(emx_TE_norm$TE_class, 1, (nchar(emx_TE_norm$TE_class)-1)),emx_TE_norm$TE_class))
emx_TE_norm <- emx_TE_norm[which(emx_TE_norm$TE_class %in% c('LINE', 'LTR', 'SINE')),]

emx_TE_signdiff_condition <- melt(emx_TE_exp$Mean)
emx_TE_signdiff_condition <- emx_TE_signdiff_condition[emx_TE_signdiff_condition$Var1 %in% rownames(subset(emx_TE_res, emx_TE_res$padj < 0.05)),]
emx_TE_signdiff_condition$value <- log2(emx_TE_signdiff_condition$value+0.5)
colnames(emx_TE_signdiff_condition) <- c('TE_subfamily', 'Condition', 'log2Mean')

p_signdiff_TE_emx <- ggplot(data=emx_TE_signdiff_condition, aes(x=TE_subfamily, y=log2Mean)) +   
  geom_bar(aes(fill = Condition, width=0.7), position = "dodge", stat="identity") + theme_classic()+labs(y="log2(mean normalized read counts)", x="TE subfamily", fill="Condition")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ scale_fill_manual(values=c("steelblue", "tomato2"))

ggsave(p_signdiff_TE_emx, file='6_TEtranscripts/invivo_bd/plots/signdiff_TE.svg', width=25, height=10, units="cm", dpi=96)


# NPC ----
npc <- fread('/Volumes/Seagate Backup /trim28/09.10.19/6_TEtranscripts/invitro_crispr/invitro_crispr.cntTable', data.table = F)
colnames(npc)[-1] <- paste(unlist(lapply(strsplit(colnames(npc)[-1], '/'), `[[`, 5)), unlist(lapply(strsplit(colnames(npc)[-1], '/'), `[[`, 4)), unlist(lapply(strsplit(colnames(npc)[-1], '/'), `[[`, 6)), sep='_')
npc_gene <- subset(npc, startsWith(npc$`gene/TE`, 'ENSMUS'))
colnames(npc_gene)[1] <- 'gene_id'
rownames(npc_gene) <- npc_gene$gene_id
npc_gene <- npc_gene[,-1]

npc_TE <- subset(npc, !startsWith(npc$`gene/TE`, 'ENSMUS'))
colnames(npc_TE)[1] <- 'TE_id'
rownames(npc_TE) <- npc_TE$TE_id
npc_TE <- npc_TE[,-1]

npc_coldata <- data.frame(samples=colnames(npc_gene), condition=unlist(lapply(strsplit(colnames(npc_gene), '_'), `[[`, 2)))
rownames(npc_coldata) <- npc_coldata$samples
npc_coldata$samples <- as.character(npc_coldata$samples)

npc_genes_dds <- DESeqDataSetFromMatrix(npc_gene[,rownames(npc_coldata)], npc_coldata, design = ~ condition)
npc_genes_dds <- DESeq(npc_genes_dds)
npc_genes_res <- results(npc_genes_dds)
npc_genes_exp <- getAverage(npc_genes_dds)

p_gene_meanplot_npc <- meanPlot_cus(npc_genes_exp$Mean, test=npc_genes_res, p=0.05, c1='ko', c2='ctrl',ttl='', repel = FALSE) + labs(title="", subtitle="")
ggsave(p_gene_meanplot_npc, file="6_TEtranscripts/invitro_crispr/plots/gene_meanplot.svg", width=20, height=20, units="cm", dpi=96)

npc_TE_dds <- DESeqDataSetFromMatrix(npc_TE[,rownames(npc_coldata)], npc_coldata, design = ~ condition)
npc_TE_dds <- DESeq(npc_TE_dds)
npc_TE_res <- results(npc_TE_dds)
npc_TE_exp <- getAverage(npc_TE_dds)

p_TE_meanplot_npc <- meanPlot_cus(npc_TE_exp$Mean, test=npc_TE_res, col2 = 'black', col3 = 'firebrick', p=0.05, c1='ko', c2='ctrl',ttl='', repel = F)  + labs(title="", subtitle="")
ggsave(p_TE_meanplot_npc, file="6_TEtranscripts/invitro_crispr/plots/TE_meanplot.svg", width=20, height=20, units="cm", dpi=96)

npc_gene_norm <- counts(npc_genes_dds, normalized = TRUE)
npc_gene_norm <- merge(npc_gene_norm, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
npc_gene_norm_trim28 <- as.data.frame(t(npc_gene_norm[which(npc_gene_norm$gene_name == 'Trim28'),as.character(npc_coldata$samples)]))
colnames(npc_gene_norm_trim28) <- 'value'
npc_gene_norm_trim28$sample <- rownames(npc_gene_norm_trim28)
npc_gene_norm_trim28$Condition <- ifelse(npc_coldata[as.character(npc_gene_norm_trim28$sample),'condition'] == 'ko', 'Knock out', 'Control')
npc_gene_norm_trim28$gRNA <- unlist(lapply(strsplit(npc_gene_norm_trim28$sample, '_'), `[[` , 1))

p_trim28_npc <- ggplot(npc_gene_norm_trim28, aes(x=Condition, y=value, fill=gRNA)) + geom_boxplot() + theme_classic() + labs(y="Median-of-ratios normalized read count", title='TRIM28 in vitro CRISPR knock out')+
  theme(
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15,margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 15,margin = margin(t = 0, r = 20, b = 0, l = 0))) + scale_y_continuous(limits = c(0, 3500))


ggsave(p_trim28_npc, file="6_TEtranscripts/invitro_crispr/plots/trim28.svg", width=20, height=20, units="cm", dpi=96)

npc_upreg <- as.data.frame(as.matrix(subset(npc_genes_res, npc_genes_res$log2FoldChange > 0 & npc_genes_res$padj < 0.05)))
npc_upreg <- merge(npc_upreg, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
colnames(npc_upreg) <- c('Gene ID', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'P-adj', 'Gene Name')

write.xlsx(npc_upreg[,c(8,1,3,7)], '6_TEtranscripts/invitro_crispr/upregulated_genes.xlsx', row.names = FALSE, col.names=TRUE)
write.table(unique(npc_upreg$`Gene Name`), 'GO_analysis/invitro_crispr/upregulated/sign_upreg_genes.txt', row.names = FALSE, col.names=FALSE, quote = F)

npc_dwnreg <- as.data.frame(as.matrix(subset(npc_genes_res, npc_genes_res$log2FoldChange < 0 & npc_genes_res$padj < 0.05)))
npc_dwnreg <- merge(npc_dwnreg, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
colnames(npc_dwnreg) <- c('Gene ID', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'P-adj', 'Gene Name')

write.xlsx(npc_dwnreg[,c(8,1,3,7)], 'multimapping/6_TEtranscripts/invitro_crispr/downregulated_genes.xlsx', row.names = FALSE, col.names=TRUE)
write.table(unique(npc_dwnreg$`Gene Name`), 'GO_analysis/invitro_crispr/downregulated/sign_downreg_genes.txt', row.names = FALSE, col.names=FALSE, quote = F)

npc_gene_expressed <- npc_gene[which(apply(npc_gene, 1, more_10) > 0),]
write(rownames(npc_gene_expressed), 'GO_analysis/invitro_crispr/background_npc.txt')

npc_TE_norm <- npc_TE
npc_TE_norm[] <- mapply('/', npc_TE_norm, npc_genes_dds$sizeFactor)
npc_TE_norm$TE_subfamily <- as.character(unlist(lapply(strsplit(rownames(npc_TE_norm), ':'), `[[`, 1)))
npc_TE_norm$TE_family <- unlist(lapply(strsplit(rownames(npc_TE_norm), ':'), `[[`, 2))
npc_TE_norm$TE_class <- unlist(lapply(strsplit(rownames(npc_TE_norm), ':'), `[[`, 3))
npc_TE_norm$TE_family <- as.character(ifelse(endsWith((npc_TE_norm$TE_family), '?'), substr(npc_TE_norm$TE_family, 1, (nchar(npc_TE_norm$TE_family)-1)),npc_TE_norm$TE_family))
npc_TE_norm$TE_class <- as.character(ifelse(endsWith((npc_TE_norm$TE_class), '?'), substr(npc_TE_norm$TE_class, 1, (nchar(npc_TE_norm$TE_class)-1)),npc_TE_norm$TE_class))
npc_TE_norm <- subset(npc_TE_norm, npc_TE_norm$TE_class %in% c('LINE', 'LTR', 'SINE'))

npc_TE_norm$ID <- rownames(npc_TE_norm)
npc_TE_signdiff_retro <- rownames(npc_TE_res[which(npc_TE_res$padj < 0.05 ),])[which(rownames(npc_TE_res[which(npc_TE_res$padj < 0.05 ),]) %in% rownames(npc_TE_norm))]
npc_TE_signdiff_retro <- npc_TE_norm[npc_TE_signdiff_retro, as.character(npc_coldata$samples)]

npc_TE_signdiff_condition <- melt(npc_TE_exp$Mean)
npc_TE_signdiff_condition <- npc_TE_signdiff_condition[npc_TE_signdiff_condition$Var1 %in% rownames(subset(npc_TE_res, npc_TE_res$padj < 0.05)),]
npc_TE_signdiff_condition$value <- log2(npc_TE_signdiff_condition$value+0.5)
colnames(npc_TE_signdiff_condition) <- c('TE_subfamily', 'Condition', 'log2Mean')

p_signdiff_TE_npc <- ggplot(data=npc_TE_signdiff_condition, aes(x=TE_subfamily, y=log2Mean)) +   
  geom_bar(aes(fill = Condition, width=0.7), position = "dodge", stat="identity") + theme_classic()+labs(y="log2(mean normalized read counts)", x="TE subfamily", fill="Condition")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ scale_fill_manual(values=c("steelblue", "tomato2"))

ggsave(p_signdiff_TE_npc, file='6_TEtranscripts/invitro_crispr/plots/signdiff_TE.svg', width=25, height=15, units="cm", dpi=96)


# CRISPR invivo ---- 
invivo_crispr <-  fread('6_TEtranscripts/invivo_crispr/invivo_crispr.cntTable', data.table = F)
colnames(invivo_crispr)[-1] <- paste(unlist(lapply(strsplit(colnames(invivo_crispr)[-1], '/'), `[[`, 5)), unlist(lapply(strsplit(colnames(invivo_crispr)[-1], '/'), `[[`, 4)), sep='_')
invivo_crispr_gene <- subset(invivo_crispr, startsWith(invivo_crispr$`gene/TE`, 'ENSMUS'))
colnames(invivo_crispr_gene)[1] <- 'gene_id'
rownames(invivo_crispr_gene) <- invivo_crispr_gene$gene_id
invivo_crispr_gene <- invivo_crispr_gene[,-1]

invivo_crispr_coldata <- data.frame(samples=colnames(invivo_crispr_gene), condition=unlist(lapply(strsplit(colnames(invivo_crispr_gene), '_'), `[[`, 3)))
invivo_crispr_coldata$type <- unlist(lapply(strsplit(colnames(invivo_crispr_gene), '_'), `[[`, 2))
rownames(invivo_crispr_coldata) <- invivo_crispr_coldata$samples
invivo_crispr_coldata$samples <- as.character(invivo_crispr_coldata$samples)

invivo_crispr_genes_dds <- DESeqDataSetFromMatrix(invivo_crispr_gene[,rownames(invivo_crispr_coldata)], invivo_crispr_coldata, design = ~ condition)
invivo_crispr_genes_dds <- DESeq(invivo_crispr_genes_dds)
invivo_crispr_genes_res <- results(invivo_crispr_genes_dds)
invivo_crispr_genes_exp <- getAverage(invivo_crispr_genes_dds)

p_gene_meanplot_invivocrispr <- meanPlot_cus(invivo_crispr_genes_exp$Mean, test=invivo_crispr_genes_res, p=0.05, c1='ko', c2='ctrl',ttl='Gene DEA in invivo CRISPR KO experiment', repel = FALSE) + labs(title="", subtitle="")
ggsave(p_gene_meanplot_invivocrispr, file="6_TEtranscripts/invivo_crispr/plots/gene_meanplot.svg", width=20, height=20, units="cm", dpi=96)

invivo_crispr_gene_norm <- counts(invivo_crispr_genes_dds, normalized = TRUE)
invivo_crispr_gene_norm <- merge(invivo_crispr_gene_norm, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
invivo_crispr_gene_norm_trim28 <- as.data.frame(t(invivo_crispr_gene_norm[which(invivo_crispr_gene_norm$gene_name == 'Trim28'),as.character(invivo_crispr_coldata$samples)]))
colnames(invivo_crispr_gene_norm_trim28) <- 'value'
invivo_crispr_gene_norm_trim28$sample <- rownames(invivo_crispr_gene_norm_trim28)
invivo_crispr_gene_norm_trim28$Condition <- ifelse(invivo_crispr_coldata[as.character(invivo_crispr_gene_norm_trim28$sample),'condition'] == 'ko', 'Knock out', 'Control')
invivo_crispr_gene_norm_trim28$Condition <- ifelse(endsWith(invivo_crispr_gene_norm_trim28$sample, 'ctrl'), 'Control', ifelse(startsWith(invivo_crispr_gene_norm_trim28$sample, 'flexcas'), 'FlexCas', 'Cas'))
invivo_crispr_gene_norm_trim28$Condition <- factor(invivo_crispr_gene_norm_trim28$Condition, levels = c('Control', 'Cas', 'FlexCas'))
invivo_crispr_gene_norm_trim28$gRNA <- ifelse(endsWith(invivo_crispr_gene_norm_trim28$sample, 'ctrl'), 'Control', ifelse(startsWith(invivo_crispr_gene_norm_trim28$sample, 'flexcas'), paste('FlexCas', unlist(lapply(str_split(invivo_crispr_gene_norm_trim28$sample, '_'), `[[`, 2))), paste('Cas', unlist(lapply(str_split(invivo_crispr_gene_norm_trim28$sample, '_'), `[[`, 2)))))
invivo_crispr_gene_norm_trim28$gRNA <- factor(invivo_crispr_gene_norm_trim28$gRNA, levels=c('Control', 'Cas g3', 'Cas g4', 'Cas g13', 'Cas g3g4g13', 'FlexCas g3', 'FlexCas g4', 'FlexCas g13', 'FlexCas g3g4g13'))

p_trim28_invivocrispr <- ggplot(invivo_crispr_gene_norm_trim28[rownames(invivo_crispr_gene_norm_trim28),], aes(x=Condition, y=value, fill=Condition)) + geom_boxplot(alpha=0.5)+ theme_classic() + labs(y="Median-of-ratios normalized read count", title='Trim28 invivo CRISPR knock out')+
  geom_point(aes(colour=gRNA), size=5) +
  guides(fill=FALSE, col = guide_legend(nrow = 5))+
  scale_colour_manual(values=c('Cas g13'="#3c913a", 'Cas g3'="#7fe87d", 'Cas g4'="#51c44f", 'Cas g3g4g13'="#255c24", 
                               Control="#eb4444", 
                               'FlexCas g13'="#418796", 'FlexCas g3'="#6ae0fc", 'FlexCas g4'="#56b2c7", 'FlexCas g3g4g13'="#234952")) +
  labs(colour='Samples') +
  
  theme(
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15,margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 15,margin = margin(t = 0, r = 20, b = 0, l = 0)))  + scale_y_continuous(limits = c(0, 1600))
ggsave(p_trim28_invivocrispr, file="6_TEtranscripts/invivo_crispr/plots/trim28.svg", width=20, height=20, units="cm", dpi=96)

invivo_crispr_upreg <- as.data.frame(as.matrix(subset(invivo_crispr_genes_res, invivo_crispr_genes_res$log2FoldChange > 0 & invivo_crispr_genes_res$padj < 0.05)))
invivo_crispr_upreg <- merge(invivo_crispr_upreg, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
colnames(invivo_crispr_upreg) <- c('Gene ID', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'P-adj', 'Gene Name')

write.xlsx(invivo_crispr_upreg[,c(8,1,3,7)], '6_TEtranscripts/invivo_crispr/upregulated_genes.xlsx', row.names = FALSE, col.names=TRUE)

invivo_crispr_TE <- subset(invivo_crispr, !startsWith(invivo_crispr$`gene/TE`, 'ENSMUS'))
colnames(invivo_crispr_TE)[1] <- 'TE_id'
rownames(invivo_crispr_TE) <- invivo_crispr_TE$TE_id
invivo_crispr_TE <- invivo_crispr_TE[,-1]

invivo_crispr_TE$TE_id <- as.character(unlist(lapply(strsplit(rownames(invivo_crispr_TE), ':'), `[[`, 1)))
invivo_crispr_TE$TE_family <- unlist(lapply(strsplit(rownames(invivo_crispr_TE), ':'), `[[`, 2))
invivo_crispr_TE$TE_class <- unlist(lapply(strsplit(rownames(invivo_crispr_TE), ':'), `[[`, 3))
invivo_crispr_TE$TE_family <- as.character(ifelse(endsWith((invivo_crispr_TE$TE_family), '?'), substr(invivo_crispr_TE$TE_family, 1, (nchar(invivo_crispr_TE$TE_family)-1)),invivo_crispr_TE$TE_family))
invivo_crispr_TE$TE_class <- as.character(ifelse(endsWith((invivo_crispr_TE$TE_class), '?'), substr(invivo_crispr_TE$TE_class, 1, (nchar(invivo_crispr_TE$TE_class)-1)),invivo_crispr_TE$TE_class))
invivo_crispr_TE <- subset(invivo_crispr_TE, invivo_crispr_TE$TE_class %in% c('LINE', 'LTR', 'SINE'))

invivo_crispr_TE_dds <- DESeqDataSetFromMatrix(invivo_crispr_TE[,names(invivo_crispr_genes_dds$sizeFactor)], invivo_crispr_coldata, design = ~ condition)
invivo_crispr_TE_dds <- DESeq(invivo_crispr_TE_dds)
invivo_crispr_TE_res <- results(invivo_crispr_TE_dds)
invivo_crispr_TE_exp <- getAverage(invivo_crispr_TE_dds)

p_TE_meanplot_invivocrispr <- meanPlot_cus(invivo_crispr_TE_exp$Mean, test=invivo_crispr_TE_res, p=0.05, c1='ko', col2='black', col3='firebrick', c2='ctrl',ttl='TE subfamilies at the CRISPR KO experiment', repel = FALSE) + labs(title="", subtitle="")
ggsave(p_TE_meanplot_invivocrispr, file="6_TEtranscripts/invivo_crispr/plots/TE_meanplot.svg", width=20, height=20, units="cm", dpi=96)


# Adult invivo ----
invivo_adult <- fread('6_TEtranscripts/invivo_adult/invivo_adult.cntTable', data.table = F)
colnames(invivo_adult)[-1] <- paste(unlist(lapply(strsplit(colnames(invivo_adult)[-1], '/'), `[[`, 5)), unlist(lapply(strsplit(colnames(invivo_adult)[-1], '/'), `[[`, 4)), sep='_')
colnames(invivo_adult)[-1] <- paste(unlist(lapply(strsplit(colnames(invivo_adult)[-1], '_'), `[[`, 1)), ifelse(unlist(lapply(strsplit(colnames(invivo_adult)[-1], '_'), `[[`, 2)) == 'ctrl', 'ko', 'ctrl'), sep='_')
invivo_adult_gene <- subset(invivo_adult, startsWith(invivo_adult$`gene/TE`, 'ENSMUS'))
colnames(invivo_adult_gene)[1] <- 'gene_id'
rownames(invivo_adult_gene) <- invivo_adult_gene$gene_id
invivo_adult_gene <- invivo_adult_gene[,-1]

invivo_adult_coldata <- data.frame(samples=colnames(invivo_adult_gene), condition=unlist(lapply(strsplit(colnames(invivo_adult_gene), '_'), `[[`, 2)))
rownames(invivo_adult_coldata) <- invivo_adult_coldata$samples
invivo_adult_coldata$samples <- as.character(invivo_adult_coldata$samples)

invivo_adult_genes_dds <- DESeqDataSetFromMatrix(invivo_adult_gene[,rownames(invivo_adult_coldata)], invivo_adult_coldata, design = ~ condition)
invivo_adult_genes_dds <- DESeq(invivo_adult_genes_dds)
invivo_adult_genes_res <- results(invivo_adult_genes_dds)
invivo_adult_genes_exp <- getAverage(invivo_adult_genes_dds)

p_gene_meanplot_invivoadult <- meanPlot_cus(invivo_adult_genes_exp$Mean, test=invivo_adult_genes_res, p=0.05, c1='ko', c2='ctrl',ttl='Gene DEA at floxed adult', repel = FALSE) + labs(title="", subtitle="")
ggsave(p_gene_meanplot_invivocrispr, file="6_TEtranscripts/invivo_adult/plots/gene_meanplot.svg", width=20, height=20, units="cm", dpi=96)

invivo_adult_gene_norm <- counts(invivo_adult_genes_dds, normalized = TRUE)
invivo_adult_gene_norm <- merge(invivo_adult_gene_norm, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
invivo_adult_gene_norm_trim28 <- as.data.frame(t(invivo_adult_gene_norm[which(invivo_adult_gene_norm$gene_name == 'Trim28'),as.character(invivo_adult_coldata$samples)]))
colnames(invivo_adult_gene_norm_trim28) <- 'value'
invivo_adult_gene_norm_trim28$sample <- rownames(invivo_adult_gene_norm_trim28)
invivo_adult_gene_norm_trim28$Condition <- ifelse(invivo_adult_coldata[as.character(invivo_adult_gene_norm_trim28$sample),'condition'] == 'ko', 'Knock out', 'Control')

p_trim28_invivoadult <- ggplot(invivo_adult_gene_norm_trim28, aes(x=Condition, y=value, fill=Condition)) + geom_boxplot()+ theme_classic() + labs(y="Median-of-ratios normalized read count", title='TRIM28 adult knock out')+
  theme(
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15,margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 15,margin = margin(t = 0, r = 20, b = 0, l = 0)))+ scale_y_continuous(limits = c(0, 3000))


ggsave(p_trim28_invivoadult, file="6_TEtranscripts/invivo_adult/plots/trim28.svg", width=20, height=20, units="cm", dpi=96)

invivo_adult_upreg <- as.data.frame(as.matrix(subset(invivo_adult_genes_res, invivo_adult_genes_res$log2FoldChange > 0 & invivo_adult_genes_res$padj < 0.05)))
invivo_adult_upreg <- merge(invivo_adult_upreg, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
colnames(invivo_adult_upreg) <- c('Gene ID', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'P-adj', 'Gene Name')

write.xlsx(invivo_adult_upreg[,c(8,1,3,7)], '6_TEtranscripts/invivo_adult/upregulated_genes.xlsx', row.names = FALSE, col.names=TRUE)

invivo_adult_dwnreg <- as.data.frame(as.matrix(subset(invivo_adult_genes_res, invivo_adult_genes_res$log2FoldChange < 0 & invivo_adult_genes_res$padj < 0.05)))
invivo_adult_dwnreg <- merge(invivo_adult_dwnreg, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
colnames(invivo_adult_dwnreg) <- c('Gene ID', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'P-adj', 'Gene Name')

write.xlsx(invivo_adult_dwnreg[,c(8,1,3,7)], '6_TEtranscripts/invivo_adult/downregulated_genes.xlsx', row.names = FALSE, col.names=TRUE)

invivo_adult_TEs <- subset(invivo_adult, !startsWith(invivo_adult$`gene/TE`, 'ENSMUS'))
colnames(invivo_adult_TEs)[1] <- 'TE_id'
rownames(invivo_adult_TEs) <- invivo_adult_TEs$TE_id
invivo_adult_TEs <- invivo_adult_TEs[,-1]

invivo_adult_TEs$TE_id <- as.character(unlist(lapply(strsplit(rownames(invivo_adult_TEs), ':'), `[[`, 1)))
invivo_adult_TEs$TE_family <- unlist(lapply(strsplit(rownames(invivo_adult_TEs), ':'), `[[`, 2))
invivo_adult_TEs$TE_class <- unlist(lapply(strsplit(rownames(invivo_adult_TEs), ':'), `[[`, 3))
invivo_adult_TEs$TE_family <- as.character(ifelse(endsWith((invivo_adult_TEs$TE_family), '?'), substr(invivo_adult_TEs$TE_family, 1, (nchar(invivo_adult_TEs$TE_family)-1)),invivo_adult_TEs$TE_family))
invivo_adult_TEs$TE_class <- as.character(ifelse(endsWith((invivo_adult_TEs$TE_class), '?'), substr(invivo_adult_TEs$TE_class, 1, (nchar(invivo_adult_TEs$TE_class)-1)),invivo_adult_TEs$TE_class))
invivo_adult_TEs <- subset(invivo_adult_TEs, invivo_adult_TEs$TE_class %in% c('LINE', 'LTR', 'SINE'))

invivo_adult_TEs_dds <- DESeqDataSetFromMatrix(invivo_adult_TEs[,rownames(invivo_adult_coldata)], invivo_adult_coldata, design = ~ condition)
invivo_adult_TEs_dds <- DESeq(invivo_adult_TEs_dds)
invivo_adult_TEs_res <- results(invivo_adult_TEs_dds)
invivo_adult_TEs_vst <- varianceStabilizingTransformation(invivo_adult_TEs_dds)
invivo_adult_TE_prcomp <- make_pca(t(assay(invivo_adult_TEs_vst))[ , apply(t(assay(invivo_adult_TEs_vst)), 2, var) != 0],invivo_adult_coldata,c1='Control',c2='Ko', folder='multimapping/plots/invivo_adult/TE_', TRUE, '\nAdult Trim28 (TEs only)')
invivo_adult_TEs_exp <- getAverage(invivo_adult_TEs_dds)

p_TE_meanplot_invivoadult <- meanPlot_cus(invivo_adult_TEs_exp$Mean, test=invivo_adult_TEs_res, p=0.05, c1='ko', c2='ctrl', col2='black', col3='firebrick', ttl='', repel = FALSE) + labs(title="", subtitle="")
ggsave(p_TE_meanplot_invivoadult, file="6_TEtranscripts/invivo_adult/plots/TE_meanplot.svg", width=20, height=20, units="cm", dpi=96)

# EMX GO analysis ----
# Invivo bd without adult 
write.xlsx(emx_upreg[which(!emx_upreg$`Gene Name` %in% invivo_adult_upreg$`Gene Name`),c(8,1,3,7)], 'GO_analysis/invivo_bd/upregulated/upregulated_genes_emx_not_invivo_adult.xlsx', col.names = T, row.names = F)
write(unique(emx_upreg[which(!emx_upreg$`Gene Name` %in% invivo_adult_upreg$`Gene Name`),'Gene Name']), 'GO_analysis/invivo_bd/upregulated/sign_upreg_genes_not_in_adult.txt')

gencode <- fread('/Volumes/Seagate Backup /annotation/mouse/gencode/gencode.vM20.annotation.bed', data.table = F)
colnames(gencode) <- c('chr', 'start', 'end', 'gene_id', 'strand')
rownames(gencode) <- gencode$gene_id
write.table(gencode[emx_upreg$`Gene ID`,], '6_TEtranscripts/invivo_bd/upregulated_genes_emx_not_invivo_adult.bed', col.names = T, row.names = F, sep='\t', quote = F)

write(unique(emx_dwreg[which(!emx_dwreg$`Gene Name` %in% invivo_adult_dwnreg$`Gene Name`),'Gene Name']), 'GO_analysis/invivo_bd/downregulated/sign_downreg_genes_not_in_adult.txt')

# Invivo bd viral defence 
viral_defence <- c("IFI16","IFI27","MS2","OAS1","IRF7","OASL","OAS2","OAS3","ISG20","MX1","IFIH1","IFIT3","IFI6","STAT1","IFIT2","ISG15","DDX58","DHX58","IFITM2","IFI35","B2M","IRF9","IFITM1","IFIT1","MX2")
viral_defence <- tools::toTitleCase(tolower(viral_defence))
viral_defence <- data.frame(gene_name=viral_defence)

viral_defence <- unique(merge(viral_defence, transcript_gene[,c(2,3)], by='gene_name'))

emx_gene_norm_viral_defence <- emx_gene_norm[which(emx_gene_norm$gene_name %in% viral_defence$gene_name),c('gene_name', as.character(emx_coldata$samples))]
rownames(emx_gene_norm_viral_defence) <- emx_gene_norm_viral_defence$gene_name
emx_gene_norm_viral_defence_more10 <- emx_gene_norm_viral_defence[which(rowSums(emx_gene_norm_viral_defence[,emx_coldata$samples]) > 10),]
emx_gene_norm_viral_defence_more10 <- merge(emx_gene_norm_viral_defence_more10, unique(transcript_gene[,c(2,3)]), by='gene_name')

emx_viral_defence_fc <- emx_genes_res_df[viral_defence$gene_id,c('log2FoldChange', 'ci_low', 'ci_high'), drop=FALSE]
emx_viral_defence_fc <- merge(emx_viral_defence_fc, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
emx_viral_defence_fc$reg <- ifelse(emx_viral_defence_fc$log2FoldChange < 0, 'logFC < 0', 'logFC > 0')
emx_viral_defence_fc$gene_name <- factor(emx_viral_defence_fc$gene_name, levels=emx_viral_defence_fc[order(emx_viral_defence_fc$log2FoldChange, decreasing=TRUE),'gene_name'])

emx_viral_defence_fc <- emx_viral_defence_fc[which(emx_viral_defence_fc$Row.names %in% emx_gene_norm_viral_defence_more10$gene_id),]

p_viral_defence <- ggplot() + 
  geom_errorbar(data=emx_viral_defence_fc, mapping=aes(x=gene_name, ymin=ci_low, ymax=ci_high, colour=reg), width=0.2) + 
  geom_point(stat='identity',data=emx_viral_defence_fc, mapping=aes(x=gene_name, y=log2FoldChange, colour=reg)) + coord_flip()+theme_bw() + scale_colour_manual(values=c("steelblue", 'tomato')) +
  guides(fill=FALSE)+
  labs(x="Gene name", y= 'log2 Fold Change', title='Change of expression on viral defence genes at CreLoxP Trim28 KO\nError lines of 95% confidence intervals', colour="KO Relative expression")+
  theme(
    plot.title = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15,margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 15,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  geom_hline(yintercept = 0, color = "black", size=0.5)+
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -0.5, ymax = 0.5), fill = 'lightpink',  alpha = 0.3,size = 2)
ggsave(p_viral_defence, file="6_TEtranscripts/invivo_bd/plots/viral_defence.svg", width=20, height=20, units="cm", dpi=96)

# EMX ERVK upregulated ----
ERVK_count <- fread('/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/ERVK_count_matrix_2.csv', data.table = F)
ERVK_coldata <- data.frame(sample=unlist(lapply(str_split(colnames(ERVK_count)[7:ncol(ERVK_count)], '/'), `[[`, 5)),
                      condition=unlist(lapply(str_split(colnames(ERVK_count)[7:ncol(ERVK_count)], '/'), `[[`, 4)))
colnames(ERVK_count)[7:ncol(ERVK_count)] <- unlist(lapply(str_split(colnames(ERVK_count)[7:ncol(ERVK_count)], '/'), `[[`, 5))

rownames(ERVK_coldata) <- ERVK_coldata$sample
rownames(ERVK_count) <- ERVK_count$Geneid

ERVK_dds <- DESeqDataSetFromMatrix(ERVK_count[,rownames(ERVK_coldata)], ERVK_coldata, design = ~ condition)
ERVK_dds <- DESeq(ERVK_dds)
ERVK_res <- results(ERVK_dds)

ERVK_upreg <- rownames(ERVK_res[which(ERVK_res$log2FoldChange > 0),])
ERVK_upreg <- ERVK_upreg[which(ERVK_upreg %in% rownames(ERVK_count))]
write.table(ERVK_count[ERVK_upreg,c('Chr', 'Start', 'End', 'Geneid', 'Strand')], col.names = F, row.names = F, quote = F, sep='\t', file='10_nearbygenes/upregulated/emx_upreg_ERVK.bed')

ERVK_not_upreg <- rownames(ERVK_res[which(!ERVK_res$log2FoldChange > 0 | is.na(ERVK_res$log2FoldChange)),])
ERVK_not_upreg <- ERVK_not_upreg[which(ERVK_not_upreg %in% rownames(ERVK_res))]
write.table(ERVK_count[ERVK_not_upreg,c('Chr', 'Start', 'End', 'Geneid', 'Strand')], col.names = F, row.names = F, quote = F, sep='\t', file='10_nearbygenes/not_upregulated/emx_notupreg_ERVK.bed')

