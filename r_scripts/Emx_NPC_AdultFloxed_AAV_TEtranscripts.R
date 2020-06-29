# Settings ----
setwd('/Volumes/Seagate Backup /trim28/09.10.19/')
library('data.table')
library('stringr')
library("RColorBrewer")
library("xlsx")
library('pheatmap')
library('ggplot2')
library('DESeq2')
library('rjson')
library('Hmisc')
# Functions ----
## getSignName ##
# Get significantly different gene names. 
# Taken from source code of the package deseqAbstraction which is no longer available on github.
# Credits to Per L. Brattås
# Parameters:
# x = results object from deseq
# p = padj threshold for significance
# l = log2FC threshold for significance
getSignName <- function(x,p,l=0) {
  up <- x[!is.na(x$padj) & x$padj < p & x$log2FoldChange > l,]
  down <- x[!is.na(x$padj) & x$padj < p & x$log2FoldChange < -l,]
  return(list(up=rownames(up),down=rownames(down)))
}
## getAverage ##
# Get average expression (normalized by median of ratios) of each of the conditions in a deseq object.
# Taken from source code of the package deseqAbstraction which is no longer available on github.
# Credits to Per L. Brattås
# Parameters:
# dds = deseq object
getAverage <- function(dds) {
  baseMeanPerLvl <- sapply( levels(dds$condition), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$condition == lvl] ) )
  baseSDPerLvl <- sapply( levels(dds$condition), function(lvl) apply( counts(dds,normalized=TRUE)[,dds$condition == lvl],1,sd ) )
  colnames(baseSDPerLvl) <- paste("st.dev:",colnames(baseSDPerLvl),sep="")
  return(list(Mean=baseMeanPerLvl,SD=baseSDPerLvl))
}
## GO_json ##
# Parses json files from Panther DB into a relatively well sized plots
# It outputs a barplot of log2FC of each of the selected GO terms (those that passes the criteria
# specified by the parameters) combined with a dot plot showing FDRs.
# Parameters:
# file_name = path to the json file
# prefix = output path or prefix of the resulting plots in pdf
# botlog2FC = filter GO terms which have more than this log2FC
# toplog2FC = filter GO terms which have less than this log2FC
# morethan = filter GO terms which have less than this number of genes supporting it
# ttl = title of the plot
GO_json <- function(file_name, prefix, botlog2fc=NULL, toplog2fc=NULL, morethan=NULL, ttl=''){
  json <- fromJSON(file=file_name)
  json$overrepresentation$group[[1]]$result[[1]]$input_list$number_in_list
  results <- vector("list",  length(json$overrepresentation$group))
  for(j in 1:length(results)){
    
    for(i in 1:length(json$overrepresentation$group[[j]]$result)){
      
      if(length(json$overrepresentation$group[[j]]$result[[i]]) > 1){
        if(!is.null(json$overrepresentation$group[[j]]$result[[i]]$input_list$number_in_list)){
          if(json$overrepresentation$group[[j]]$result[[i]]$input_list$number_in_list > 0 ){
            tmp <- data.frame(GO_label= json$overrepresentation$group[[j]]$result[[i]]$term$label,
                              GO_id = json$overrepresentation$group[[j]]$result[[i]]$term$id,
                              gene_id=paste(json$overrepresentation$group[[j]]$result[[i]]$input_list$mapped_id_list$mapped_id, collapse=', '),
                              fold_enrichment=json$overrepresentation$group[[j]]$result[[i]]$input_list$fold_enrichment,
                              sign=json$overrepresentation$group[[j]]$result[[i]]$input_list$plus_minus,
                              expected=json$overrepresentation$group[[j]]$result[[i]]$input_list$expected,
                              fdr=json$overrepresentation$group[[j]]$result[[i]]$input_list$fdr,
                              pval = json$overrepresentation$group[[j]]$result[[i]]$input_list$pValue)
            
            
            results[[j]] <- rbind(results[[j]], tmp)
          }
        }
      }
    }
    
  }
  
  results <- results[-which(sapply(results, is.null))]
  results <- results[which(sapply(results, nrow) != 0)]
  results_df <- do.call(rbind, results)
  
  results_df$log2fold_enrichment <- log2(results_df$fold_enrichment+0.5)
  
  results_df$gene_id <- as.character(results_df$gene_id)
  
  if(startsWith(results_df$gene_id[1], 'ENS')){
    transcript_gene <- fread('annotation/mm10/gencode/gencode.vM20.annotation.transc.gene.tab', data.table = FALSE, header=FALSE)
    colnames(transcript_gene) <- c('transcript_id', 'gene_id', 'gene_name', 'gene_type')
    
    transcript_gene$gene_id_panther <- unlist(lapply(str_split(transcript_gene$gene_id, '[.]'), `[[` , 1))
    results_df$gene_name <- rep(NULL, nrow(results_df))
    
    for (i in 1:length(results_df$gene_id)){
      tmp <- as.data.frame(str_split(results_df$gene_id, ', ')[[i]])
      colnames(tmp) <- 'gene_id_panther'
      tmp2 <- merge(tmp, unique(transcript_gene[,c(3,5)]), by='gene_id_panther')
      tmp2 <- tmp2[match(tmp$gene_id_panther, tmp2$gene_id_panther),]
      colnames(tmp2) <- c('gene_id', 'gene_name')
      results_df[i,'gene_name'] <- paste(tmp2$gene_name, collapse=', ')
    }
    
    
  }
  
  write.xlsx(x=results_df[order(results_df$fdr), ], file=paste(prefix, '.xlsx', sep = ''), col.names = T, row.names = F)
  
  results_df_filtered <- results_df
  if(!is.null(toplog2fc)){
    if(!is.null(botlog2fc)){
      results_df_filtered <- subset(results_df, results_df$log2fold_enrichment > toplog2fc | results_df$log2fold_enrichment < botlog2fc)
    }
    else{
      results_df_filtered <- subset(results_df, results_df$log2fold_enrichment > toplog2fc)
    }
  }else{
    if(!is.null(botlog2fc)){
      results_df_filtered <- subset(results_df, results_df$log2fold_enrichment < botlog2fc)
    }
  }
  if(!is.null(morethan)){
    results_df_filtered <- results_df_filtered[which(unlist(lapply(str_split(results_df_filtered$gene_id, ', '), length)) > morethan),]
  }
  # results_df_filtered <- subset(results_df_filtered, results_df_filtered$pval < pvalue)
  
  results_df_filtered$GO_label <- as.character(results_df_filtered$GO_label)
  for(i in 1:nrow(results_df_filtered)){
    if(length(unlist(str_split(as.character(results_df_filtered[i,'GO_label']), " "))) > 10){
      newline <- ceiling(length(unlist(str_split(as.character(results_df_filtered[i,'GO_label']), " ")))/2)
      no_words <- floor(length(unlist(str_split(as.character(results_df_filtered[i,'GO_label']), " ")))) 
      results_df_filtered[i,'GO_label'] <- paste(c(unlist(str_split(as.character(results_df_filtered[i,'GO_label']), " "))[1:newline], "\n", unlist(str_split(as.character(results_df_filtered[i,'GO_label']), " "))[(newline+2):no_words]), collapse = " ")
    }  
  }
  results_df_filtered$fdr_log10 <- -log10(results_df_filtered$fdr)
  results_df_filtered$num_of_genes <- unlist(lapply(str_split(results_df_filtered$gene_id, ', '), FUN=length))
  
  steps <- ceiling((max(results_df_filtered$num_of_genes) - min(results_df_filtered$num_of_genes))/5)
  results_df_filtered$num_of_genes <- cut2(results_df_filtered$num_of_genes, seq(min(results_df_filtered$num_of_genes), max(results_df_filtered$num_of_genes), steps), digits=0, oneval=FALSE)
  
  results_df_filtered$num_of_genes <- sub(',', '-', gsub('^.|.$', '', as.character(results_df_filtered$num_of_genes)))
  
  library(RColorBrewer)
  getPalette = brewer.pal(length(unique(results_df_filtered$num_of_genes))+1, "Blues")[-1]
  results_df_filtered$GO_label <- factor(results_df_filtered$GO_label, levels=results_df_filtered[order(results_df_filtered$fdr, decreasing = T),'GO_label'])
  
  plot <- ggplot(data=results_df_filtered, aes(x=GO_label, y=log2fold_enrichment)) +
    geom_bar(stat="identity", aes(fill=num_of_genes))+
    expand_limits(y=c(min(results_df_filtered$log2fold_enrichment), (max(results_df_filtered$log2fold_enrichment)+0.5)))+
    xlab("GO term") + ylab("\nlog2(Fold Change)")+ theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip() +
    geom_hline(yintercept=0, color='black')+
    geom_point(aes(y=fdr_log10), size=3) + scale_y_continuous(sec.axis = sec_axis(~.*1, name = "-log10(padj)"))+
    geom_hline(yintercept=-log10(0.05), linetype='dashed', color='red')+ scale_fill_manual(values=getPalette) + labs(fill="Num. of genes") +
    theme(text = element_text(size=18), axis.text.x = element_text(angle=90, hjust=1)) 
  
  
  
  ggsave(plot, file=paste(prefix, '_fdr.svg', sep = ''), width=32, height=length(results_df_filtered$GO_label)*2, units="cm")
  ggsave(plot, file=paste(prefix, '_fdr.png', sep = ''), width=32, height=length(results_df_filtered$GO_label)*2, units="cm")
  
  return(plot)
}
## more_10 ##
# How many items in this row (vector) have a value greater than 10?
more_10 <- function(row){
  return(length(which(row > 10)))
}
## meanPlot_cus ##
# Create a mean plot. The function is mostly taken from the package deseqAbstraction which is no longer
# available on gitHub.
# Credits to Per L. Brattås
# The plot is generated with ggplot instead of basic R to add the posibility of having highlights or 
# labels.
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
# Annotations ----

# Read a parsed file from the original GTF file given by TEtranscripts authors to have the TE classification in hand
# The tabulated file includes TE id (e.g. "MMERVK10C-int_dup104"), TE subfamily (e.g. "MMERVK10C-int")
# TE family (e.g. "ERVK") and TE class (e.g. "LTR")
TE_classification <- fread('/Volumes/Seagate Backup /annotation/mouse/repeatmasker/mm10_rmsk_TE_classification.tab', data.table = FALSE, header=FALSE, fill=TRUE)
colnames(TE_classification) <- c('TE_id', 'TE_subfamily', 'TE_family', 'TE_class')

# Parsed tabulated file of gencode annotation vM20 including transcript ids, gene ids, gene name and biotype of the gene
transcript_gene <- fread('/Volumes/Seagate Backup /annotation/mouse/gencode/gencode.vM20.annotation.transc.gene.tab', data.table = FALSE, header=FALSE)
colnames(transcript_gene) <- c('transcript_id', 'gene_id', 'gene_name', 'gene_type')

# EMX animals ----
# Read output from TEtranscripts
emx <- fread('6_TEtranscripts/emx/emx.cntTable', data.table = F)
colnames(emx)[-1] <- paste(unlist(lapply(strsplit(colnames(emx)[-1], '/'), `[[`, 5)), unlist(lapply(strsplit(colnames(emx)[-1], '/'), `[[`, 4)), sep='_')
# Gene counts
emx_gene <- subset(emx, startsWith(emx$`gene/TE`, 'ENSMUS'))
colnames(emx_gene)[1] <- 'gene_id'
rownames(emx_gene) <- emx_gene$gene_id
emx_gene <- emx_gene[,-1]
# TE counts (where TE_id, in difference to TEclassification, is TE subfamily:TE family:TE class)
emx_TE <- subset(emx, !startsWith(emx$`gene/TE`, 'ENSMUS'))
colnames(emx_TE)[1] <- 'TE_id'
rownames(emx_TE) <- emx_TE$TE_id
emx_TE <- emx_TE[,-1]

# For this experiment we focus on cortical samples
emx_ctx <- colnames(emx_gene)[grepl('ctx', colnames(emx_gene))]

# Create the metadata 
emx_coldata <- data.frame(samples=emx_ctx, condition=unlist(lapply(strsplit(emx_ctx, '_'), `[[`, 2)))
rownames(emx_coldata) <- emx_coldata$samples
emx_coldata$samples <- as.character(emx_coldata$samples)

# Gene differential expression analysis testing for condition (Trim28 KO vs control)
emx_genes_dds <- DESeqDataSetFromMatrix(emx_gene[,rownames(emx_coldata)], emx_coldata, design = ~ condition)
emx_genes_dds <- DESeq(emx_genes_dds)
emx_genes_res <- results(emx_genes_dds)
# Calculate average expression on normalized reads on genes
emx_genes_exp <- getAverage(emx_genes_dds)
emx_genes_vst <- varianceStabilizingTransformation(emx_genes_dds)
# Calculate log2FC confidence intervals
emx_genes_res_df <- as.data.frame(emx_genes_res)
emx_genes_res_df$ci_low <- emx_genes_res_df$log2FoldChange - (qnorm(0.05)*emx_genes_res_df$lfcSE)
emx_genes_res_df$ci_high <- emx_genes_res_df$log2FoldChange + (qnorm(0.05)*emx_genes_res_df$lfcSE)

# Plot PCA based on gene expression (variance stabilized) 
p_gene_pca <- plotPCA(emx_genes_vst) + theme_classic() + ggtitle("PCA of gene expression")
ggsave(p_gene_pca, file="6_TEtranscripts/emx/plots/gene_pca_emx.png", width=20, height=20, units="cm", dpi=320)
ggsave(p_gene_pca, file="6_TEtranscripts/emx/plots/gene_pca_emx.svg", width=20, height=20, units="cm", dpi=320)

# Gene mean plot
p_gene_meanplot_emx <- meanPlot_cus(emx_genes_exp$Mean, test=emx_genes_res, l=0.5, p=0.05, c1='ko', c2='ctrl',ttl='Gene DEA in Cre Loxp experiment', repel = FALSE) + labs(title="", subtitle="")
ggsave(p_gene_meanplot_emx, file="6_TEtranscripts/emx/plots/gene_meanplot_0.5.png", width=20, height=20, units="cm", dpi=320)

# Boxplot showing expression of Trim28 (normalized expression)
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

ggsave(p_trim28_emx, file='6_TEtranscripts/emx/plots/trim28.svg', width=20, height=20, units="cm", dpi=96)

# Subset of upregulated genes to check for GO term enrichment
emx_upreg <- as.data.frame(as.matrix(subset(emx_genes_res, emx_genes_res$log2FoldChange > 0 & emx_genes_res$padj < 0.05)))
emx_upreg <- merge(emx_upreg, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
colnames(emx_upreg) <- c('Gene ID', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'P-adj', 'Gene Name')

write.xlsx(emx_upreg[,c(8,1,3,7)], '6_TEtranscripts/emx/upregulated_genes.xlsx', row.names = FALSE, col.names=TRUE)
write.table(unique(emx_upreg$`Gene Name`), 'GO_analysis/emx/upregulated/sign_upreg_genes.txt', row.names = FALSE, col.names=FALSE, quote = F)

# Subset of downregulated genes to check for GO term enrichment
emx_dwreg <- as.data.frame(as.matrix(subset(emx_genes_res, emx_genes_res$log2FoldChange < 0 & emx_genes_res$padj < 0.05)))
emx_dwreg <- merge(emx_dwreg, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
colnames(emx_dwreg) <- c('Gene ID', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'P-adj', 'Gene Name')

write.table(unique(emx_dwreg$`Gene Name`), 'GO_analysis/emx/downregulated/sign_downreg_genes.txt', row.names = FALSE, col.names=FALSE, quote=F)

# Subset of expressed genes to use as background on GO term enrichment tests
emx_gene_expressed <- emx_gene[which(apply(emx_gene, 1, more_10) > 0),]
write(rownames(emx_gene_expressed), 'GO_analysis/emx/background_emx.txt')

# emx_TE$TE_class <- sapply(str_split(rownames(emx_TE), ':'), `[[`, 3)
# emx_TE <- subset(emx_TE, emx_TE$TE_class %in% c('LINE', 'LTR', 'SINE'))
# TE differential expression analysis testing for differences between conditions (Trim28 KO vs control)
emx_TE_dds <- DESeqDataSetFromMatrix(emx_TE[,rownames(emx_coldata)], emx_coldata, design = ~ condition)
emx_TE_dds <- DESeq(emx_TE_dds)
emx_TE_res <- results(emx_TE_dds)
emx_TE_exp <- getAverage(emx_TE_dds)
emx_TE_vst <- varianceStabilizingTransformation(emx_TE_dds)

# TE mean plot
p_TE_meanplot_emx <- meanPlot_cus(emx_TE_exp$Mean, test=emx_TE_res, l=0.5, p=0.05, c1='ko', c2='ctrl',ttl='TE DEA in Cre Loxp experiment', repel = F, col3='firebrick', col2='black') + labs(title="", subtitle="")
ggsave(p_TE_meanplot_emx, file="6_TEtranscripts/emx/plots/TE_meanplot_0.5.png", width=20, height=20, units="cm", dpi=320)

# PCA plot based on TE expression (variance stabilized)
p_TE_pca <- plotPCA(emx_TE_vst) + theme_classic() + ggtitle("PCA of TE expression")
ggsave(p_TE_pca, file="6_TEtranscripts/emx/plots/TE_pca_emx.png", width=20, height=20, units="cm", dpi=320)
ggsave(p_TE_pca, file="6_TEtranscripts/emx/plots/TE_pca_emx.svg", width=20, height=20, units="cm", dpi=320)

# Normalize TE expression by the resulting sizeFactors of the gene expression data
emx_TE_norm <- emx_TE
emx_TE_norm[] <- mapply('/', emx_TE_norm, emx_genes_dds$sizeFactor)
# Add TE classification to just show retrotransposons
emx_TE_norm$TE_subfamily <- as.character(unlist(lapply(strsplit(rownames(emx_TE_norm), ':'), `[[`, 1)))
emx_TE_norm$TE_family <- unlist(lapply(strsplit(rownames(emx_TE_norm), ':'), `[[`, 2))
emx_TE_norm$TE_class <- unlist(lapply(strsplit(rownames(emx_TE_norm), ':'), `[[`, 3))
emx_TE_norm$TE_family <- as.character(ifelse(endsWith((emx_TE_norm$TE_family), '?'), substr(emx_TE_norm$TE_family, 1, (nchar(emx_TE_norm$TE_family)-1)),emx_TE_norm$TE_family))
emx_TE_norm$TE_class <- as.character(ifelse(endsWith((emx_TE_norm$TE_class), '?'), substr(emx_TE_norm$TE_class, 1, (nchar(emx_TE_norm$TE_class)-1)),emx_TE_norm$TE_class))
emx_TE_norm <- emx_TE_norm[which(emx_TE_norm$TE_class %in% c('LINE', 'LTR', 'SINE')),]

# Expression of significantly different TE subfamilies
emx_TE_signdiff_condition <- melt(emx_TE_exp$Mean)
emx_TE_signdiff_condition <- emx_TE_signdiff_condition[emx_TE_signdiff_condition$Var1 %in% rownames(subset(emx_TE_res, emx_TE_res$padj < 0.05)),]
emx_TE_signdiff_condition$value <- log2(emx_TE_signdiff_condition$value+0.5)
colnames(emx_TE_signdiff_condition) <- c('TE_subfamily', 'Condition', 'log2Mean')
emx_TE_signdiff_condition$TE_class <- unlist(lapply(str_split(emx_TE_signdiff_condition$TE_subfamily, ':'), `[[`, 3))
emx_TE_signdiff_condition$TE_family <- unlist(lapply(str_split(emx_TE_signdiff_condition$TE_subfamily, ':'), `[[`, 2))
emx_TE_signdiff_condition$TE_subfamily <- unlist(lapply(str_split(emx_TE_signdiff_condition$TE_subfamily, ':'), `[[`, 1))

# Mean expression per condition of significantly different TE subfamilies
emx_TE_signdiff_condition_mean <- emx_TE_exp$Mean[rownames(subset(emx_TE_res, emx_TE_res$padj < 0.05)),]
emx_TE_signdiff_condition_mean <- as.data.frame(emx_TE_signdiff_condition_mean)
emx_TE_signdiff_condition_mean$TE_subfamily <- unlist(lapply(str_split(rownames(emx_TE_signdiff_condition_mean), ":"), `[[`, 1))
emx_TE_signdiff_condition_mean$TE_family <- unlist(lapply(str_split(rownames(emx_TE_signdiff_condition_mean), ":"), `[[`, 2))
emx_TE_signdiff_condition_mean$TE_class <- unlist(lapply(str_split(rownames(emx_TE_signdiff_condition_mean), ":"), `[[`, 3))
rownames(emx_TE_signdiff_condition_mean) <- emx_TE_signdiff_condition_mean$TE_subfamily

TE_class_list <- colorRampPalette(rev(brewer.pal(n = 7, name = "GnBu")))(120)[seq(20,100,20)]
names(TE_class_list) <- unique(emx_TE_signdiff_condition_mean$TE_class)
TE_family_list <- colorRampPalette(rev(brewer.pal(n = 7, name = "Purples")))(140)[seq(20,120,20)]
names(TE_family_list) <- unique(emx_TE_signdiff_condition_mean$TE_family)

classification_colors = list(TE_class = TE_class_list, TE_family = TE_family_list)

p_signdiff_TE_emx_annotation <- pheatmap(log2(emx_TE_signdiff_condition_mean[,c(1,2), drop=F]+0.5), 
                                         cluster_cols = F, 
                                         cluster_rows = F,
                                         annotation_row = emx_TE_signdiff_condition_mean[, c('TE_class', 'TE_family')], 
                                         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                         annotation_colors=classification_colors,
                                         main="Significantly different TE subfamilies in Emx animals\nTrim28 (padj < 0.05, log2FC > 0)")
ggsave(p_signdiff_TE_emx_annotation, file='6_TEtranscripts/emx/plots/signdiff_TE_heatmap.svg', width=20, height=20, units="cm", dpi=96)

# Barplot of significantly upregulated TE subfamilies in Trim28 KO
emx_TE_res_df <-  as.data.frame(emx_TE_res)
emx_TE_res_df$TE_subfamily <- unlist(lapply(str_split(rownames(emx_TE_res_df), ':'), `[[`, 1))
emx_TE_signdiff_condition <- subset(emx_TE_signdiff_condition, emx_TE_signdiff_condition$Condition == 'ko')
emx_TE_signdiff_condition <- merge(emx_TE_signdiff_condition, emx_TE_res_df[,c('TE_subfamily', 'padj', 'log2FoldChange')], by='TE_subfamily')
emx_TE_signdiff_condition$FoldChange <- 2^(emx_TE_signdiff_condition$log2FoldChange)
emx_TE_signdiff_condition <- emx_TE_signdiff_condition[order(emx_TE_signdiff_condition$FoldChange, decreasing = T),]
emx_TE_signdiff_condition$TE_subfamily <- factor(emx_TE_signdiff_condition$TE_subfamily, levels=unique(emx_TE_signdiff_condition[order(emx_TE_signdiff_condition$FoldChange, decreasing = T),'TE_subfamily']))

p_signdiff_TE_emx <- ggplot(data=emx_TE_signdiff_condition, aes(x=TE_subfamily, y=FoldChange)) +   
  geom_bar(aes(width=0.7), position = "dodge", stat="identity", fill = "tomato2") + theme_classic()+labs(y="FoldChange", x="TE subfamily", fill="Condition")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
  geom_hline(yintercept=1, linetype="dashed")

ggsave(p_signdiff_TE_emx, file='6_TEtranscripts/emx/plots/signdiff_TE.svg', width=15, height=15, units="cm", dpi=96)

# emx_iap_norm <- subset(emx_TE_norm, startsWith(emx_TE_norm$TE_subfamily, 'IAP'))
# emx_iap_norm <- emx_iap_norm[, emx_coldata$samples]
# 
# emx_iap_norm_means <- merge(data.frame(ctrl_mean=rowMeans(emx_iap_norm[,subset(emx_coldata, emx_coldata$condition == 'ctrl')$samples])), data.frame(ko_mean=rowMeans(emx_iap_norm[,subset(emx_coldata, emx_coldata$condition == 'ko')$samples])), by='row.names')
# write.xlsx(emx_iap_norm_means, '/Volumes/Seagate Backup /trim28/09.10.19/6_TEtranscripts/emx/IAP_expression.xlsx', row.names = F)
# 
# emx_iap_norm <- subset(emx_TE_norm, startsWith(emx_TE_norm$TE_subfamily, 'IAP'))
# emx_iap_norm <- emx_iap_norm[, c(emx_coldata$samples, "TE_subfamily")]
# emx_iap_norm <- merge(melt(emx_iap_norm), emx_coldata, by.x='variable', by.y='samples')
# p <- ggplot(emx_iap_norm, aes(y=value, x=condition, fill=condition)) + geom_boxplot() + theme_classic() 
# p + facet_wrap( ~ TE_subfamily, scales="free")

# NPC ----
# Read output from TEtranscripts
npc <- fread('/Volumes/Seagate Backup /trim28/09.10.19/6_TEtranscripts/npc/npc.cntTable', data.table = F)
colnames(npc)[-1] <- paste(unlist(lapply(strsplit(colnames(npc)[-1], '/'), `[[`, 5)), unlist(lapply(strsplit(colnames(npc)[-1], '/'), `[[`, 4)), unlist(lapply(strsplit(colnames(npc)[-1], '/'), `[[`, 6)), sep='_')
# Gene counts
npc_gene <- subset(npc, startsWith(npc$`gene/TE`, 'ENSMUS'))
colnames(npc_gene)[1] <- 'gene_id'
rownames(npc_gene) <- npc_gene$gene_id
npc_gene <- npc_gene[,-1]
# TE counts (where TE_id, in difference to TEclassification, is TE subfamily:TE family:TE class)
npc_TE <- subset(npc, !startsWith(npc$`gene/TE`, 'ENSMUS'))
colnames(npc_TE)[1] <- 'TE_id'
rownames(npc_TE) <- npc_TE$TE_id
npc_TE <- npc_TE[,-1]

# Create the metadata 
npc_coldata <- data.frame(samples=colnames(npc_gene), condition=unlist(lapply(strsplit(colnames(npc_gene), '_'), `[[`, 2)))
rownames(npc_coldata) <- npc_coldata$samples
npc_coldata$samples <- as.character(npc_coldata$samples)

# Gene differential expression analysis testing for condition (Trim28 KO vs control)
npc_genes_dds <- DESeqDataSetFromMatrix(npc_gene[,rownames(npc_coldata)], npc_coldata, design = ~ condition)
npc_genes_dds <- DESeq(npc_genes_dds)
npc_genes_res <- results(npc_genes_dds)
# Calculate average expression on normalized reads on genes
npc_genes_exp <- getAverage(npc_genes_dds)

npc_genes_vst <- varianceStabilizingTransformation(npc_genes_dds)
npc_gene_pca <- plotPCA(npc_genes_vst) + ylim(c(-30,20)) + theme_classic() + ggtitle("PCA - In vitro CRISPR gene expression")
ggsave(npc_gene_pca, file="6_TEtranscripts/npc/plots/gene_pca.svg", width=20, height=20, units="cm", dpi=96)


# PCA plot based on protein coding genes
npc_gene_protein <- merge(npc_gene, unique(transcript_gene[,-1]), by.x='row.names', by.y='gene_id')
npc_gene_protein <- subset(npc_gene_protein, npc_gene_protein$gene_type == 'protein_coding')
colnames(npc_gene_protein)[1] <- 'gene_id'
npc_genes_vst_protein <- varianceStabilizingTransformation(npc_genes_dds[npc_gene_protein$gene_id, ])
npc_gene_pca_protein <- plotPCA(npc_genes_vst_protein) + ylim(c(-30,20)) + theme_classic() + ggtitle("PCA - In vitro CRISPR gene expression")
ggsave(npc_gene_pca_protein, file="6_TEtranscripts/npc/plots/gene_protein_pca.svg", width=20, height=20, units="cm", dpi=96)

# Mean plot based on gene expression
p_gene_meanplot_npc <- meanPlot_cus(npc_genes_exp$Mean, test=npc_genes_res, l=0.5, p=0.05, c1='ko', c2='ctrl',ttl='', repel = FALSE) + labs(title="", subtitle="")
ggsave(p_gene_meanplot_npc, file="6_TEtranscripts/npc/plots/gene_meanplot_0.5.png", width=20, height=20, units="cm", dpi=320)

# npc_TE$TE_class <- sapply(str_split(rownames(npc_TE), ':'), `[[`, 3)
# npc_TE <- subset(npc_TE, npc_TE$TE_class %in% c('LINE', 'LTR', 'SINE'))
# TE differential expression analysis testing for differences in condition (Trim28 vs control)
npc_TE_dds <- DESeqDataSetFromMatrix(npc_TE[,rownames(npc_coldata)], npc_coldata, design = ~ condition)
npc_TE_dds <- DESeq(npc_TE_dds)
npc_TE_res <- results(npc_TE_dds)
# Calculate average expression on normalized reads on TEs
npc_TE_exp <- getAverage(npc_TE_dds)
npc_TE_vst <- varianceStabilizingTransformation(npc_TE_dds)

# Plot PCA based on TE expression (variance stabilized) 
npc_TE_pca <- plotPCA(npc_TE_vst) + theme_classic() + ggtitle("PCA - In vitro CRISPR TE expression")
ggsave(npc_TE_pca, file="6_TEtranscripts/npc/plots/TE_pca.svg", width=20, height=20, units="cm", dpi=96)

# TE mean plot
p_TE_meanplot_npc <- meanPlot_cus(npc_TE_exp$Mean, test=npc_TE_res, col2 = 'black', col3 = 'firebrick', p=0.05, c1='ko', c2='ctrl',ttl='', repel = F, l=0.5)  + labs(title="", subtitle="")
ggsave(p_TE_meanplot_npc, file="6_TEtranscripts/npc/plots/TE_meanplot_0.5.png", width=20, height=20, units="cm", dpi=320)

# Boxplot showing expression of Trim28 (normalized expression)
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


ggsave(p_trim28_npc, file="6_TEtranscripts/npc/plots/trim28.svg", width=20, height=20, units="cm", dpi=96)

# Subset of upregulated genes to check for GO term enrichment
npc_upreg <- as.data.frame(as.matrix(subset(npc_genes_res, npc_genes_res$log2FoldChange > 0 & npc_genes_res$padj < 0.05)))
npc_upreg <- merge(npc_upreg, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
colnames(npc_upreg) <- c('Gene ID', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'P-adj', 'Gene Name')

write.xlsx(npc_upreg[,c(8,1,3,7)], '6_TEtranscripts/npc/upregulated_genes.xlsx', row.names = FALSE, col.names=TRUE)
write.table(unique(npc_upreg$`Gene Name`), 'GO_analysis/npc/upregulated/sign_upreg_genes.txt', row.names = FALSE, col.names=FALSE, quote = F)

# Subset of downregulated genes to check for GO term enrichment
npc_dwnreg <- as.data.frame(as.matrix(subset(npc_genes_res, npc_genes_res$log2FoldChange < 0 & npc_genes_res$padj < 0.05)))
npc_dwnreg <- merge(npc_dwnreg, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
colnames(npc_dwnreg) <- c('Gene ID', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'P-adj', 'Gene Name')

write.xlsx(npc_dwnreg[,c(8,1,3,7)], 'multimapping/6_TEtranscripts/npc/downregulated_genes.xlsx', row.names = FALSE, col.names=TRUE)
write.table(unique(npc_dwnreg$`Gene Name`), 'GO_analysis/npc/downregulated/sign_downreg_genes.txt', row.names = FALSE, col.names=FALSE, quote = F)

# Subset of expressed genes to use as background on GO term enrichment tests
npc_gene_expressed <- npc_gene[which(apply(npc_gene, 1, more_10) > 0),]
write(rownames(npc_gene_expressed), 'GO_analysis/npc/background_npc.txt')

# Normalize TE expression by the resulting sizeFactors of gene expression
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

# Mean TE expression per condition of signifiantly different TEs
npc_TE_signdiff_condition <- melt(npc_TE_exp$Mean)
npc_TE_signdiff_condition <- npc_TE_signdiff_condition[npc_TE_signdiff_condition$Var1 %in% rownames(subset(npc_TE_res, npc_TE_res$padj < 0.05)),]
npc_TE_signdiff_condition$value <- log2(npc_TE_signdiff_condition$value+0.5)
colnames(npc_TE_signdiff_condition) <- c('TE_subfamily', 'Condition', 'log2Mean')
npc_TE_signdiff_condition$TE_subfamily <- unlist(lapply(str_split(npc_TE_signdiff_condition$TE_subfamily, ":"), `[[`, 1))

npc_TE_signdiff_condition_mean <- npc_TE_exp$Mean[rownames(subset(npc_TE_res, npc_TE_res$padj < 0.05)),]
npc_TE_signdiff_condition_mean <- as.data.frame(npc_TE_signdiff_condition_mean)
npc_TE_signdiff_condition_mean$TE_subfamily <- unlist(lapply(str_split(rownames(npc_TE_signdiff_condition_mean), ":"), `[[`, 1))
npc_TE_signdiff_condition_mean$TE_family <- unlist(lapply(str_split(rownames(npc_TE_signdiff_condition_mean), ":"), `[[`, 2))
npc_TE_signdiff_condition_mean$TE_class <- unlist(lapply(str_split(rownames(npc_TE_signdiff_condition_mean), ":"), `[[`, 3))
rownames(npc_TE_signdiff_condition_mean) <- npc_TE_signdiff_condition_mean$TE_subfamily

TE_class_list <- c(LTR="#21B6A8")
TE_family_list <- colorRampPalette(rev(brewer.pal(n = 7, name = "Purples")))(100)[seq(10,80,20)]
names(TE_family_list) <- unique(npc_TE_signdiff_condition_mean$TE_family)

classification_colors = list(TE_class = TE_class_list, TE_family = TE_family_list)

p_signdiff_TE_npc_annotation <- pheatmap(log2(npc_TE_signdiff_condition_mean[,c(1,2), drop=F]+0.5), 
                                         cluster_cols = F, 
                                         cluster_rows = F,
                                         annotation_row = npc_TE_signdiff_condition_mean[, c('TE_class', 'TE_family')], 
                                         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
                                         annotation_colors=classification_colors,
                                         main="Significantly different TE subfamilies in npc animals\nTrim28 (padj < 0.05, log2FC > 0)")
ggsave(p_signdiff_TE_npc_annotation, file='6_TEtranscripts/npc/plots/signdiff_TE_heatmap.svg', width=20, height=25, units="cm", dpi=96)


# Barplot of significantly upregulated TE subfamilies in Trim28 KO
npc_TE_res_df <-  as.data.frame(npc_TE_res)
npc_TE_res_df$TE_subfamily <- unlist(lapply(str_split(rownames(npc_TE_res_df), ':'), `[[`, 1))
npc_TE_signdiff_condition <- subset(npc_TE_signdiff_condition, npc_TE_signdiff_condition$Condition == 'ko')
npc_TE_signdiff_condition <- merge(npc_TE_signdiff_condition, npc_TE_res_df[,c('TE_subfamily', 'padj', 'log2FoldChange')], by='TE_subfamily')
npc_TE_signdiff_condition$FoldChange <- 2^(npc_TE_signdiff_condition$log2FoldChange)
npc_TE_signdiff_condition <- npc_TE_signdiff_condition[order(npc_TE_signdiff_condition$FoldChange, decreasing = T),]

npc_TE_signdiff_condition$TE_subfamily <- factor(npc_TE_signdiff_condition$TE_subfamily, levels=unique(npc_TE_signdiff_condition[order(npc_TE_signdiff_condition$FoldChange, decreasing = T),'TE_subfamily']))

p_signdiff_TE_npc <- ggplot(data=npc_TE_signdiff_condition, aes(x=TE_subfamily, y=FoldChange)) +   
  geom_bar(aes(width=0.7), position = "dodge", stat="identity", fill = "tomato2") + theme_classic()+labs(y="FoldChange", x="TE subfamily", fill="Condition")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
  geom_hline(yintercept=1, linetype="dashed")

ggsave(p_signdiff_TE_npc, file='6_TEtranscripts/npc/plots/signdiff_TE.svg', width=25, height=15, units="cm", dpi=96)


# CRISPR invivo ---- 
# Read output from TEtranscripts
invivo_crispr <-  fread('6_TEtranscripts/invivo_crispr/invivo_crispr.cntTable', data.table = F)
colnames(invivo_crispr)[-1] <- paste(unlist(lapply(strsplit(colnames(invivo_crispr)[-1], '/'), `[[`, 5)), unlist(lapply(strsplit(colnames(invivo_crispr)[-1], '/'), `[[`, 4)), sep='_')

# Gene counts
invivo_crispr_gene <- subset(invivo_crispr, startsWith(invivo_crispr$`gene/TE`, 'ENSMUS'))
colnames(invivo_crispr_gene)[1] <- 'gene_id'
rownames(invivo_crispr_gene) <- invivo_crispr_gene$gene_id
invivo_crispr_gene <- invivo_crispr_gene[,-1]

# Create the metadata 
invivo_crispr_coldata <- data.frame(samples=colnames(invivo_crispr_gene), condition=unlist(lapply(strsplit(colnames(invivo_crispr_gene), '_'), `[[`, 3)))
invivo_crispr_coldata$type <- unlist(lapply(strsplit(colnames(invivo_crispr_gene), '_'), `[[`, 2))
rownames(invivo_crispr_coldata) <- invivo_crispr_coldata$samples
invivo_crispr_coldata$samples <- as.character(invivo_crispr_coldata$samples)

# Gene differential expression analysis testing for condition (Trim28 KO vs control)
invivo_crispr_genes_dds <- DESeqDataSetFromMatrix(invivo_crispr_gene[,rownames(invivo_crispr_coldata)], invivo_crispr_coldata, design = ~ condition)
invivo_crispr_genes_dds <- DESeq(invivo_crispr_genes_dds)
invivo_crispr_genes_res <- results(invivo_crispr_genes_dds)
# Calculate average expression on normalized reads on genes
invivo_crispr_genes_exp <- getAverage(invivo_crispr_genes_dds)

# Mean plot based on gene expression
p_gene_meanplot_invivocrispr <- meanPlot_cus(invivo_crispr_genes_exp$Mean, test=invivo_crispr_genes_res, l=0.5,p=0.05, c1='ko', c2='ctrl',ttl='Gene DEA in invivo CRISPR KO experiment', repel = FALSE) + labs(title="", subtitle="")
ggsave(p_gene_meanplot_invivocrispr, file="6_TEtranscripts/invivo_crispr/plots/gene_meanplot_0.5.png", width=20, height=20, units="cm", dpi=320)

# Boxplot showing expression of Trim28 (normalized expression)
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

# Subset of upregulated genes to check for GO term enrichment
invivo_crispr_upreg <- as.data.frame(as.matrix(subset(invivo_crispr_genes_res, invivo_crispr_genes_res$log2FoldChange > 0 & invivo_crispr_genes_res$padj < 0.05)))
invivo_crispr_upreg <- merge(invivo_crispr_upreg, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
colnames(invivo_crispr_upreg) <- c('Gene ID', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'P-adj', 'Gene Name')
write.xlsx(invivo_crispr_upreg[,c(8,1,3,7)], '6_TEtranscripts/invivo_crispr/upregulated_genes.xlsx', row.names = FALSE, col.names=TRUE)

# TE counts (where TE_id, in difference to TEclassification, is TE subfamily:TE family:TE class)
invivo_crispr_TE <- subset(invivo_crispr, !startsWith(invivo_crispr$`gene/TE`, 'ENSMUS'))
colnames(invivo_crispr_TE)[1] <- 'TE_id'
rownames(invivo_crispr_TE) <- invivo_crispr_TE$TE_id
invivo_crispr_TE <- invivo_crispr_TE[,-1]

invivo_crispr_TE$TE_id <- as.character(unlist(lapply(strsplit(rownames(invivo_crispr_TE), ':'), `[[`, 1)))
invivo_crispr_TE$TE_family <- unlist(lapply(strsplit(rownames(invivo_crispr_TE), ':'), `[[`, 2))
invivo_crispr_TE$TE_class <- unlist(lapply(strsplit(rownames(invivo_crispr_TE), ':'), `[[`, 3))
invivo_crispr_TE$TE_family <- as.character(ifelse(endsWith((invivo_crispr_TE$TE_family), '?'), substr(invivo_crispr_TE$TE_family, 1, (nchar(invivo_crispr_TE$TE_family)-1)),invivo_crispr_TE$TE_family))
invivo_crispr_TE$TE_class <- as.character(ifelse(endsWith((invivo_crispr_TE$TE_class), '?'), substr(invivo_crispr_TE$TE_class, 1, (nchar(invivo_crispr_TE$TE_class)-1)),invivo_crispr_TE$TE_class))
# invivo_crispr_TE <- subset(invivo_crispr_TE, invivo_crispr_TE$TE_class %in% c('LINE', 'LTR', 'SINE'))

# invivo_crispr_TE$TE_class <- sapply(str_split(rownames(invivo_crispr_TE), ':'), `[[`, 3)
# invivo_crispr_TE <- subset(invivo_crispr_TE, invivo_crispr_TE$TE_class %in% c('LINE', 'LTR', 'SINE'))
# TE differential expression analysis testing for differences in condition (Trim28 vs control)
invivo_crispr_TE_dds <- DESeqDataSetFromMatrix(invivo_crispr_TE[,names(invivo_crispr_genes_dds$sizeFactor)], invivo_crispr_coldata, design = ~ condition)
invivo_crispr_TE_dds <- DESeq(invivo_crispr_TE_dds)
invivo_crispr_TE_res <- results(invivo_crispr_TE_dds)
# Calculate average expression on normalized reads on TEs
invivo_crispr_TE_exp <- getAverage(invivo_crispr_TE_dds)

# TE mean plot
p_TE_meanplot_invivocrispr <- meanPlot_cus(invivo_crispr_TE_exp$Mean, test=invivo_crispr_TE_res, l=0.5, p=0.05, c1='ko', col2='black', col3='firebrick', c2='ctrl',ttl='TE subfamilies at the CRISPR KO experiment', repel = FALSE) + labs(title="", subtitle="")
ggsave(p_TE_meanplot_invivocrispr, file="6_TEtranscripts/invivo_crispr/plots/TE_meanplot_0.5.png", width=20, height=20, units="cm", dpi=320)

# Create the metadata for the Cas experiment
invivo_crispr_cas_coldata <- subset(invivo_crispr_coldata, startsWith(invivo_crispr_coldata$samples, 'cas') | invivo_crispr_coldata$condition == 'ctrl')
# TE DEA testing for differences in condition (Trim28 vs control) for the cas experiment
invivo_crispr_cas_TE_dds <- DESeqDataSetFromMatrix(invivo_crispr_TE[,rownames(invivo_crispr_cas_coldata)], invivo_crispr_cas_coldata, design = ~ condition)
invivo_crispr_cas_TE_dds <- DESeq(invivo_crispr_cas_TE_dds)
invivo_crispr_cas_TE_res <- results(invivo_crispr_cas_TE_dds)
invivo_crispr_cas_TE_exp <- getAverage(invivo_crispr_cas_TE_dds)
# TE mean plot (cas experiment)
p_TE_meanplot_invivocrispr_cas <- meanPlot_cus(invivo_crispr_cas_TE_exp$Mean, test=invivo_crispr_cas_TE_res, l=0.5,p=0.05, c1='ko', c2='ctrl',ttl='TE DEA in invivo CRISPR KO experiment - Cas', repel = FALSE, col2 = 'black') + labs(title="", subtitle="")

# Create the metadata for the FlexCas experiment
invivo_crispr_flexcas_coldata <- subset(invivo_crispr_coldata, startsWith(invivo_crispr_coldata$samples, 'flexcas') | invivo_crispr_coldata$condition == 'ctrl')
# TE DEA testing for differences in condition (Trim28 vs control) for the cas experiment
invivo_crispr_flexcas_TE_dds <- DESeqDataSetFromMatrix(invivo_crispr_TE[,rownames(invivo_crispr_flexcas_coldata)], invivo_crispr_flexcas_coldata, design = ~ condition)
invivo_crispr_flexcas_TE_dds <- DESeq(invivo_crispr_flexcas_TE_dds)
invivo_crispr_flexcas_TE_res <- results(invivo_crispr_flexcas_TE_dds)
invivo_crispr_flexcas_TE_exp <- getAverage(invivo_crispr_flexcas_TE_dds)
invivo_crispr_flexcas_TE_res_df <- as.data.frame(invivo_crispr_flexcas_TE_res)
View(invivo_crispr_TE[,rownames(invivo_crispr_flexcas_coldata)])
p_TE_meanplot_invivocrispr_flexcas <- meanPlot_cus(invivo_crispr_flexcas_TE_exp$Mean, test=invivo_crispr_flexcas_TE_res, l=0.5,p=0.05, c1='ko', c2='ctrl',ttl='TE DEA in invivo CRISPR KO experiment - flexcas', repel = F, col3 = 'firebrick3', col2 = 'black') + labs(title="", subtitle="")


# Adult invivo: FLOXED ----
# Read output from TEtranscripts
invivo_crispr_trim28fl <- fread('6_TEtranscripts/invivo_crispr_trim28fl/invivo_crispr_trim28fl.cntTable', data.table = F)
colnames(invivo_crispr_trim28fl)[-1] <- paste(unlist(lapply(strsplit(colnames(invivo_crispr_trim28fl)[-1], '/'), `[[`, 5)), unlist(lapply(strsplit(colnames(invivo_crispr_trim28fl)[-1], '/'), `[[`, 4)), sep='_')
colnames(invivo_crispr_trim28fl)[-1] <- paste(unlist(lapply(strsplit(colnames(invivo_crispr_trim28fl)[-1], '_'), `[[`, 1)), ifelse(unlist(lapply(strsplit(colnames(invivo_crispr_trim28fl)[-1], '_'), `[[`, 2)) == 'ctrl', 'ko', 'ctrl'), sep='_')
# Gene counts
invivo_crispr_trim28fl_gene <- subset(invivo_crispr_trim28fl, startsWith(invivo_crispr_trim28fl$`gene/TE`, 'ENSMUS'))
colnames(invivo_crispr_trim28fl_gene)[1] <- 'gene_id'
rownames(invivo_crispr_trim28fl_gene) <- invivo_crispr_trim28fl_gene$gene_id
invivo_crispr_trim28fl_gene <- invivo_crispr_trim28fl_gene[,-1]
# Create the metadata 
invivo_crispr_trim28fl_coldata <- data.frame(samples=colnames(invivo_crispr_trim28fl_gene), condition=unlist(lapply(strsplit(colnames(invivo_crispr_trim28fl_gene), '_'), `[[`, 2)))
rownames(invivo_crispr_trim28fl_coldata) <- invivo_crispr_trim28fl_coldata$samples
invivo_crispr_trim28fl_coldata$samples <- as.character(invivo_crispr_trim28fl_coldata$samples)

# Gene differential expression analysis testing for condition (Trim28 KO vs control)
invivo_crispr_trim28fl_genes_dds <- DESeqDataSetFromMatrix(invivo_crispr_trim28fl_gene[,rownames(invivo_crispr_trim28fl_coldata)], invivo_crispr_trim28fl_coldata, design = ~ condition)
invivo_crispr_trim28fl_genes_dds <- DESeq(invivo_crispr_trim28fl_genes_dds)
invivo_crispr_trim28fl_genes_res <- results(invivo_crispr_trim28fl_genes_dds)
# Calculate average expression on normalized reads on genes
invivo_crispr_trim28fl_genes_exp <- getAverage(invivo_crispr_trim28fl_genes_dds)

# Mean plot based on gene expression
p_gene_meanplot_invivoadult <- meanPlot_cus(invivo_crispr_trim28fl_genes_exp$Mean, test=invivo_crispr_trim28fl_genes_res, l=0.5,p=0.05, c1='ko', c2='ctrl',ttl='Gene DEA at floxed adult', repel = FALSE) + labs(title="", subtitle="")
ggsave(p_gene_meanplot_invivoadult, file="6_TEtranscripts/invivo_crispr_trim28fl/plots/gene_meanplot_0.5.png", width=20, height=20, units="cm", dpi=320)

# Boxplot showing expression of Trim28 (normalized expression)
invivo_crispr_trim28fl_gene_norm <- counts(invivo_crispr_trim28fl_genes_dds, normalized = TRUE)
invivo_crispr_trim28fl_gene_norm <- merge(invivo_crispr_trim28fl_gene_norm, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
invivo_crispr_trim28fl_gene_norm_trim28 <- as.data.frame(t(invivo_crispr_trim28fl_gene_norm[which(invivo_crispr_trim28fl_gene_norm$gene_name == 'Trim28'),as.character(invivo_crispr_trim28fl_coldata$samples)]))
colnames(invivo_crispr_trim28fl_gene_norm_trim28) <- 'value'
invivo_crispr_trim28fl_gene_norm_trim28$sample <- rownames(invivo_crispr_trim28fl_gene_norm_trim28)
invivo_crispr_trim28fl_gene_norm_trim28$Condition <- ifelse(invivo_crispr_trim28fl_coldata[as.character(invivo_crispr_trim28fl_gene_norm_trim28$sample),'condition'] == 'ko', 'Knock out', 'Control')

p_trim28_invivoadult <- ggplot(invivo_crispr_trim28fl_gene_norm_trim28, aes(x=Condition, y=value, fill=Condition)) + geom_boxplot()+ theme_classic() + labs(y="Median-of-ratios normalized read count", title='TRIM28 adult knock out')+
  theme(
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15,margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 15,margin = margin(t = 0, r = 20, b = 0, l = 0)))+ scale_y_continuous(limits = c(0, 3000))


ggsave(p_trim28_invivoadult, file="6_TEtranscripts/invivo_crispr_trim28fl/plots/trim28.svg", width=20, height=20, units="cm", dpi=96)

# Subset of upregulated genes to check for GO term enrichment
invivo_crispr_trim28fl_upreg <- as.data.frame(as.matrix(subset(invivo_crispr_trim28fl_genes_res, invivo_crispr_trim28fl_genes_res$log2FoldChange > 0 & invivo_crispr_trim28fl_genes_res$padj < 0.05)))
invivo_crispr_trim28fl_upreg <- merge(invivo_crispr_trim28fl_upreg, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
colnames(invivo_crispr_trim28fl_upreg) <- c('Gene ID', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'P-adj', 'Gene Name')
write.xlsx(invivo_crispr_trim28fl_upreg[,c(8,1,3,7)], '6_TEtranscripts/invivo_crispr_trim28fl/upregulated_genes.xlsx', row.names = FALSE, col.names=TRUE)

# Subset of downregulated genes to check for GO term enrichment
invivo_crispr_trim28fl_dwnreg <- as.data.frame(as.matrix(subset(invivo_crispr_trim28fl_genes_res, invivo_crispr_trim28fl_genes_res$log2FoldChange < 0 & invivo_crispr_trim28fl_genes_res$padj < 0.05)))
invivo_crispr_trim28fl_dwnreg <- merge(invivo_crispr_trim28fl_dwnreg, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
colnames(invivo_crispr_trim28fl_dwnreg) <- c('Gene ID', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'P-adj', 'Gene Name')
write.xlsx(invivo_crispr_trim28fl_dwnreg[,c(8,1,3,7)], '6_TEtranscripts/invivo_crispr_trim28fl/downregulated_genes.xlsx', row.names = FALSE, col.names=TRUE)

# TE counts (where TE_id, in difference to TEclassification, is TE subfamily:TE family:TE class)
invivo_crispr_trim28fl_TEs <- subset(invivo_crispr_trim28fl, !startsWith(invivo_crispr_trim28fl$`gene/TE`, 'ENSMUS'))
colnames(invivo_crispr_trim28fl_TEs)[1] <- 'TE_id'
rownames(invivo_crispr_trim28fl_TEs) <- invivo_crispr_trim28fl_TEs$TE_id
invivo_crispr_trim28fl_TEs <- invivo_crispr_trim28fl_TEs[,-1]

invivo_crispr_trim28fl_TEs$TE_id <- as.character(unlist(lapply(strsplit(rownames(invivo_crispr_trim28fl_TEs), ':'), `[[`, 1)))
invivo_crispr_trim28fl_TEs$TE_family <- unlist(lapply(strsplit(rownames(invivo_crispr_trim28fl_TEs), ':'), `[[`, 2))
invivo_crispr_trim28fl_TEs$TE_class <- unlist(lapply(strsplit(rownames(invivo_crispr_trim28fl_TEs), ':'), `[[`, 3))
invivo_crispr_trim28fl_TEs$TE_family <- as.character(ifelse(endsWith((invivo_crispr_trim28fl_TEs$TE_family), '?'), substr(invivo_crispr_trim28fl_TEs$TE_family, 1, (nchar(invivo_crispr_trim28fl_TEs$TE_family)-1)),invivo_crispr_trim28fl_TEs$TE_family))
invivo_crispr_trim28fl_TEs$TE_class <- as.character(ifelse(endsWith((invivo_crispr_trim28fl_TEs$TE_class), '?'), substr(invivo_crispr_trim28fl_TEs$TE_class, 1, (nchar(invivo_crispr_trim28fl_TEs$TE_class)-1)),invivo_crispr_trim28fl_TEs$TE_class))

# invivo_crispr_trim28fl_TEs <- subset(invivo_crispr_trim28fl_TEs, invivo_crispr_trim28fl_TEs$TE_class %in% c('LINE', 'LTR', 'SINE'))
# TE differential expression analysis testing for differences in condition (Trim28 vs control)
invivo_crispr_trim28fl_TEs_dds <- DESeqDataSetFromMatrix(invivo_crispr_trim28fl_TEs[,rownames(invivo_crispr_trim28fl_coldata)], invivo_crispr_trim28fl_coldata, design = ~ condition)
invivo_crispr_trim28fl_TEs_dds <- DESeq(invivo_crispr_trim28fl_TEs_dds)
invivo_crispr_trim28fl_TEs_res <- results(invivo_crispr_trim28fl_TEs_dds)
invivo_crispr_trim28fl_TEs_vst <- varianceStabilizingTransformation(invivo_crispr_trim28fl_TEs_dds)
invivo_crispr_trim28fl_TEs_exp <- getAverage(invivo_crispr_trim28fl_TEs_dds)

p_TE_meanplot_invivoadult <- meanPlot_cus(invivo_crispr_trim28fl_TEs_exp$Mean, test=invivo_crispr_trim28fl_TEs_res, p=0.05, l=0.5, c1='ko', c2='ctrl', col2='black', col3='firebrick', ttl='', repel = FALSE) + labs(title="", subtitle="")
ggsave(p_TE_meanplot_invivoadult, file="6_TEtranscripts/invivo_crispr_trim28fl/plots/TE_meanplot_0.5.png", width=20, height=20, units="cm", dpi=320)

# EMX dysregulated gene intersection with adult invivo floxed ----
emx_upreg <- as.data.frame(as.matrix(subset(emx_genes_res, emx_genes_res$log2FoldChange > 0.5 & emx_genes_res$padj < 0.05)))
emx_dwreg <- as.data.frame(as.matrix(subset(emx_genes_res, emx_genes_res$log2FoldChange < -0.5 & emx_genes_res$padj < 0.05)))

invivo_crispr_trim28fl_upreg <- as.data.frame(as.matrix(subset(invivo_crispr_trim28fl_genes_res, invivo_crispr_trim28fl_genes_res$log2FoldChange > 0.5 & invivo_crispr_trim28fl_genes_res$padj < 0.05)))
invivo_crispr_trim28fl_dwnreg <- as.data.frame(as.matrix(subset(invivo_crispr_trim28fl_genes_res, invivo_crispr_trim28fl_genes_res$log2FoldChange < -0.5 & invivo_crispr_trim28fl_genes_res$padj < 0.05)))

upreg_intersect <- table(rownames(invivo_crispr_trim28fl_upreg) %in% rownames(emx_upreg))["TRUE"]

dwreg_intersect <- table(rownames(invivo_crispr_trim28fl_dwnreg) %in% rownames(emx_dwreg))["TRUE"]

library(VennDiagram)
grid.newpage()
emx_invivo_crispr_trim28fl_upreg_venn <- draw.pairwise.venn(nrow(invivo_crispr_trim28fl_upreg), nrow(emx_upreg), upreg_intersect, category = c("AAV-Cre", "Emx1-Cre"), 
                                                            lty = rep("blank", 2), fill = c("light blue", "pink"), 
                                                            alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))

ggsave(emx_invivo_crispr_trim28fl_upreg_venn, file="/Volumes/Seagate Backup /trim28/09.10.19/6_TEtranscripts/emx_invivo_crispr_trim28fl_upreg_venn.svg", width=20, height=20, units="cm", dpi=320)

grid.newpage()
emx_invivo_crispr_trim28fl_dwreg_venn <- draw.pairwise.venn(nrow(invivo_crispr_trim28fl_dwnreg), nrow(emx_dwreg), dwreg_intersect, category = c("AAV-Cre", "Emx1-Cre"), 
                                                            lty = rep("blank", 2), fill = c("light blue", "pink"), 
                                                            alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
ggsave(emx_invivo_crispr_trim28fl_dwreg_venn, file="/Volumes/Seagate Backup /trim28/09.10.19/6_TEtranscripts/emx_invivo_crispr_trim28fl_dwreg_venn.svg", width=20, height=20, units="cm", dpi=320)

# EMX GO analysis ----
# Upregulated genes in Emx that are not upregulated in adult invivo crispr (floxed) KO 
emx_upreg <- merge(emx_upreg, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
colnames(emx_upreg)[1] <- 'Gene id'
colnames(emx_upreg)[ncol(emx_upreg)] <- 'Gene Name'
paste(emx_upreg[which(!emx_upreg$`Gene Name` %in% invivo_crispr_trim28fl_upreg$`Gene Name`),c(8)], collapse = '|')
write.xlsx(emx_upreg[which(!emx_upreg$`Gene Name` %in% invivo_crispr_trim28fl_upreg$`Gene Name`),c(8,1,3,7)], 'GO_analysis/emx/upregulated/upregulated_genes_emx_not_invivo_crispr_trim28fl.xlsx', col.names = T, row.names = F)
write(unique(emx_upreg[which(!emx_upreg$`Gene Name` %in% invivo_crispr_trim28fl_upreg$`Gene Name`),'Gene Name']), 'GO_analysis/emx/upregulated/sign_upreg_genes_not_in_adult.txt')

# Coordinates of those
gencode <- fread('/Volumes/Seagate Backup /annotation/mouse/gencode/gencode.vM20.annotation.bed', data.table = F)
colnames(gencode) <- c('chr', 'start', 'end', 'gene_id', 'strand')
rownames(gencode) <- gencode$gene_id
write.table(gencode[emx_upreg$`Gene ID`,], '6_TEtranscripts/emx/upregulated_genes_emx_not_invivo_crispr_trim28fl.bed', col.names = T, row.names = F, sep='\t', quote = F)

# Downregulated genes in Emx that are not downregulated in adult invivo crispr (floxed) KO 
write(unique(emx_dwreg[which(!emx_dwreg$`Gene Name` %in% invivo_crispr_trim28fl_dwnreg$`Gene Name`),'Gene Name']), 'GO_analysis/emx/downregulated/sign_downreg_genes_not_in_adult.txt')

# All Emx upregulated genes
emx_upreg_out <- merge(emx_upreg, unique(transcript_gene[,c('gene_id', 'gene_name')]), by.x='row.names', by.y='gene_id', all.x=T)
colnames(emx_upreg_out)[1] <- 'gene_id'
write.xlsx(emx_upreg_out, 'GO_analysis/emx/upregulated/upregulated_genes_emx.xlsx', col.names = T, row.names = F)
# All Emx downregulated genes
emx_dwreg_out <- merge(emx_dwreg, unique(transcript_gene[,c('gene_id', 'gene_name')]), by.x='row.names', by.y='gene_id', all.x=T)
colnames(emx_dwreg_out)[1] <- 'gene_id'
write.xlsx(emx_dwreg_out, 'GO_analysis/emx/downregulated/downregulated_genes_emx.xlsx', col.names = T, row.names = F)

# All upregulated genes at in vivo adult floxed KO
invivo_crispr_trim28fl_upreg_out <- merge(invivo_crispr_trim28fl_upreg, unique(transcript_gene[,c('gene_id', 'gene_name')]), by.x='row.names', by.y='gene_id', all.x=T)
colnames(invivo_crispr_trim28fl_upreg_out)[1] <- 'gene_id'
write.xlsx(invivo_crispr_trim28fl_upreg_out, 'GO_analysis/invivo_crispr_trim28fl/upregulated/upregulated_genes_invivo_crispr_trim28fl.xlsx', col.names = T, row.names = F)
# All downregulated genes at in vivo adult floxed KO
invivo_crispr_trim28fl_dwreg_out <- merge(invivo_crispr_trim28fl_dwnreg, unique(transcript_gene[,c('gene_id', 'gene_name')]), by.x='row.names', by.y='gene_id', all.x=T)
colnames(invivo_crispr_trim28fl_dwreg_out)[1] <- 'gene_id'
write.xlsx(invivo_crispr_trim28fl_dwreg_out, 'GO_analysis/invivo_crispr_trim28fl/downregulated/downregulated_genes_invivo_crispr_trim28fl.xlsx', col.names = T, row.names = F)

# EMX viral defence ----
viral_defence <- c("IFI16","IFI27","MS2","OAS1","IRF7","OASL","OAS2","OAS3","ISG20","MX1","IFIH1","IFIT3","IFI6","STAT1","IFIT2","ISG15","DDX58","DHX58","IFITM2","IFI35","B2M","IRF9","IFITM1","IFIT1","MX2")
viral_defence <- tools::toTitleCase(tolower(viral_defence))
viral_defence <- data.frame(gene_name=viral_defence)

viral_defence <- unique(merge(viral_defence, transcript_gene[,c(2,3)], by='gene_name'))

# Normalized expression of viral defence genes in Emx animals
emx_gene_norm_viral_defence <- emx_gene_norm[which(emx_gene_norm$gene_name %in% viral_defence$gene_name),c('gene_name', as.character(emx_coldata$samples))]
rownames(emx_gene_norm_viral_defence) <- emx_gene_norm_viral_defence$gene_name
# Normalized expression of viral defence genes in Emx animals with more than 10 reads (sum of all samples)
emx_gene_norm_viral_defence_more10 <- emx_gene_norm_viral_defence[which(rowSums(emx_gene_norm_viral_defence[,emx_coldata$samples]) > 10),]
emx_gene_norm_viral_defence_more10 <- merge(emx_gene_norm_viral_defence_more10, unique(transcript_gene[,c(2,3)]), by='gene_name')

# Log2 fold change confidence intervals of viral defence genes in Emx animals
emx_viral_defence_fc <- emx_genes_res_df[viral_defence$gene_id,c('log2FoldChange', 'ci_low', 'ci_high', 'padj'), drop=FALSE]
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
ggsave(p_viral_defence, file="6_TEtranscripts/emx/plots/viral_defence.svg", width=20, height=20, units="cm", dpi=96)

# EMX immune response ----
immune_response <- fread('/Volumes/Seagate Backup /trim28/09.10.19/GO_analysis/emx/inmmune_response_GO0006954.tab', data.table = F)
immune_response <- immune_response$V1
immune_response <- tools::toTitleCase(tolower(immune_response))
immune_response <- data.frame(gene_name=immune_response)

immune_response <- unique(merge(immune_response, transcript_gene[,c(2,3)], by='gene_name'))

# Normalized gene expression of inmmune response related genes
emx_gene_norm_immune_response <- emx_gene_norm[which(emx_gene_norm$gene_name %in% immune_response$gene_name),c('gene_name', as.character(emx_coldata$samples))]
rownames(emx_gene_norm_immune_response) <- emx_gene_norm_immune_response$gene_name
# Normalized gene expression of inmmune response related genes with more than 10 reads (sum of samples)
emx_gene_norm_immune_response_more10 <- emx_gene_norm_immune_response[which(rowSums(emx_gene_norm_immune_response[,emx_coldata$samples]) > 10),]
emx_gene_norm_immune_response_more10 <- merge(emx_gene_norm_immune_response_more10, unique(transcript_gene[,c(2,3)]), by='gene_name')

# Log2 confidence intervals of immune response related genes in Emx animals
emx_immune_response_fc <- emx_genes_res_df[immune_response$gene_id,c('log2FoldChange', 'ci_low', 'ci_high', 'padj'), drop=FALSE]
emx_immune_response_fc <- merge(emx_immune_response_fc, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
emx_immune_response_fc$reg <- ifelse(emx_immune_response_fc$log2FoldChange < 0, 'logFC < 0', 'logFC > 0')
emx_immune_response_fc$gene_name <- factor(emx_immune_response_fc$gene_name, levels=emx_immune_response_fc[order(emx_immune_response_fc$log2FoldChange, decreasing=TRUE),'gene_name'])

emx_immune_response_fc <- emx_immune_response_fc[which(emx_immune_response_fc$Row.names %in% emx_gene_norm_immune_response_more10$gene_id),]
emx_immune_response_fc_goodconfint <- subset(emx_immune_response_fc, (emx_immune_response_fc$ci_low > 0 & emx_immune_response_fc$ci_high > 0) | (emx_immune_response_fc$ci_low < 0 & emx_immune_response_fc$ci_high < 0))

p_immune_response <- ggplot() + 
  geom_errorbar(data=emx_immune_response_fc_goodconfint, mapping=aes(x=gene_name, ymin=ci_low, ymax=ci_high, colour=reg), width=0.2) + 
  geom_point(stat='identity',data=emx_immune_response_fc_goodconfint, mapping=aes(x=gene_name, y=log2FoldChange, colour=reg)) + coord_flip()+theme_bw() + scale_colour_manual(values=c("steelblue", 'tomato')) +
  guides(fill=FALSE)+
  labs(x="Gene name", y= 'log2 Fold Change', title='Change of expression on immune response genes at CreLoxP Trim28 KO\nError lines of 95% confidence intervals', colour="KO Relative expression")+
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

p_immune_response_complete <- ggplot() + 
  geom_errorbar(data=emx_immune_response_fc, mapping=aes(x=gene_name, ymin=ci_low, ymax=ci_high, colour=reg), width=0.2) + 
  geom_point(stat='identity',data=emx_immune_response_fc, mapping=aes(x=gene_name, y=log2FoldChange, colour=reg)) + coord_flip()+theme_bw() + scale_colour_manual(values=c("steelblue", 'tomato')) +
  guides(fill=FALSE)+
  labs(x="Gene name", y= 'log2 Fold Change', title='Change of expression on immune response genes at CreLoxP Trim28 KO\nError lines of 95% confidence intervals', colour="KO Relative expression")+
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

ggsave(p_immune_response_complete, file="6_TEtranscripts/emx/plots/immune_response_complete.svg", width=40, height=100, units="cm", dpi=96)
ggsave(p_immune_response, file="6_TEtranscripts/emx/plots/immune_response.svg", width=30, height=20, units="cm", dpi=96)

# EMX microglia ----
microglia <- c('Cd9', 'Cd11c',
               'Clec7a', 'Cd63',
               'Cst7', 'Csf1',
               'B2m', 'Bcl2',
               'Cd68', 'Hif1a',
               'C1qa', 'Ank',
               'Lgals3', 'Ly86',
               'Lamp1', 'Ccl2',
               'Ccl6', 'Rela',
               'Tnf', 'Lc3',
               'Myd88', 'Tlr1',
               'Tlr2', 'Tlr4',
               'Lamtor1', 'Mapk14',
               'Ncf1', 'Rage',
               'Sox1', 'Sox2',
               'Jak1',
               'Trim28')

microglia <- data.frame(gene_name=microglia)
microglia <- unique(merge(microglia, transcript_gene[,c(2,3)], by='gene_name'))

# Normalized expression of microglia related genes in Emx animals
emx_gene_norm_microglia <- emx_gene_norm[which(emx_gene_norm$gene_name %in% microglia$gene_name),c('gene_name', as.character(emx_coldata$samples))]
rownames(emx_gene_norm_microglia) <- emx_gene_norm_microglia$gene_name
# Normalized expression of microglia related genes in Emx animals with more than 10 reads (sum of samples)
emx_gene_norm_microglia_more10 <- emx_gene_norm_microglia[which(rowSums(emx_gene_norm_microglia[,emx_coldata$samples]) > 10),]
emx_gene_norm_microglia_more10 <- merge(emx_gene_norm_microglia_more10, unique(transcript_gene[,c(2,3)]), by='gene_name')

# Log2 confidence intervals of microglia activation genes
emx_genes_res_df <- as.data.frame(emx_genes_res)
emx_genes_res_df$ci_low <- emx_genes_res_df$log2FoldChange - (qnorm(0.05)*emx_genes_res_df$lfcSE)
emx_genes_res_df$ci_high <- emx_genes_res_df$log2FoldChange + (qnorm(0.05)*emx_genes_res_df$lfcSE)

emx_microglia_fc <- emx_genes_res_df[microglia$gene_id,c('log2FoldChange', 'ci_low', 'ci_high'), drop=FALSE]
emx_microglia_fc <- merge(emx_microglia_fc, unique(transcript_gene[,c(2,3)]), by.x='row.names', by.y='gene_id')
emx_microglia_fc$reg <- ifelse(emx_microglia_fc$log2FoldChange < 0, 'logFC < 0', 'logFC > 0')
emx_microglia_fc$gene_name <- factor(emx_microglia_fc$gene_name, levels=emx_microglia_fc[order(emx_microglia_fc$log2FoldChange, decreasing=TRUE),'gene_name'])

emx_microglia_fc <- emx_microglia_fc[which(emx_microglia_fc$Row.names %in% emx_gene_norm_microglia_more10$gene_id),]

png('multimapping/plots/emx/inflammatory_genes_fc.png', width = 12, height = 7, units = 'in', res=200)
ggplot() + 
  geom_errorbar(data=emx_microglia_fc, mapping=aes(x=gene_name, ymin=ci_low, ymax=ci_high, colour=reg), width=0.2) + 
  geom_point(stat='identity',data=emx_microglia_fc, mapping=aes(x=gene_name, y=log2FoldChange, colour=reg)) + coord_flip()+theme_bw() + scale_colour_manual(values=c("steelblue", 'tomato')) +
  guides(fill=FALSE)+
  labs(x="Gene name", y= 'log2 Fold Change', title='Change of expression on inflammatory genes at CreLoxP Trim28 KO\nError lines of 95% confidence intervals', colour="KO Relative expression")+
  theme(
    plot.title = element_text(size = 18),
    axis.text.x = element_text(size = 15),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15,margin = margin(t = 20, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 15,margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  geom_hline(yintercept = 0, color = "black", size=0.5)
dev.off()

emx_gene_norm_microglia_more10 <- melt(emx_gene_norm_microglia_more10, by='gene_name')
emx_gene_norm_microglia_more10$Condition <- ifelse(emx_coldata[as.character(emx_gene_norm_microglia_more10$variable),'condition'] == 'ko', 'Knock out', 'Control')

emx_gene_norm_microglia_more10$gene_name <- factor(emx_gene_norm_microglia_more10$gene_name, levels=emx_microglia_fc[order(emx_microglia_fc$log2FoldChange, decreasing=TRUE),'gene_name'])
emx_gene_norm_microglia_more10 <- subset(emx_gene_norm_microglia_more10, !is.na(emx_gene_norm_microglia_more10$gene_name))


# EMX GO Biological process ----
GO_json(file_name='/Volumes/Seagate Backup /trim28/GO_analysis/invivo_brain_dev/upregulated/slim_biological_process/not_in_adult/slim_biological_process_upreg_emx_not_in_adult.json', prefix='/Volumes/Seagate Backup /trim28/GO_analysis/invivo_brain_dev/upregulated/slim_biological_process/slim_biological_process_upreg_emx', morethan = 4,  botlog2fc = -0.5, toplog2fc = 0.5)
GO_json(file_name='/Volumes/Seagate Backup /trim28/GO_analysis/invivo_brain_dev/downregulated/slim_biological_process/slim_biological_process_downreg_emx.json', prefix='/Volumes/Seagate Backup /trim28/GO_analysis/invivo_brain_dev/downregulated/slim_biological_process/slim_biological_process_downreg_emx', morethan = 4,  botlog2fc = -0.5, toplog2fc = 0.5)

# EMX GO Biological process without genes that are also dysregulated in adult invivo crispr ----
GO_json(file_name='/Volumes/Seagate Backup /trim28/GO_analysis/invivo_brain_dev/upregulated/slim_biological_process/not_in_adult/slim_biological_process_upreg_emx_not_in_adult.json', prefix='/Volumes/Seagate Backup /trim28/GO_analysis/invivo_brain_dev/upregulated/slim_biological_process/not_in_adult/slim_biological_process_upreg_emx_not_in_adult', morethan = 4, botlog2fc = -0.5, toplog2fc = 0.5)
GO_json(file_name='/Volumes/Seagate Backup /trim28/GO_analysis/invivo_brain_dev/downregulated/slim_biological_process/not_in_adult/slim_biological_process_downreg_emx_not_in_adult.json', prefix='/Volumes/Seagate Backup /trim28/GO_analysis/invivo_brain_dev/downregulated/slim_biological_process/not_in_adult/slim_biological_process_downreg_emx_not_in_adult', morethan = 4, botlog2fc = -0.5, toplog2fc = 0.5)

# EMX GO Molecular function ----
GO_json(file_name='/Volumes/Seagate Backup /trim28/GO_analysis/invivo_brain_dev/upregulated/slim_molecular_function/slim_molecular_function_upreg_emx.json', prefix='/Volumes/Seagate Backup /trim28/GO_analysis/invivo_brain_dev/upregulated/slim_molecular_function/slim_molecular_function_upreg_emx', morethan = 4, botlog2fc = -0.5, toplog2fc = 0.5)
View(GO_json(file_name='/Volumes/Seagate Backup /trim28/GO_analysis/invivo_brain_dev/downregulated/slim_molecular_function/slim_molecular_function_downreg_emx.json', prefix='/Volumes/Seagate Backup /trim28/GO_analysis/invivo_brain_dev/downregulated/slim_molecular_function/slim_molecular_function_downreg_emx', morethan = 4, botlog2fc = -0.5, toplog2fc = 0.5))

# EMX GO Molecular function without genes that are also dysregulated in adult invivo crispr ----
GO_json(file_name='/Volumes/Seagate Backup /trim28/GO_analysis/invivo_brain_dev/upregulated/slim_molecular_function/not_in_adult/slim_molecular_function_upreg_emx_not_in_adult.json', prefix='/Volumes/Seagate Backup /trim28/GO_analysis/invivo_brain_dev/upregulated/slim_molecular_function/not_in_adult/slim_molecular_function_upreg_emx_not_in_adult', morethan = 4, toplog2fc = 0.5)
GO_json(file_name='/Volumes/Seagate Backup /trim28/GO_analysis/invivo_brain_dev/downregulated/slim_molecular_function/not_in_adult/slim_molecular_function_downreg_emx_not_in_adult.json', prefix='/Volumes/Seagate Backup /trim28/GO_analysis/invivo_brain_dev/downregulated/slim_molecular_function/not_in_adult/slim_molecular_function_downreg_emx_not_in_adult', morethan = 4, botlog2fc = -0.5, toplog2fc = 0.5)

# EMX MMERVK10C upregulated ----
# Read quantification of full length MMERVK10C elements (unique mapping)
MMERVK10C_count <- fread('/Volumes/Seagate Backup /trim28/09.10.19/1_uniqmapping/2_readcount/fullMMERVK10C_count_matrix_2.csv', data.table = F)
# Create the metadata 
MMERVK10C_coldata <- data.frame(sample=paste(sapply(str_split(colnames(MMERVK10C_count)[7:ncol(MMERVK10C_count)], '/'), `[[`, 5), sapply(str_split(colnames(MMERVK10C_count)[7:ncol(MMERVK10C_count)], '/'), `[[`, 6), sep='_'),
                                experiment=sapply(str_split(colnames(MMERVK10C_count)[7:ncol(MMERVK10C_count)], '/'), `[[`, 3),
                                condition=sapply(str_split(colnames(MMERVK10C_count)[7:ncol(MMERVK10C_count)], '/'), `[[`, 4))
MMERVK10C_coldata$sample <- as.character(MMERVK10C_coldata$sample)
MMERVK10C_coldata$sample <- ifelse(MMERVK10C_coldata$experiment == 'invivo_crispr', 
                                   paste(sapply(str_split(MMERVK10C_coldata$sample, '_'),`[[`, 1), sapply(str_split(MMERVK10C_coldata$sample, '_'),`[[`, 2), sep='_'),
                                   ifelse(MMERVK10C_coldata$experiment != 'npc', sapply(str_split(MMERVK10C_coldata$sample, '_'),`[[`, 1), MMERVK10C_coldata$sample))

colnames(MMERVK10C_count)[7:ncol(MMERVK10C_count)] <- ifelse(sapply(str_split(colnames(MMERVK10C_count)[7:ncol(MMERVK10C_count)], '/'), `[[`, 3) == 'invivo_crispr', 
                                                             paste(sapply(str_split(paste(sapply(str_split(colnames(MMERVK10C_count)[7:ncol(MMERVK10C_count)], '/'), `[[`, 5), sapply(str_split(colnames(MMERVK10C_count)[7:ncol(MMERVK10C_count)], '/'), `[[`, 6), sep='_'), '_'),`[[`, 1), sapply(str_split(paste(sapply(str_split(colnames(MMERVK10C_count)[7:ncol(MMERVK10C_count)], '/'), `[[`, 5), sapply(str_split(colnames(MMERVK10C_count)[7:ncol(MMERVK10C_count)], '/'), `[[`, 6), sep='_'), '_'),`[[`, 2), sep='_'),
                                                             ifelse(sapply(str_split(colnames(MMERVK10C_count)[7:ncol(MMERVK10C_count)], '/'), `[[`, 3) != 'npc', sapply(str_split(paste(sapply(str_split(colnames(MMERVK10C_count)[7:ncol(MMERVK10C_count)], '/'), `[[`, 5), sapply(str_split(colnames(MMERVK10C_count)[7:ncol(MMERVK10C_count)], '/'), `[[`, 6), sep='_'), '_'),`[[`, 1), paste(sapply(str_split(colnames(MMERVK10C_count)[7:ncol(MMERVK10C_count)], '/'), `[[`, 5), sapply(str_split(colnames(MMERVK10C_count)[7:ncol(MMERVK10C_count)], '/'), `[[`, 6), sep='_')))

rownames(MMERVK10C_coldata) <- MMERVK10C_coldata$sample
rownames(MMERVK10C_count) <- MMERVK10C_count$Geneid

# Subset quantification of NPCs
MMERVK10C_count_npc <- MMERVK10C_count[,c('Chr', 'Start', 'End', 'Geneid', 'Strand', subset(MMERVK10C_coldata, MMERVK10C_coldata$experiment == 'npc')$sample)]
MMERVK10C_coldata_npc <- subset(MMERVK10C_coldata, MMERVK10C_coldata$experiment == 'npc')
# Differential MMERVK10C expression analysis testing for difference in condition in NPCs (Trim28 KO vs control)
MMERVK10C_dds_npc <- DESeqDataSetFromMatrix(MMERVK10C_count_npc[,rownames(MMERVK10C_coldata_npc)], MMERVK10C_coldata_npc, design = ~ condition)
MMERVK10C_dds_npc <- DESeq(MMERVK10C_dds_npc)
MMERVK10C_res_npc <- results(MMERVK10C_dds_npc)

# Subset of upregulated MMERVK10C elements
MMERVK10C_upreg_npc <- rownames(MMERVK10C_res_npc[which(MMERVK10C_res_npc$log2FoldChange > 0),])
MMERVK10C_upreg_npc <- MMERVK10C_upreg_npc[which(MMERVK10C_upreg_npc %in% rownames(MMERVK10C_count_npc))]
write.table(MMERVK10C_count_npc[MMERVK10C_upreg_npc,c('Chr', 'Start', 'End', 'Geneid', 'Strand')], col.names = F, row.names = F, quote = F, sep='\t', file='/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/upregulated/npc_upreg_MMERVK10C.bed')
# Subset of not upregulated MMERVK10C elements (downregulated and not dysregulated)
# This is to test if there is a difference in nearby gene expression between the genes close to
# upregulated MMERVK10C elements and the not dysregulated MMERVK10C
MMERVK10C_not_upreg_npc <- rownames(MMERVK10C_res_npc[which(!MMERVK10C_res_npc$log2FoldChange > 0 | is.na(MMERVK10C_res_npc$log2FoldChange)),])
MMERVK10C_not_upreg_npc <- MMERVK10C_not_upreg_npc[which(MMERVK10C_not_upreg_npc %in% rownames(MMERVK10C_res_npc))]
write.table(MMERVK10C_count_npc[MMERVK10C_not_upreg_npc,c('Chr', 'Start', 'End', 'Geneid', 'Strand')], col.names = F, row.names = F, quote = F, sep='\t', file='/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/not_upregulated/npc_notupreg_MMERVK10C.bed')

# Subset quantification of Emx animals
MMERVK10C_count_emx <- MMERVK10C_count[,c('Chr', 'Start', 'End', 'Geneid', 'Strand', subset(MMERVK10C_coldata, MMERVK10C_coldata$experiment == 'emx')$sample)]
MMERVK10C_coldata_emx <- subset(MMERVK10C_coldata, MMERVK10C_coldata$experiment == 'emx')
# Differential MMERVK10C expression analysis testing for difference in condition in Emx (Trim28 KO vs control)
MMERVK10C_dds_emx <- DESeqDataSetFromMatrix(MMERVK10C_count_emx[,rownames(MMERVK10C_coldata_emx)], MMERVK10C_coldata_emx, design = ~ condition)
MMERVK10C_dds_emx <- DESeq(MMERVK10C_dds_emx)
MMERVK10C_res_emx <- results(MMERVK10C_dds_emx)

# Subset of upregulated MMERVK10C elements
MMERVK10C_upreg_emx <- rownames(MMERVK10C_res_emx[which(MMERVK10C_res_emx$log2FoldChange > 0),])
MMERVK10C_upreg_emx <- MMERVK10C_upreg_emx[which(MMERVK10C_upreg_emx %in% rownames(MMERVK10C_count_emx))]
write.table(MMERVK10C_count_emx[MMERVK10C_upreg_emx,c('Chr', 'Start', 'End', 'Geneid', 'Strand')], col.names = F, row.names = F, quote = F, sep='\t', file='/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/upregulated/emx_upreg_MMERVK10C.bed')
# Subset of not upregulated MMERVK10C elements (downregulated and not dysregulated)
# Same idea : This to test if there is a difference in nearby gene expression between the genes close to
# upregulated MMERVK10C elements and the not dysregulated MMERVK10C
MMERVK10C_not_upreg_emx <- rownames(MMERVK10C_res_emx[which(!MMERVK10C_res_emx$log2FoldChange > 0 | is.na(MMERVK10C_res_emx$log2FoldChange)),])
MMERVK10C_not_upreg_emx <- MMERVK10C_not_upreg_emx[which(MMERVK10C_not_upreg_emx %in% rownames(MMERVK10C_res_emx))]
write.table(MMERVK10C_count_emx[MMERVK10C_not_upreg_emx,c('Chr', 'Start', 'End', 'Geneid', 'Strand')], col.names = F, row.names = F, quote = F, sep='\t', file='/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/not_upregulated/emx_notupreg_MMERVK10C.bed')

# EMX - MMERVK10C nearby gene expression ----
# Read intersection of upregulated MMERVK10C elements and genes in 10kb, 25kb and 50kb distance
emx_upreg10kb <- fread('/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/upregulated/emx/emx_upreg_MMERVK10C_10kb_windows_intersect_genes.bed', data.table = F, header = F)
emx_upreg25kb <- fread('/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/upregulated/emx/emx_upreg_MMERVK10C_25kb_windows_intersect_genes.bed', data.table = F, header = F)
emx_upreg50kb <- fread('/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/upregulated/emx/emx_upreg_MMERVK10C_50kb_windows_intersect_genes.bed', data.table = F, header = F)

# Are these genes dysregulated?
emx_upreg10kb_res <- subset(emx_genes_res, rownames(emx_genes_res) %in% as.character(emx_upreg10kb$V4))
emx_upreg25kb_res <- subset(emx_genes_res, rownames(emx_genes_res) %in% as.character(emx_upreg25kb$V4))
emx_upreg50kb_res <- subset(emx_genes_res, rownames(emx_genes_res) %in% as.character(emx_upreg50kb$V4))
emx_nearby_genes_upregMMERVK10C <- data.frame(log2FC=emx_upreg10kb_res[,'log2FoldChange', drop=F], set=rep('10kb', nrow(emx_upreg10kb_res)))
emx_nearby_genes_upregMMERVK10C <- rbind(emx_nearby_genes_upregMMERVK10C, data.frame(log2FC=emx_upreg25kb_res[,'log2FoldChange', drop=F], set=rep('25kb', nrow(emx_upreg25kb_res))))
emx_nearby_genes_upregMMERVK10C <- rbind(emx_nearby_genes_upregMMERVK10C, data.frame(log2FC=emx_upreg50kb_res[,'log2FoldChange', drop=F], set=rep('50kb', nrow(emx_upreg50kb_res))))
emx_nearby_genes_upregMMERVK10C$disregulation <- rep("Upregulated", nrow(emx_nearby_genes_upregMMERVK10C))

# Read intersection of NOT upregulated MMERVK10C elements and genes in 10kb, 25kb and 50kb distance
emx_notupreg10kb <- fread('/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/not_upregulated/emx/emx_notupreg_MMERVK10C_10kb_windows_intersect_genes.bed', data.table = F, header = F)
emx_notupreg25kb <- fread('/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/not_upregulated/emx/emx_notupreg_MMERVK10C_25kb_windows_intersect_genes.bed', data.table = F, header = F)
emx_notupreg50kb <- fread('/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/not_upregulated/emx/emx_notupreg_MMERVK10C_50kb_windows_intersect_genes.bed', data.table = F, header = F)

# Are these genes dysregulated?
emx_notupreg10kb_res <- subset(emx_genes_res, rownames(emx_genes_res) %in% as.character(emx_notupreg10kb$V4))
emx_notupreg25kb_res <- subset(emx_genes_res, rownames(emx_genes_res) %in% as.character(emx_notupreg25kb$V4))
emx_notupreg50kb_res <- subset(emx_genes_res, rownames(emx_genes_res) %in% as.character(emx_notupreg50kb$V4))
emx_nearby_genes_notupregMMERVK10C <- data.frame(log2FC=emx_notupreg10kb_res[,'log2FoldChange', drop=F], set=rep('10kb', nrow(emx_notupreg10kb_res)))
emx_nearby_genes_notupregMMERVK10C <- rbind(emx_nearby_genes_notupregMMERVK10C, data.frame(log2FC=emx_notupreg25kb_res[,'log2FoldChange', drop=F], set=rep('25kb', nrow(emx_notupreg25kb_res))))
emx_nearby_genes_notupregMMERVK10C <- rbind(emx_nearby_genes_notupregMMERVK10C, data.frame(log2FC=emx_notupreg50kb_res[,'log2FoldChange', drop=F], set=rep('50kb', nrow(emx_notupreg50kb_res))))
emx_nearby_genes_notupregMMERVK10C$disregulation <- rep("Not Upregulated", nrow(emx_nearby_genes_notupregMMERVK10C))

emx_nearby_genes_MMERVK10C <- rbind(emx_nearby_genes_upregMMERVK10C, emx_nearby_genes_notupregMMERVK10C)

library(ggpubr)

emx_nearby_genes_MMERVK10C_plot <- ggplot(emx_nearby_genes_MMERVK10C, aes(y=log2FoldChange, x=set, fill=disregulation)) + 
  geom_boxplot() + theme_classic() + scale_fill_brewer(palette="Set2") + 
  labs(fill = "Distance to\nfull MMERVK10C", x="") + ggtitle("Emx - Genes nearby upregulated MMERVK10Cs")+ 
  ylim(c(-3,6)) + stat_compare_means(method ='t.test')

ggsave(emx_nearby_genes_MMERVK10C_plot, file='/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/emx/emx_nearby_gene_MMERVK10C.svg', width=20, height=20, units="cm", dpi=96)

# NPC - MMERVK10C nearby gene expression ----
# Read intersection of upregulated MMERVK10C elements and genes in 10kb, 25kb and 50kb distance
npc_upreg10kb <- fread('/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/upregulated/npc/npc_upreg_MMERVK10C_10kb_windows_intersect_genes.bed', data.table = F, header = F)
npc_upreg25kb <- fread('/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/upregulated/npc/npc_upreg_MMERVK10C_25kb_windows_intersect_genes.bed', data.table = F, header = F)
npc_upreg50kb <- fread('/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/upregulated/npc/npc_upreg_MMERVK10C_50kb_windows_intersect_genes.bed', data.table = F, header = F)

# Are these genes dysregulated?
npc_upreg10kb_res <- subset(npc_genes_res, rownames(npc_genes_res) %in% as.character(npc_upreg10kb$V4))
npc_upreg25kb_res <- subset(npc_genes_res, rownames(npc_genes_res) %in% as.character(npc_upreg25kb$V4))
npc_upreg50kb_res <- subset(npc_genes_res, rownames(npc_genes_res) %in% as.character(npc_upreg50kb$V4))
# plotMA(npc_upreg10kb_res, cex=1.5, ylim=c(-4,4))
# plotMA(npc_upreg25kb_res, cex=1.5, ylim=c(-4,4))
# plotMA(npc_upreg50kb_res, cex=1.5, ylim=c(-4,4))
npc_nearby_genes_upregMMERVK10C <- data.frame(log2FC=npc_upreg10kb_res[,'log2FoldChange', drop=F], set=rep('10kb', nrow(npc_upreg10kb_res)))
npc_nearby_genes_upregMMERVK10C <- rbind(npc_nearby_genes_upregMMERVK10C, data.frame(log2FC=npc_upreg25kb_res[,'log2FoldChange', drop=F], set=rep('25kb', nrow(npc_upreg25kb_res))))
npc_nearby_genes_upregMMERVK10C <- rbind(npc_nearby_genes_upregMMERVK10C, data.frame(log2FC=npc_upreg50kb_res[,'log2FoldChange', drop=F], set=rep('50kb', nrow(npc_upreg50kb_res))))
npc_nearby_genes_upregMMERVK10C$disregulation <- rep("Upregulated", nrow(npc_nearby_genes_upregMMERVK10C))

# Read intersection of NOT upregulated MMERVK10C elements and genes in 10kb, 25kb and 50kb distance
npc_notupreg10kb <- fread('/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/not_upregulated/npc/npc_notupreg_MMERVK10C_10kb_windows_intersect_genes.bed', data.table = F, header = F)
npc_notupreg25kb <- fread('/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/not_upregulated/npc/npc_notupreg_MMERVK10C_25kb_windows_intersect_genes.bed', data.table = F, header = F)
npc_notupreg50kb <- fread('/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/not_upregulated/npc/npc_notupreg_MMERVK10C_50kb_windows_intersect_genes.bed', data.table = F, header = F)

# Are these genes dysregulated?
npc_notupreg10kb_res <- subset(npc_genes_res, rownames(npc_genes_res) %in% as.character(npc_notupreg10kb$V4))
npc_notupreg25kb_res <- subset(npc_genes_res, rownames(npc_genes_res) %in% as.character(npc_notupreg25kb$V4))
npc_notupreg50kb_res <- subset(npc_genes_res, rownames(npc_genes_res) %in% as.character(npc_notupreg50kb$V4))
# plotMA(npc_notupreg10kb_res, cex=1.5, ylim=c(-4,4))
# plotMA(npc_notupreg25kb_res, cex=1.5, ylim=c(-4,4))
# plotMA(npc_notupreg50kb_res, cex=1.5, ylim=c(-4,4))
npc_nearby_genes_notupregMMERVK10C <- data.frame(log2FC=npc_notupreg10kb_res[,'log2FoldChange', drop=F], set=rep('10kb', nrow(npc_notupreg10kb_res)))
npc_nearby_genes_notupregMMERVK10C <- rbind(npc_nearby_genes_notupregMMERVK10C, data.frame(log2FC=npc_notupreg25kb_res[,'log2FoldChange', drop=F], set=rep('25kb', nrow(npc_notupreg25kb_res))))
npc_nearby_genes_notupregMMERVK10C <- rbind(npc_nearby_genes_notupregMMERVK10C, data.frame(log2FC=npc_notupreg50kb_res[,'log2FoldChange', drop=F], set=rep('50kb', nrow(npc_notupreg50kb_res))))
npc_nearby_genes_notupregMMERVK10C$disregulation <- rep("Not Upregulated", nrow(npc_nearby_genes_notupregMMERVK10C))

npc_nearby_genes_MMERVK10C <- rbind(npc_nearby_genes_upregMMERVK10C, npc_nearby_genes_notupregMMERVK10C)

npc_nearby_genes_MMERVK10C_plot <- ggplot(npc_nearby_genes_MMERVK10C, aes(y=log2FoldChange, x=set, fill=disregulation)) + 
  geom_boxplot() + theme_classic() + scale_fill_brewer(palette="Set2") + 
  labs(fill = "Distance to\nfull MMERVK10C", x="") + ggtitle("NPC - Genes nearby upregulated MMERVK10Cs")+ 
  ylim(c(-3,6)) + stat_compare_means(method ='t.test')

ggsave(npc_nearby_genes_MMERVK10C_plot, file='/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/npc/npc_nearby_gene_MMERVK10C.svg', width=20, height=20, units="cm", dpi=96)





