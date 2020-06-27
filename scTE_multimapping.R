# Settings ----
library(data.table)
library(ggplot2)
library(stringr)
# Reading and coldata ----
scTE_default <- fread('/Volumes/Seagate Backup /scTE/May2020/2_stdmapping/2_counts/1_featureCounts/0_gene/gene_count_matrix_1.csv', data.table=F)
colnames(scTE_default)[7:ncol(scTE_default)] <- sapply(str_split(colnames(scTE_default)[7:ncol(scTE_default)], '/'), `[[`, 9)

scTE_multi <- fread('/Volumes/Seagate Backup /scTE/May2020/3_multimapping/2_counts/0_TEtranscripts/scTE_trim28.cntTable', data.table=F)
colnames(scTE_multi)[-1] <- sapply(str_split(colnames(scTE_multi)[-1], '/'), `[[`, 9)

coldata <- data.frame(sample=colnames(scTE_multi)[-1],
                      cell_type=sapply(str_split(colnames(scTE_multi)[-1], '_'), `[[`, 3),
                      condition=sapply(str_split(colnames(scTE_multi)[-1], '_'), `[[`, 2))
coldata$sample <- as.character(coldata$sample)
coldata$cell_type <- as.character(coldata$cell_type)
coldata$condition <- as.character(coldata$condition)
rownames(coldata) <- coldata$sample

rownames(scTE_multi) <- scTE_multi$`gene/TE`
rownames(scTE_default) <- scTE_default$Geneid

# Normalization by gene expression default mapping ----
scTE_default_EX <- scTE_default[,subset(coldata, coldata$cell_type == 'EX')$sample]
scTE_default_IN <- scTE_default[,subset(coldata, coldata$cell_type == 'IN')$sample]
scTE_default_Astro <- scTE_default[,subset(coldata, coldata$cell_type == 'Astro')$sample]
scTE_default_Oligo <- scTE_default[,subset(coldata, coldata$cell_type == 'Oligo')$sample]
scTE_default_OligoPre <- scTE_default[,subset(coldata, coldata$cell_type == 'OligoPre')$sample]
scTE_default_Vas <- scTE_default[,subset(coldata, coldata$cell_type == 'Vas')$sample]
scTE_default_MGlia <- scTE_default[,subset(coldata, coldata$cell_type == 'MGlia')$sample]

coldata_EX <- subset(coldata, coldata$cell_type == 'EX')
coldata_IN <- subset(coldata, coldata$cell_type == 'IN')
coldata_Astro <- subset(coldata, coldata$cell_type == 'Astro')
coldata_Oligo <- subset(coldata, coldata$cell_type == 'Oligo')
coldata_OligoPre <- subset(coldata, coldata$cell_type == 'OligoPre')
coldata_Vas <- subset(coldata, coldata$cell_type == 'Vas')
coldata_MGlia <- subset(coldata, coldata$cell_type == 'MGlia')

library(DESeq2)
scTE_default_EX_dds <- DESeqDataSetFromMatrix(scTE_default_EX[,rownames(coldata_EX)], coldata_EX, design = ~ condition)
scTE_default_EX_dds <- DESeq(scTE_default_EX_dds)

scTE_default_IN_dds <- DESeqDataSetFromMatrix(scTE_default_IN[,rownames(coldata_IN)], coldata_IN, design = ~ condition)
scTE_default_IN_dds <- DESeq(scTE_default_IN_dds)

scTE_default_MGlia_dds <- DESeqDataSetFromMatrix(scTE_default_MGlia[,rownames(coldata_MGlia)], coldata_MGlia, design = ~ condition)
scTE_default_MGlia_dds <- DESeq(scTE_default_MGlia_dds)

scTE_default_Vas_dds <- DESeqDataSetFromMatrix(scTE_default_Vas[,rownames(coldata_Vas)], coldata_Vas, design = ~ condition)
scTE_default_Vas_dds <- DESeq(scTE_default_Vas_dds)

scTE_default_Astro_dds <- DESeqDataSetFromMatrix(scTE_default_Astro[,rownames(coldata_Astro)], coldata_Astro, design = ~ condition)
scTE_default_Astro_dds <- DESeq(scTE_default_Astro_dds)

scTE_default_Oligo_dds <- DESeqDataSetFromMatrix(scTE_default_Oligo[,rownames(coldata_Oligo)], coldata_Oligo, design = ~ condition)
scTE_default_Oligo_dds <- DESeq(scTE_default_Oligo_dds)

scTE_default_OligoPre_dds <- DESeqDataSetFromMatrix(scTE_default_OligoPre[,rownames(coldata_OligoPre)], coldata_OligoPre, design = ~ condition)
scTE_default_OligoPre_dds <- DESeq(scTE_default_OligoPre_dds)

# DEA on TE expression : EX ----
scTE_multi_EX <- scTE_multi[which(!startsWith(rownames(scTE_multi), 'ENS')),subset(coldata, coldata$cell_type == 'EX')$sample]
scTE_multi_EX <- scTE_multi_EX[which(rowSums(scTE_multi_EX[,coldata_EX$sample]) > 0),]
scTE_multi_EX$TE_subfamily <- sapply(str_split(rownames(scTE_multi_EX), ":"), `[[`, 1)
scTE_multi_EX$TE_family <- sapply(str_split(rownames(scTE_multi_EX), ":"), `[[`, 2)
scTE_multi_EX$TE_class <- sapply(str_split(rownames(scTE_multi_EX), ":"), `[[`, 3)
scTE_multi_EX <- scTE_multi_EX[which(scTE_multi_EX$TE_class %in% c('LINE', 'SINE', 'LTR', 'Retroposon')),]
scTE_multi_EX <- scTE_multi_EX[,coldata_EX$sample]

scTE_multi_EX_dds <- DESeqDataSetFromMatrix(scTE_multi_EX[,rownames(coldata_EX)], coldata_EX, design = ~ condition)
scTE_multi_EX_dds <- DESeq(scTE_multi_EX_dds)
scTE_multi_EX_res <- results(scTE_multi_EX_dds)
scTE_multi_EX_res_df <- as.data.frame(scTE_multi_EX_res)
scTE_multi_EX_res_df$TE_id <- rownames(scTE_multi_EX_res_df)

svg('/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/MAplots/EX_maplot.svg', height = 5, width = 5)
plotMA(scTE_multi_EX_dds, ylim=c(-4,4), alpha=0.05)
abline(h=1, lty=2, col='mistyrose3')
abline(h=-1, lty=2, col='mistyrose3')
dev.off()

png('/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/MAplots/EX_maplot.png', height = 1300, width = 1300, units = 'px', res=320)
plotMA(scTE_multi_EX_dds, ylim=c(-4,4), alpha=0.05)
abline(h=1, lty=2, col='mistyrose3')
abline(h=-1, lty=2, col='mistyrose3')
dev.off()

# IN ----
scTE_multi_IN <- scTE_multi[which(!startsWith(rownames(scTE_multi), 'ENS')),subset(coldata, coldata$cell_type == 'IN')$sample]
scTE_multi_IN <- scTE_multi_IN[which(rowSums(scTE_multi_IN[,coldata_IN$sample]) > 0),]
scTE_multi_IN$TE_subfamily <- sapply(str_split(rownames(scTE_multi_IN), ":"), `[[`, 1)
scTE_multi_IN$TE_family <- sapply(str_split(rownames(scTE_multi_IN), ":"), `[[`, 2)
scTE_multi_IN$TE_class <- sapply(str_split(rownames(scTE_multi_IN), ":"), `[[`, 3)
scTE_multi_IN <- scTE_multi_IN[which(scTE_multi_IN$TE_class %in% c('LINE', 'SINE', 'LTR', 'Retroposon')),]
scTE_multi_IN <- scTE_multi_IN[,coldata_IN$sample]

scTE_multi_IN_dds <- DESeqDataSetFromMatrix(scTE_multi_IN[,rownames(coldata_IN)], coldata_IN, design = ~ condition)
scTE_multi_IN_dds <- DESeq(scTE_multi_IN_dds)
scTE_multi_IN_res <- results(scTE_multi_IN_dds)
scTE_multi_IN_res_df <- as.data.frame(scTE_multi_IN_res)
scTE_multi_IN_res_df$TE_id <- rownames(scTE_multi_IN_res_df)

svg('/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/MAplots/IN_maplot.svg', height = 1300, width = 1300, units = 'px', res=320)
plotMA(scTE_multi_IN_dds, ylim=c(-4,4), alpha=0.05)
abline(h=1, lty=2, col='mistyrose3')
abline(h=-1, lty=2, col='mistyrose3')
dev.off()

png('/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/MAplots/IN_maplot.png', height = 1300, width = 1300, units = 'px', res=320)
plotMA(scTE_multi_IN_dds, ylim=c(-4,4), alpha=0.05)
abline(h=1, lty=2, col='mistyrose3')
abline(h=-1, lty=2, col='mistyrose3')
dev.off()

# Oligo ----
scTE_multi_Oligo <- scTE_multi[which(!startsWith(rownames(scTE_multi), 'ENS')),subset(coldata, coldata$cell_type == 'Oligo')$sample]
scTE_multi_Oligo <- scTE_multi_Oligo[which(rowSums(scTE_multi_Oligo[,coldata_Oligo$sample]) > 0),]
scTE_multi_Oligo$TE_subfamily <- sapply(str_split(rownames(scTE_multi_Oligo), ":"), `[[`, 1)
scTE_multi_Oligo$TE_family <- sapply(str_split(rownames(scTE_multi_Oligo), ":"), `[[`, 2)
scTE_multi_Oligo$TE_class <- sapply(str_split(rownames(scTE_multi_Oligo), ":"), `[[`, 3)
scTE_multi_Oligo <- scTE_multi_Oligo[which(scTE_multi_Oligo$TE_class %in% c('LINE', 'SINE', 'LTR', 'Retroposon')),]
scTE_multi_Oligo <- scTE_multi_Oligo[,coldata_Oligo$sample]

scTE_multi_Oligo_dds <- DESeqDataSetFromMatrix(scTE_multi_Oligo[,rownames(coldata_Oligo)], coldata_Oligo, design = ~ condition)
scTE_multi_Oligo_dds <- DESeq(scTE_multi_Oligo_dds)
scTE_multi_Oligo_res <- results(scTE_multi_Oligo_dds)
scTE_multi_Oligo_res_df <- as.data.frame(scTE_multi_Oligo_res)
scTE_multi_Oligo_res_df$TE_id <- rownames(scTE_multi_Oligo_res_df)

svg('/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/MAplots/Oligo_maplot.svg', height = 5, width = 5)
plotMA(scTE_multi_Oligo_dds, ylim=c(-4,4), alpha=0.05)
abline(h=1, lty=2, col='mistyrose3')
abline(h=-1, lty=2, col='mistyrose3')
dev.off()

png('/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/MAplots/Oligo_maplot.png', height = 480, width = 480, units = 'px', res=200)
plotMA(scTE_multi_Oligo_dds, ylim=c(-4,4), alpha=0.05)
abline(h=1, lty=2, col='mistyrose3')
abline(h=-1, lty=2, col='mistyrose3')
dev.off()

# Astro ----
scTE_multi_Astro <- scTE_multi[which(!startsWith(rownames(scTE_multi), 'ENS')),subset(coldata, coldata$cell_type == 'Astro')$sample]
scTE_multi_Astro <- scTE_multi_Astro[which(rowSums(scTE_multi_Astro[,coldata_Astro$sample]) > 0),]
scTE_multi_Astro$TE_subfamily <- sapply(str_split(rownames(scTE_multi_Astro), ":"), `[[`, 1)
scTE_multi_Astro$TE_family <- sapply(str_split(rownames(scTE_multi_Astro), ":"), `[[`, 2)
scTE_multi_Astro$TE_class <- sapply(str_split(rownames(scTE_multi_Astro), ":"), `[[`, 3)
scTE_multi_Astro <- scTE_multi_Astro[which(scTE_multi_Astro$TE_class %in% c('LINE', 'SINE', 'LTR', 'Retroposon')),]
scTE_multi_Astro <- scTE_multi_Astro[,coldata_Astro$sample]

scTE_multi_Astro_dds <- DESeqDataSetFromMatrix(scTE_multi_Astro[,rownames(coldata_Astro)], coldata_Astro, design = ~ condition)
scTE_multi_Astro_dds <- DESeq(scTE_multi_Astro_dds)
scTE_multi_Astro_res <- results(scTE_multi_Astro_dds)
scTE_multi_Astro_res_df <- as.data.frame(scTE_multi_Astro_res)
scTE_multi_Astro_res_df$TE_id <- rownames(scTE_multi_Astro_res_df)

svg('/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/MAplots/Astro_maplot.svg', height = 5, width = 5)
plotMA(scTE_multi_Astro_dds, ylim=c(-4,4), alpha=0.05)
abline(h=1, lty=2, col='mistyrose3')
abline(h=-1, lty=2, col='mistyrose3')
dev.off()

png('/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/MAplots/Astro_maplot.png', height = 480, width = 480, units = 'px', res=200)
plotMA(scTE_multi_Astro_dds, ylim=c(-4,4), alpha=0.05)
abline(h=1, lty=2, col='mistyrose3')
abline(h=-1, lty=2, col='mistyrose3')
dev.off()

# Vas ----
scTE_multi_Vas <- scTE_multi[which(!startsWith(rownames(scTE_multi), 'ENS')),subset(coldata, coldata$cell_type == 'Vas')$sample]
scTE_multi_Vas <- scTE_multi_Vas[which(rowSums(scTE_multi_Vas[,coldata_Vas$sample]) > 0),]
scTE_multi_Vas$TE_subfamily <- sapply(str_split(rownames(scTE_multi_Vas), ":"), `[[`, 1)
scTE_multi_Vas$TE_family <- sapply(str_split(rownames(scTE_multi_Vas), ":"), `[[`, 2)
scTE_multi_Vas$TE_class <- sapply(str_split(rownames(scTE_multi_Vas), ":"), `[[`, 3)
scTE_multi_Vas <- scTE_multi_Vas[which(scTE_multi_Vas$TE_class %in% c('LINE', 'SINE', 'LTR', 'Retroposon')),]
scTE_multi_Vas <- scTE_multi_Vas[,coldata_Vas$sample]

scTE_multi_Vas_dds <- DESeqDataSetFromMatrix(scTE_multi_Vas[,rownames(coldata_Vas)], coldata_Vas, design = ~ condition)
scTE_multi_Vas_dds <- DESeq(scTE_multi_Vas_dds)
scTE_multi_Vas_res <- results(scTE_multi_Vas_dds)
scTE_multi_Vas_res_df <- as.data.frame(scTE_multi_Vas_res)
scTE_multi_Vas_res_df$TE_id <- rownames(scTE_multi_Vas_res_df)

svg('/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/MAplots/Vas_maplot.svg', height = 5, width = 5)
plotMA(scTE_multi_Vas_dds, ylim=c(-4,4), alpha=0.05)
abline(h=1, lty=2, col='mistyrose3')
abline(h=-1, lty=2, col='mistyrose3')
dev.off()

png('/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/MAplots/Vas_maplot.png', height = 480, width = 480, units = 'px', res=200)
plotMA(scTE_multi_Vas_dds, ylim=c(-4,4), alpha=0.05)
abline(h=1, lty=2, col='mistyrose3')
abline(h=-1, lty=2, col='mistyrose3')
dev.off()

# MGlia ----
scTE_multi_MGlia <- scTE_multi[which(!startsWith(rownames(scTE_multi), 'ENS')),subset(coldata, coldata$cell_type == 'MGlia')$sample]
scTE_multi_MGlia <- scTE_multi_MGlia[which(rowSums(scTE_multi_MGlia[,coldata_MGlia$sample]) > 0),]
scTE_multi_MGlia$TE_subfamily <- sapply(str_split(rownames(scTE_multi_MGlia), ":"), `[[`, 1)
scTE_multi_MGlia$TE_family <- sapply(str_split(rownames(scTE_multi_MGlia), ":"), `[[`, 2)
scTE_multi_MGlia$TE_class <- sapply(str_split(rownames(scTE_multi_MGlia), ":"), `[[`, 3)
scTE_multi_MGlia <- scTE_multi_MGlia[which(scTE_multi_MGlia$TE_class %in% c('LINE', 'SINE', 'LTR', 'Retroposon')),]
scTE_multi_MGlia <- scTE_multi_MGlia[,coldata_MGlia$sample]

scTE_multi_MGlia_dds <- DESeqDataSetFromMatrix(scTE_multi_MGlia[,rownames(coldata_MGlia)], coldata_MGlia, design = ~ condition)
scTE_multi_MGlia_dds <- DESeq(scTE_multi_MGlia_dds)
scTE_multi_MGlia_res <- results(scTE_multi_MGlia_dds)
scTE_multi_MGlia_res_df <- as.data.frame(scTE_multi_MGlia_res)
scTE_multi_MGlia_res_df$TE_id <- rownames(scTE_multi_MGlia_res_df)

svg('/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/MAplots/MGlia_maplot.svg', height = 5, width = 5)
plotMA(scTE_multi_MGlia_dds, ylim=c(-4,4), alpha=0.05)
abline(h=1, lty=2, col='mistyrose3')
abline(h=-1, lty=2, col='mistyrose3')
dev.off()

png('/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/MAplots/MGlia_maplot.png', height = 480, width = 480, units = 'px', res=200)
plotMA(scTE_multi_MGlia_dds, ylim=c(-4,4), alpha=0.05)
abline(h=1, lty=2, col='mistyrose3')
abline(h=-1, lty=2, col='mistyrose3')
dev.off()


# OligoPre ----
scTE_multi_OligoPre <- scTE_multi[which(!startsWith(rownames(scTE_multi), 'ENS')),subset(coldata, coldata$cell_type == 'OligoPre')$sample]
scTE_multi_OligoPre <- scTE_multi_OligoPre[which(rowSums(scTE_multi_OligoPre[,coldata_OligoPre$sample]) > 0),]
scTE_multi_OligoPre$TE_subfamily <- sapply(str_split(rownames(scTE_multi_OligoPre), ":"), `[[`, 1)
scTE_multi_OligoPre$TE_family <- sapply(str_split(rownames(scTE_multi_OligoPre), ":"), `[[`, 2)
scTE_multi_OligoPre$TE_class <- sapply(str_split(rownames(scTE_multi_OligoPre), ":"), `[[`, 3)
scTE_multi_OligoPre <- scTE_multi_OligoPre[which(scTE_multi_OligoPre$TE_class %in% c('LINE', 'SINE', 'LTR', 'Retroposon')),]
scTE_multi_OligoPre <- scTE_multi_OligoPre[,coldata_OligoPre$sample]

scTE_multi_OligoPre_dds <- DESeqDataSetFromMatrix(scTE_multi_OligoPre[,rownames(coldata_OligoPre)], coldata_OligoPre, design = ~ condition)
scTE_multi_OligoPre_dds <- DESeq(scTE_multi_OligoPre_dds)
scTE_multi_OligoPre_res <- results(scTE_multi_OligoPre_dds)
scTE_multi_OligoPre_res_df <- as.data.frame(scTE_multi_OligoPre_res)
scTE_multi_OligoPre_res_df$TE_id <- rownames(scTE_multi_OligoPre_res_df)

svg('/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/MAplots/OligoPre_maplot.svg', height = 5, width = 5)
plotMA(scTE_multi_OligoPre_dds, ylim=c(-4,4), alpha=0.05)
abline(h=1, lty=2, col='mistyrose3')
abline(h=-1, lty=2, col='mistyrose3')
dev.off()

png('/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/MAplots/OligoPre_maplot.png', height = 480, width = 480, units = 'px', res=200)
plotMA(scTE_multi_OligoPre_dds, ylim=c(-4,4), alpha=0.05)
abline(h=1, lty=2, col='mistyrose3')
abline(h=-1, lty=2, col='mistyrose3')
dev.off()

scTE_multi_EX_res_df[is.na(scTE_multi_EX_res_df$padj),'padj'] <- 1
scTE_multi_EX_res_df$significant <- as.factor(ifelse(as.numeric(as.character(scTE_multi_EX_res_df$padj)) < 0.05, 1, 0))
scTE_multi_OligoPre_res_df$significant <- as.factor(ifelse(as.numeric(as.character(scTE_multi_OligoPre_res_df$padj)) < 0.05, 1, 0))
scTE_multi_IN_res_df$significant <- as.factor(ifelse(as.numeric(as.character(scTE_multi_IN_res_df$padj)) < 0.05, 1, 0))
scTE_multi_Oligo_res_df$significant <- as.factor(ifelse(as.numeric(as.character(scTE_multi_Oligo_res_df$padj)) < 0.05, 1, 0))
scTE_multi_Astro_res_df$significant <- as.factor(ifelse(as.numeric(as.character(scTE_multi_Astro_res_df$padj)) < 0.05, 1, 0))
scTE_multi_Vas_res_df$significant <- as.factor(ifelse(as.numeric(as.character(scTE_multi_Vas_res_df$padj)) < 0.05, 1, 0))
scTE_multi_MGlia_res_df$significant <- as.factor(ifelse(as.numeric(as.character(scTE_multi_MGlia_res_df$padj)) < 0.05, 1, 0))

scTE_multi_OligoPre_res_df$TE_id <- sapply(str_split(scTE_multi_OligoPre_res_df$TE_id, ':'), `[[`, 1)
scTE_multi_OligoPre_res_df$TE_id <- ifelse(scTE_multi_OligoPre_res_df$significant == 1, scTE_multi_OligoPre_res_df$TE_id, "")
scTE_multi_EX_res_df$TE_id <- sapply(str_split(scTE_multi_EX_res_df$TE_id, ':'), `[[`, 1)
scTE_multi_EX_res_df$TE_id <- ifelse(scTE_multi_EX_res_df$significant == 1, scTE_multi_EX_res_df$TE_id, "")
scTE_multi_IN_res_df$TE_id <- sapply(str_split(scTE_multi_IN_res_df$TE_id, ':'), `[[`, 1)
scTE_multi_IN_res_df$TE_id <- ifelse(scTE_multi_IN_res_df$significant == 1, scTE_multi_IN_res_df$TE_id, "")
scTE_multi_Oligo_res_df$TE_id <- sapply(str_split(scTE_multi_Oligo_res_df$TE_id, ':'), `[[`, 1)
scTE_multi_Oligo_res_df$TE_id <- ifelse(scTE_multi_Oligo_res_df$significant == 1, scTE_multi_Oligo_res_df$TE_id, "")
scTE_multi_Astro_res_df$TE_id <- sapply(str_split(scTE_multi_Astro_res_df$TE_id, ':'), `[[`, 1)
scTE_multi_Astro_res_df$TE_id <- ifelse(scTE_multi_Astro_res_df$significant == 1, scTE_multi_Astro_res_df$TE_id, "")
scTE_multi_Vas_res_df$TE_id <- sapply(str_split(scTE_multi_Vas_res_df$TE_id, ':'), `[[`, 1)
scTE_multi_Vas_res_df$TE_id <- ifelse(scTE_multi_Vas_res_df$significant == 1, scTE_multi_Vas_res_df$TE_id, "")
scTE_multi_MGlia_res_df$TE_id <- sapply(str_split(scTE_multi_MGlia_res_df$TE_id, ':'), `[[`, 1)
scTE_multi_MGlia_res_df$TE_id <- ifelse(scTE_multi_MGlia_res_df$significant == 1, scTE_multi_MGlia_res_df$TE_id, "")

scTE_multi_Astro_dds


library(ggrepel)
oligoprep <- ggplot(scTE_multi_OligoPre_res_df, aes(x=log2(baseMean), y=log2FoldChange, colour=significant, labels=TE_id)) + 
  geom_point() + theme_classic()+ theme(legend.position = "none") + scale_colour_manual(values=c("0"="black", "1"="firebrick2")) +
  geom_label_repel(aes(label = TE_id),box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50')

ggsave('/Volumes/Seagate Backup /scTE/April2020/3_multimapping/plots/MAplots/OligoPre_MAplot.svg', oligoprep)

exp <- ggplot(scTE_multi_EX_res_df, aes(x=log2(baseMean), y=log2FoldChange, colour=significant, labels=TE_id)) + 
  geom_point() + theme_classic()+ theme(legend.position = "none") + scale_colour_manual(values=c("0"="black", "1"="firebrick2")) +
  geom_label_repel(aes(label = TE_id),box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50')

ggsave('/Volumes/Seagate Backup /scTE/April2020/3_multimapping/plots/MAplots/EX_MAplot.svg', exp)

inp <- ggplot(scTE_multi_IN_res_df, aes(x=log2(baseMean), y=log2FoldChange, colour=significant, labels=TE_id)) + 
  geom_point() + theme_classic()+ theme(legend.position = "none") + scale_colour_manual(values=c("0"="black", "1"="firebrick2")) +
  geom_label_repel(aes(label = TE_id),box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50')

ggsave('/Volumes/Seagate Backup /scTE/April2020/3_multimapping/plots/MAplots/IN_MAplot.svg', inp)

oligop <- ggplot(scTE_multi_Oligo_res_df, aes(x=log2(baseMean), y=log2FoldChange, colour=significant, labels=TE_id)) + 
  geom_point() + theme_classic()+ theme(legend.position = "none") + scale_colour_manual(values=c("0"="black", "1"="firebrick2")) + 
  geom_label_repel(aes(label = TE_id),box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50')

ggsave('/Volumes/Seagate Backup /scTE/April2020/3_multimapping/plots/MAplots/Oligo_MAplot.svg', oligop)

astrop <- ggplot(scTE_multi_Astro_res_df, aes(x=log2(baseMean), y=log2FoldChange, colour=significant, labels=TE_id)) + 
  geom_point() + theme_classic()+ theme(legend.position = "none") + scale_colour_manual(values=c("0"="black", "1"="firebrick2")) + 
  geom_label_repel(aes(label = TE_id),box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50')

ggsave('/Volumes/Seagate Backup /scTE/April2020/3_multimapping/plots/MAplots/Astro_MAplot.svg', astrop)

vasp <- ggplot(scTE_multi_Vas_res_df, aes(x=log2(baseMean), y=log2FoldChange, colour=significant, labels=TE_id)) + 
  geom_point() + theme_classic()+ theme(legend.position = "none") + scale_colour_manual(values=c("0"="black", "1"="firebrick2")) + 
  geom_label_repel(aes(label = TE_id),box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50')

ggsave('/Volumes/Seagate Backup /scTE/April2020/3_multimapping/plots/MAplots/Vas_MAplot.svg', vasp)

mgliap <- ggplot(scTE_multi_MGlia_res_df, aes(x=log2(baseMean), y=log2FoldChange, colour=significant, labels=TE_id)) + 
  geom_point() + theme_classic()+ theme(legend.position = "none") + scale_colour_manual(values=c("0"="black", "1"="firebrick2")) + 
  geom_label_repel(aes(label = TE_id),box.padding   = 0.35, point.padding = 0.5, segment.color = 'grey50')

ggsave('/Volumes/Seagate Backup /scTE/April2020/3_multimapping/plots/MAplots/MGlia_MAplot.svg', mgliap)

# Mean plots ----
produce_mean_plot <- function(counts, dds, res_df, coldata, celltype, filename){
  counts_norm <- counts
  counts_norm[] <- mapply('/', counts_norm[,names(dds$sizeFactor)], dds$sizeFactor)
  res_df$FoldChange <- 2^(res_df$log2FoldChange)
  counts_norm$type <- ifelse(rownames(counts_norm) %in% rownames(res_df[which(res_df$padj < 0.01 & res_df$log2FoldChange > 3),]), 'Upregulated', ifelse(rownames(counts_norm) %in% rownames(res_df[which(res_df$padj < 0.01 & res_df$log2FoldChange < -3),]), 'Downregulated', 'Not significant'))
  counts_norm$colours <- ifelse(counts_norm$type == 'Upregulated', 'firebrick2', ifelse(counts_norm$type == 'Downregulated', 'steelblue', 'black'))
  counts_norm$cexs <- ifelse(counts_norm$type != 'Not significant', 1, 0.6)
  
  signdiff <- res_df[rownames(res_df[which(res_df$padj < 0.01 & res_df$log2FoldChange > 3),]), 'FoldChange', drop=F]
  signdiff$TE_id <- sapply(str_split(rownames(signdiff), ':'), `[[`, 1)
  signdiff$TE_id <- factor(signdiff$TE_id, levels=signdiff[order(signdiff$FoldChange),'TE_id'])
  
  if(nrow(signdiff) > 0){
    barplot <- ggplot(signdiff, aes(x=TE_id, y=FoldChange)) + geom_bar(stat='identity') + theme_classic() + coord_flip() + labs(x='') 
    svg(paste(filename, 'barplot.svg', sep=''))
    print(barplot)
    dev.off()
  }
  
  svg(paste(filename, 'meanplot.svg', sep=''))
  plot(log2(rowMeans(counts_norm[,subset(coldata, coldata$condition == 'Ctl')$sample])+0.5),
       log2(rowMeans(counts_norm[,subset(coldata, coldata$condition == 'KO')$sample])+0.5),
       xlab = 'log2(mean Control expression)',
       ylab = 'log2(mean KO expression)',
       pch=16, cex=counts_norm$cexs,
       col = counts_norm$colours, main=celltype, xlim=c(-1,20), ylim=c(-1,20))
  
  legend("bottomright", legend = c(paste("up (",as.numeric(table(counts_norm$type)["Upregulated"]),")",sep=""),
                                   paste("down (",as.numeric(table(counts_norm$type)["Downregulated"]),")",sep = ""),
                                   paste("not significant (",as.numeric(table(counts_norm$type)["Not significant"]),")",sep = "")),
         pch=16,col=c("firebrick3","steelblue","black"),cex=1)
  dev.off()
  # text(log2(rowMeans(counts_norm[,subset(coldata, coldata$condition == 'Ctl')$sample])+0.5),
  #      log2(rowMeans(counts_norm[,subset(coldata, coldata$condition == 'KO')$sample])+0.5),
  #      labels=counts_norm$labels, cex= 0.7, pos=3)
  # 
  # counts_norm$labels <- ifelse(counts_norm$type == 'Upregulated', rownames(counts_norm), "")
  return(signdiff)
}

EX_signdiff <- produce_mean_plot(counts=scTE_multi_EX, dds=scTE_default_EX_dds, res_df=scTE_multi_EX_res_df, coldata=coldata_EX, celltype = 'Excitatory neurons', filename='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/meanplots/EX_')
IN_signdiff <- produce_mean_plot(scTE_multi_IN, scTE_default_IN_dds, scTE_multi_IN_res_df, coldata_IN, 'Inhibitory neurons', '/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/meanplots/IN_')
Astro_signdiff <- produce_mean_plot(scTE_multi_Astro, scTE_default_Astro_dds, scTE_multi_Astro_res_df, coldata_Astro, 'Astrocytes', '/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/meanplots/Astro_')
Oligo_signdiff <- produce_mean_plot(counts=scTE_multi_Oligo, dds=scTE_default_Oligo_dds, res_df=scTE_multi_Oligo_res_df, coldata=coldata_Oligo, celltype = 'Oligodendrocytes', filename = '/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/meanplots/Oligo_')
OligoPre_signdiff <- produce_mean_plot(scTE_multi_OligoPre, scTE_default_OligoPre_dds, scTE_multi_OligoPre_res_df, coldata_OligoPre, 'Oligodendrocyte precursors', '/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/meanplots/OligoPre_')
MGlia_signdiff <- produce_mean_plot(scTE_multi_MGlia, scTE_default_MGlia_dds, scTE_multi_MGlia_res_df, coldata_MGlia, 'Microglia', '/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/meanplots/MGlia_')
Vas_signdiff <- produce_mean_plot(scTE_multi_Vas, scTE_default_Vas_dds, scTE_multi_Vas_res_df, coldata_Vas, 'Vascular', '/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/meanplots/Vas_')

# 
# EX_signdiff$type <- rep('EX',nrow(EX_signdiff))
# IN_signdiff$type <- rep('IN',nrow(IN_signdiff))
# Astro_signdiff$type <- rep('Astro',nrow(Astro_signdiff))
# Oligo_signdiff$type <- rep('Oligo',nrow(Oligo_signdiff))
# OligoPre_signdiff$type <- rep('OligoPre',nrow(OligoPre_signdiff))
# MGlia_signdiff$type <- rep('MGlia',nrow(MGlia_signdiff))
# Vas_signdiff$type <- rep('Vas',nrow(Vas_signdiff))
# 
# signdiff_allgroups <- unique(c(rownames(EX_signdiff), 
#        rownames(IN_signdiff),
#        rownames(Astro_signdiff),
#        rownames(Oligo_signdiff),
#        rownames(OligoPre_signdiff),
#        rownames(MGlia_signdiff),
#        rownames(Vas_signdiff)))
# 
# scTE_multi_EX_res_df_signdiff_allgroups <- scTE_multi_EX_res_df[signdiff_allgroups[which(signdiff_allgroups %in% rownames(scTE_multi_EX_res_df))], c('log2FoldChange'), drop=F]
# scTE_multi_IN_res_df_signdiff_allgroups <- scTE_multi_IN_res_df[signdiff_allgroups[which(signdiff_allgroups %in% rownames(scTE_multi_IN_res_df))], c('log2FoldChange'), drop=F]
# scTE_multi_Astro_res_df_signdiff_allgroups <- scTE_multi_Astro_res_df[signdiff_allgroups[which(signdiff_allgroups %in% rownames(scTE_multi_Astro_res_df))], c('log2FoldChange'), drop=F]
# scTE_multi_Oligo_res_df_signdiff_allgroups <- scTE_multi_Oligo_res_df[signdiff_allgroups[which(signdiff_allgroups %in% rownames(scTE_multi_Oligo_res_df))], c('log2FoldChange'), drop=F]
# scTE_multi_OligoPre_res_df_signdiff_allgroups <- scTE_multi_OligoPre_res_df[signdiff_allgroups[which(signdiff_allgroups %in% rownames(scTE_multi_OligoPre_res_df))], c('log2FoldChange'), drop=F]
# scTE_multi_MGlia_res_df_signdiff_allgroups <- scTE_multi_MGlia_res_df[signdiff_allgroups[which(signdiff_allgroups %in% rownames(scTE_multi_MGlia_res_df))], c('log2FoldChange'), drop=F]
# scTE_multi_Vas_res_df_signdiff_allgroups <- scTE_multi_Vas_res_df[signdiff_allgroups[which(signdiff_allgroups %in% rownames(scTE_multi_Vas_res_df))], c('log2FoldChange'), drop=F]
# 
# scTE_multi_EX_res_df_signdiff_allgroups$TE_id <- sapply(str_split(rownames(scTE_multi_EX_res_df_signdiff_allgroups), ':'), `[[`, 1)
# scTE_multi_IN_res_df_signdiff_allgroups$TE_id <- sapply(str_split(rownames(scTE_multi_IN_res_df_signdiff_allgroups), ':'), `[[`, 1)
# scTE_multi_Astro_res_df_signdiff_allgroups$TE_id <- sapply(str_split(rownames(scTE_multi_Astro_res_df_signdiff_allgroups), ':'), `[[`, 1)
# scTE_multi_Oligo_res_df_signdiff_allgroups$TE_id <- sapply(str_split(rownames(scTE_multi_Oligo_res_df_signdiff_allgroups), ':'), `[[`, 1)
# scTE_multi_OligoPre_res_df_signdiff_allgroups$TE_id <- sapply(str_split(rownames(scTE_multi_OligoPre_res_df_signdiff_allgroups), ':'), `[[`, 1)
# scTE_multi_MGlia_res_df_signdiff_allgroups$TE_id <- sapply(str_split(rownames(scTE_multi_MGlia_res_df_signdiff_allgroups), ':'), `[[`, 1)
# scTE_multi_Vas_res_df_signdiff_allgroups$TE_id <- sapply(str_split(rownames(scTE_multi_Vas_res_df_signdiff_allgroups), ':'), `[[`, 1)
# 
# scTE_multi_EX_res_df_signdiff_allgroups$cell_type <- rep("EX", nrow(scTE_multi_EX_res_df_signdiff_allgroups))
# scTE_multi_IN_res_df_signdiff_allgroups$cell_type <- rep("IN", nrow(scTE_multi_IN_res_df_signdiff_allgroups))
# scTE_multi_Astro_res_df_signdiff_allgroups$cell_type <- rep("Astro", nrow(scTE_multi_Astro_res_df_signdiff_allgroups))
# scTE_multi_Oligo_res_df_signdiff_allgroups$cell_type <- rep("Oligo", nrow(scTE_multi_Oligo_res_df_signdiff_allgroups))
# scTE_multi_OligoPre_res_df_signdiff_allgroups$cell_type <- rep("OligoPre", nrow(scTE_multi_OligoPre_res_df_signdiff_allgroups))
# scTE_multi_MGlia_res_df_signdiff_allgroups$cell_type <- rep("MGlia", nrow(scTE_multi_MGlia_res_df_signdiff_allgroups))
# scTE_multi_Vas_res_df_signdiff_allgroups$cell_type <- rep("Vas", nrow(scTE_multi_Vas_res_df_signdiff_allgroups))
# 
# signdiff_allgroups <- rbind(scTE_multi_EX_res_df_signdiff_allgroups,
#       scTE_multi_IN_res_df_signdiff_allgroups,
#       scTE_multi_Astro_res_df_signdiff_allgroups,
#       scTE_multi_Oligo_res_df_signdiff_allgroups,
#       scTE_multi_OligoPre_res_df_signdiff_allgroups,
#       scTE_multi_MGlia_res_df_signdiff_allgroups,
#       scTE_multi_Vas_res_df_signdiff_allgroups)
# 
# ggplot(data=signdiff_allgroups, aes(x=cell_type, y=log2FoldChange, fill=cell_type)) + 
#   theme_classic() +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank()) +
#   geom_bar(stat = 'identity') + 
#   facet_grid(TE_id ~ ., scales = 'free') 

# Volcano plots ---- 
library(EnhancedVolcano)
ex_volcano <- EnhancedVolcano(scTE_multi_EX_res,
                              lab = rownames(scTE_multi_EX_res),
                              x = 'log2FoldChange',
                              y = 'pvalue', xlim=c(-6,10), ylim=c(0,170), title = "Excitatory neurons", subtitle = "", subtitleLabSize = 0)

ggsave(plot=ex_volcano, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/volcanos/ex_volcano.svg', height = 7, width = 7)
ggsave(plot=ex_volcano, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/volcanos/ex_volcano.png', height = 7, width = 7)

oligo_volcano <- EnhancedVolcano(scTE_multi_Oligo_res,
                                 lab = rownames(scTE_multi_Oligo_res),
                                 x = 'log2FoldChange',
                                 y = 'pvalue', xlim=c(-6,10), ylim=c(0,170), title = "Oligodendrocytes", subtitle = "", subtitleLabSize = 0)+ theme(legend.position = "none")
ggsave(plot=oligo_volcano, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/volcanos/oligo_volcano.svg', height = 7, width = 7)
ggsave(plot=oligo_volcano, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/volcanos/oligo_volcano.png', height = 7, width = 7)

astro_volcano <- EnhancedVolcano(scTE_multi_Astro_res,
                                 lab = rownames(scTE_multi_Astro_res),
                                 x = 'log2FoldChange',
                                 y = 'pvalue', xlim=c(-6,10), ylim=c(0,170), title = "Astrocytes", subtitle = "", subtitleLabSize = 0)+ theme(legend.position = "none")
ggsave(plot=astro_volcano, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/volcanos/astro_volcano.svg', height = 7, width = 7)
ggsave(plot=astro_volcano, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/volcanos/astro_volcano.png', height = 7, width = 7)

oligopre_volcano <- EnhancedVolcano(scTE_multi_OligoPre_res,
                                    lab = rownames(scTE_multi_OligoPre_res),
                                    x = 'log2FoldChange',
                                    y = 'pvalue', xlim=c(-6,10), ylim=c(0,170), title = "Oligoprecursors", subtitle = "", subtitleLabSize = 0)+ theme(legend.position = "none")
ggsave(plot=oligopre_volcano, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/volcanos/oligopre_volcano.svg', height = 7, width = 7)
ggsave(plot=oligopre_volcano, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/volcanos/oligopre_volcano.png', height = 7, width = 7)

in_volcano <- EnhancedVolcano(scTE_multi_IN_res,
                              lab = rownames(scTE_multi_IN_res),
                              x = 'log2FoldChange',
                              y = 'pvalue', xlim=c(-6,10), ylim=c(0,170), title = "Interneurons", subtitle = "", subtitleLabSize = 0)+ theme(legend.position = "none")
ggsave(plot=in_volcano, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/volcanos/in_volcano.svg', height = 7, width = 7)
ggsave(plot=in_volcano, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/volcanos/in_volcano.png', height = 7, width = 7)

mglia_volcano <- EnhancedVolcano(scTE_multi_MGlia_res,
                                 lab = rownames(scTE_multi_MGlia_res),
                                 x = 'log2FoldChange',
                                 y = 'pvalue', xlim=c(-6,10), ylim=c(0,170), title = "Microglia", subtitle = "", subtitleLabSize = 0)+ theme(legend.position = "none")
ggsave(plot=mglia_volcano, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/volcanos/mglia_volcano.svg', height = 7, width = 7)
ggsave(plot=mglia_volcano, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/volcanos/mglia_volcano.png', height = 7, width = 7)

vas_volcano <- EnhancedVolcano(scTE_multi_Vas_res,
                               lab = rownames(scTE_multi_Vas_res),
                               x = 'log2FoldChange',
                               y = 'pvalue', xlim=c(-6,10), ylim=c(0,170), title = "Vascular", subtitle = "", subtitleLabSize = 0)+ theme(legend.position = "none")
ggsave(plot=vas_volcano, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/volcanos/vas_volcano.svg', height = 7, width = 7)
ggsave(plot=vas_volcano, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/volcanos/vas_volcano.png', height = 7, width = 7)


# Normalization ----
library(DESeq2)
scgene <- scTE_multi[which(startsWith(rownames(scTE_multi), 'ENS')),coldata$sample]
scgene <- scgene[which(rowSums(scgene[,coldata$sample]) > 0),]
scgene <- scgene[,coldata$sample]

scgene_dds <- DESeqDataSetFromMatrix(scgene[,rownames(coldata)], coldata, design = ~ condition)
scgene_dds <- DESeq(scgene_dds)

scTE_multi_norm <- scTE_multi[,names(scgene_dds$sizeFactor)]
scTE_multi_norm[] <- mapply('/', scTE_multi_norm[,names(scgene_dds$sizeFactor)], scgene_dds$sizeFactor)
scTE_multi_norm <- scTE_multi_norm[which(!startsWith(rownames(scTE_multi_norm), 'ENS')),]
scTE_multi_norm$TE_family <- sapply(str_split(rownames(scTE_multi_norm), ':'), `[[`, 2)
scTE_multi_norm$TE_subfamily <- sapply(str_split(rownames(scTE_multi_norm), ':'), `[[`, 1)

library(pheatmap)
coldata <- coldata[c('Emx60ctx_KO_EX', 'Emx73ctx_KO_EX', 'Emx64ctx_Ctl_EX', 'Emx65ctx_Ctl_EX',
                     'Emx60ctx_KO_IN', 'Emx73ctx_KO_IN', 'Emx64ctx_Ctl_IN', 'Emx65ctx_Ctl_IN',
                     'Emx60ctx_KO_Oligo', 'Emx73ctx_KO_Oligo', 'Emx64ctx_Ctl_Oligo', 'Emx65ctx_Ctl_Oligo',
                     'Emx60ctx_KO_OligoPre','Emx73ctx_KO_OligoPre',  'Emx64ctx_Ctl_OligoPre', 'Emx65ctx_Ctl_OligoPre',
                     'Emx60ctx_KO_Astro', 'Emx73ctx_KO_Astro', 'Emx64ctx_Ctl_Astro', 'Emx65ctx_Ctl_Astro',
                     'Emx60ctx_KO_MGlia', 'Emx73ctx_KO_MGlia', 'Emx64ctx_Ctl_MGlia', 'Emx65ctx_Ctl_MGlia',
                     'Emx60ctx_KO_Vas', 'Emx73ctx_KO_Vas', 'Emx64ctx_Ctl_Vas', 'Emx65ctx_Ctl_Vas'),]

order_ltrs <- c('IAP1-MM_I-int', 'IAPEY5_I-int', 'IAPLTR4_I', 'IAPLTR4',  'Zaphod2','RLTR46A',  
                'IAPEY3C_LTR','RLTR34B_MM',  'MER92-int',  'MMERVK10C-int',  'LTRIS4B',  'IAP1-MM_LTR',  
                'RLTR3_Mm',  'IAPLTR3-int',  'RLTR10-int',  'RLTR10',  'RLTR10C',  'RLTR44B',  'IAPEY3-int',  
                'RLTR10B2', 'IAP-d-int')
iap_mmervk <- subset(scTE_multi_norm, scTE_multi_norm$TE_subfamily %in% order_ltrs)
rownames(iap_mmervk) <- iap_mmervk$TE_subfamily
iap_mmervk <- iap_mmervk[order_ltrs,coldata$sample]
gaps <- rep(which(!duplicated(coldata$cell_type))[-1]-1, 3)

cell_type_order <- c('EX', 'Oligo', 'OligoPre', 'Astro', 'MGlia', 'IN', 'Vas')

coldata <- coldata[order(match(coldata$cell_type,cell_type_order)),]

iap_mmervk_heatmap <- pheatmap(log2(iap_mmervk[,coldata$sample]+0.5), cluster_cols = F, annotation_col = coldata[,-1],
                               gaps_col = gaps, show_colnames = F, cluster_rows = F)
ggsave(plot=iap_mmervk_heatmap, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/iap_mmervk.svg', height = 6, width = 10)
ggsave(plot=iap_mmervk_heatmap, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/iap_mmervk.png', height = 6, width = 10)

gene_transcript <- fread('/Volumes/Seagate Backup /annotation/mouse/gencode/gencode.vM20.annotation.transc.gene.tab', data.table=F, header=F)
colnames(gene_transcript) <- c('transcript_id', 'gene_id', 'gene_name', 'gene_type')

library(reshape)

scgene_norm <- counts(scgene_dds, normalize=T)
markers <- c('Grin1', 'Rbfox3','Gad2', 'Gad1', 'Olig1', 'Olig2', 'Pdgfra','Cspg4', 'Gfap', 'Ndrg2', 'Aqp4', 'Cd14', 'Cd68', 'Tek')
markers <- data.frame(scgene_norm[unique(gene_transcript[which(gene_transcript$gene_name %in% markers),'gene_id']),coldata$sample, drop=F])
markers$gene_id <- rownames(markers)
markers <- melt(markers, by=list('gene_id'))
markers <- merge(markers, unique(gene_transcript[,c('gene_id', 'gene_name')]))
markers <- merge(markers, coldata, by.x='variable', by.y='sample')
markers$gene_name <- factor(markers$gene_name, levels= c('Grin1', 'Rbfox3','Gad2', 'Gad1', 'Olig1', 'Olig2', 'Pdgfra','Cspg4', 'Gfap', 'Ndrg2', 'Aqp4', 'Cd14', 'Cd68', 'Tek'))

markers <- markers[order(match(markers$cell_type,cell_type_order)),]

markers_plot <- ggplot(data=markers, aes(x=variable, y=value, fill=cell_type)) + 
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_bar(stat = 'identity') + 
  facet_grid(gene_name ~ ., scales = 'free') 

ggsave(plot=markers_plot, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/markers.svg', height = 10, width = 7)
ggsave(plot=markers_plot, file='/Volumes/Seagate Backup /scTE/May2020/3_multimapping/plots/markers.png', height = 10, width = 7)


