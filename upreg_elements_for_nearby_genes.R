fullMMERVK <- fread('/Volumes/Seagate Backup /trim28/22.07.19_wchip/mm10_fullERVs_merge_repeatmasker.bed', data.table=F)
fullMMERVK <- fullMMERVK[-(which(duplicated(fullMMERVK[order(fullMMERVK$V7),c(1,2,3,6)]))),]
fullMMERVK$V11 <- make.unique(fullMMERVK$V11)
fullMMERVK <- fullMMERVK[,c(1,2,3,11,4,5,6)]
fullMMERVK <- as.data.frame(fullMMERVK)
rownames(fullMMERVK) <- fullMMERVK$V4

write.table(fullMMERVK, '/Volumes/Seagate Backup /trim28/22.07.19_wchip/mm10_fullERVs_merge_repeatmasker.bed', sep = '\t', col.names = F, row.names = F, quote = F)

fullMMERVK_count <- fread('/Volumes/Seagate Backup /trim28/22.07.19_wchip/fullMMERVK10C_count_matrix_2.csv', data.table = F)
coldata <- data.frame(sample=unlist(lapply(str_split(colnames(fullMMERVK_count)[7:ncol(fullMMERVK_count)], '/'), `[[`, 5)),
                      condition=unlist(lapply(str_split(colnames(fullMMERVK_count)[7:ncol(fullMMERVK_count)], '/'), `[[`, 4)))
colnames(fullMMERVK_count)[7:ncol(fullMMERVK_count)] <- unlist(lapply(str_split(colnames(fullMMERVK_count)[7:ncol(fullMMERVK_count)], '/'), `[[`, 5))

rownames(coldata) <- coldata$sample
rownames(fullMMERVK_count) <- fullMMERVK_count$Geneid

TE_dds <- DESeqDataSetFromMatrix(fullMMERVK_count[,rownames(coldata)], coldata, design = ~ condition)
TE_dds <- DESeq(TE_dds)
TE_res <- results(TE_dds)

upreg <- rownames(TE_res[which(TE_res$log2FoldChange > 0),])
upreg <- upreg[which(upreg %in% rownames(fullMMERVK))]
write.table(fullMMERVK[upreg,], col.names = F, row.names = F, quote = F, sep='\t', file='/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/upregulated/emx_upreg_fullMMERVK10C.bed')

not_upreg <- rownames(TE_res[which(!TE_res$log2FoldChange > 0 | is.na(TE_res$log2FoldChange)),])
not_upreg <- not_upreg[which(not_upreg %in% rownames(fullMMERVK))]
write.table(fullMMERVK[not_upreg,], col.names = F, row.names = F, quote = F, sep='\t', file='/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/not_upregulated/emx_notupreg_fullMMERVK10C.bed')

# ERVK ----
ERVK_count <- fread('/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/ERVK_count_matrix_2.csv', data.table = F)
coldata <- data.frame(sample=unlist(lapply(str_split(colnames(ERVK_count)[7:ncol(ERVK_count)], '/'), `[[`, 5)),
                      condition=unlist(lapply(str_split(colnames(ERVK_count)[7:ncol(ERVK_count)], '/'), `[[`, 4)))
colnames(ERVK_count)[7:ncol(ERVK_count)] <- unlist(lapply(str_split(colnames(ERVK_count)[7:ncol(ERVK_count)], '/'), `[[`, 5))

rownames(coldata) <- coldata$sample
rownames(ERVK_count) <- ERVK_count$Geneid

ERVK_dds <- DESeqDataSetFromMatrix(ERVK_count[,rownames(coldata)], coldata, design = ~ condition)
ERVK_dds <- DESeq(ERVK_dds)
ERVK_res <- results(ERVK_dds)

upreg <- rownames(ERVK_res[which(ERVK_res$log2FoldChange > 0),])
upreg <- upreg[which(upreg %in% rownames(ERVK))]
write.table(ERVK_count[upreg,c('Chr', 'Start', 'End', 'Geneid', 'Strand')], col.names = F, row.names = F, quote = F, sep='\t', file='/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/upregulated/emx_upreg_ERVK.bed')

not_upreg <- rownames(ERVK_res[which(!ERVK_res$log2FoldChange > 0 | is.na(ERVK_res$log2FoldChange)),])
not_upreg <- not_upreg[which(not_upreg %in% rownames(ERVK_res))]
write.table(ERVK_count[not_upreg,c('Chr', 'Start', 'End', 'Geneid', 'Strand')], col.names = F, row.names = F, quote = F, sep='\t', file='/Volumes/Seagate Backup /trim28/09.10.19/10_nearbygenes/not_upregulated/emx_notupreg_ERVK.bed')


