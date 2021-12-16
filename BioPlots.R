#If not installed, install packages. DESeq2 requires Bioconductor: https://www.bioconductor.org/install/
#Current version of Bioconductor install commands below (as of 12/16/2022)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

library(DESeq2)
library(RColorBrewer)
library(tidyverse)
library(gridExtra)
library(EnhancedVolcano)
library(pheatmap)
library(reshape2)
library(genefilter)
library(ggplot2)
library(dplyr)

#load and sort count file
data <- read.table("/Volumes/directory/gene.cntTable", header = T)

#filter your data
min_read <- 9
data <- data[apply(data,1,function(x){max(x)}) > min_read,]

#organize my data
data <- remove_rownames(data)
data <- column_to_rownames(data, var = "gene.TE")
colnames(data) <- sapply(str_split(as.character(colnames(data)), "\\."), '[', 7)

#load and sort sample information
sampleInfo <- read.csv("/Volumes/directory/sampledata.csv",header=T, row.names = 1)

#Get rid of NA and Other from my sample information
sampleInfo <- sampleInfo[complete.cases(sampleInfo$samplecolumn), ]
sampleInfo <- sampleInfo[!grepl("Other",sampleInfo$samplecolumn),]

#order to match
countsSub <- data[,colnames(data) %in% sort(rownames(sampleInfo))]
i <- match(rownames(sampleInfo), colnames(countsSub), nomatch = 0)
countsSub <- countsSub[,i]
identical(colnames(countsSub), as.character(rownames(sampleInfo))) #double checking
setdiff(union(colnames(countsSub), as.character(rownames(sampleInfo))), intersect(colnames(countsSub), as.character(rownames(sampleInfo))))

#Run DESeq2 on the filtered/organized data
suppressPackageStartupMessages(library(DESeq2))
dds <- DESeqDataSetFromMatrix(countData = countsSub, colData = sampleInfo, design = ~var1 + var2 + var3)
dds <- DESeq(dds)

#normalized
norm_dds <- counts(dds, normalized = TRUE)
res <- results(dds)
res.df <- rownames_to_column(as.data.frame(res), var = "Id")
write_tsv(res.df, "/Volumes/directory/name.tsv")

#Contrast example
caseres <- results(dds, contrast=c('var1', 'Case', 'Control'), )
rld <- vst(dds)

#print results
caseres.df <- rownames_to_column(as.data.frame(caseres), var = "Id")
write_tsv(caseres.df, "/Volumes/directory/name.tsv")

#Run and print LRT
lrtdds <- DESeq(dds, test="LRT", reduced=~var1 + var2)
lrtres <- results(lrtdds)
lrtres.df <- rownames_to_column(as.data.frame(lrtres), var = "Id")
write_tsv(lrtres.df, "/Volumes/directory/LRT.tsv")

#Interaction condition example
dds$group <- factor(paste0(dds$var1, dds$var2))
design(dds) <- ~ group
dds2 <- DESeq(dds)

caseres2 <- results(dds2, contrast=c("group", "Case1", "Case2"))
caseres2.df <- rownames_to_column(as.data.frame(casesexres), var = "TE_id")
write_tsv(casesexres.df, "/Volumes/directory/interaction.tsv")

#boxplot comparison of gene counts by padj (can adjust gene)
df <- plotCounts(dds, normalized = TRUE, gene='LSAU:Satellite:Satellite', intgroup= c("var1", "var2"), 
                returnData=TRUE)
ggplot(df, aes(x=case_control_other_at_baseline, y=count)) + 
  geom_boxplot(alpha = 0, color = 'black') +
  geom_point(position=position_jitter(w=0.1,h=0), color = "red") +
  scale_y_log10(breaks=c(25,100,400,1000,9000)) +
  xlab('Condition') +
  ylab('Normalized Counts') +
  ggtitle('LSAU') +
  theme(axis.text = element_text(size = rel(1.75))) +
  theme(axis.title = element_text(size = rel(1.75))) +
  theme(legend.text = element_text(size = rel(1.75))) +
  theme(legend.title = element_text(size = rel(1.75))) +
  theme(title = element_text(size = rel(1.5))) +
  theme_bw(base_size = 32)+
  theme(axis.text.y = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  theme(axis.title.y = element_text(angle = 45, vjust = 0.75, hjust=1))
ggsave("/Volumes/directory/LSAU.png", width = 16, height = 8)

##To filter by two different groups follow the instructions
#separate condition only, for example

Ca <- filter(sampleInfo, var1 == 'Case')
CacountsSub <- data[,colnames(data) %in% sort(rownames(Ca))]
i <- match(rownames(Ca), colnames(CacountsSub), nomatch = 0)
CacountsSub <- CacountsSub[,i]
identical(colnames(CacountsSub), as.character(rownames(Ca))) #double checking
Cadds <- DESeqDataSetFromMatrix(countData = CacountsSub, colData = Ca, design = ~var2 + var3)
Cadds <- DESeq(Cadds)

#export results

Cares <- results(Cadds)
Cares.df <- rownames_to_column(as.data.frame(Cares), var = "TE_id")
write_tsv(Cares.df, "/Volumes/directory/Caseresultsonly.tsv")

#create heatmap using top twenty normalized counts
select <- order(rowMeans(assay(rld)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("var1", "var2")])
colnames(df) <- c("Condition", "Sex")
pheatmap(scale(assay(rld))[select,], cluster_rows=FALSE, show_colnames=FALSE,
         cluster_cols=TRUE, annotation_col=df, annotation_legend = TRUE)
ggsave("/Volumes/directory/HeatmapTopTwenty.png", width = 16, height = 8)

#create volcano plot using normalized counts
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                ylim = c(0, 2),
                title = 'Expression Levels',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 5.0,
                ylab = '-log(Adjusted p-value)',
                drawConnectors = TRUE,
                legendLabels=c('Not sig.','Log (base 2) FC','-log(p-adjusted)',
                               '-log(p-adjusted) & Log (base 2) FC'),
                )
ggsave("/Volumes/directory/VolcanoPlotAll.png", width = 16, height = 8)

#Count your participants by applicable variables and make a pretty image
participants <- table(sampleInfo$var1, sampleInfo$var2)
participants <- as.data.frame(participants, header = TRUE)
names(participants)[names(participants) == 'Var1'] <- 'Sex'
names(participants)[names(participants) == 'Var2'] <- 'Condition'
names(participants)[names(participants) == 'Freq'] <- 'Participant Number'
grid.table(participants)