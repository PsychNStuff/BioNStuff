library(UpSetR)
library(dplyr)
library(tidyverse)

##This will make an upset plot by different biofluids (in this case, plasma and cerebrospinal fluid). Use different results from DESeq2 for each biofluid.

#Filter my data by up/down regulated
res.df <- filter(res.df, padj < 0.05)
CSF.counts.up <- filter(res.df, log2FoldChange > log2(1.5))
CSF.upregulated <- filter(countsSub, rownames(countsSub) %in% CSF.counts.up$gene_id)
CSF.up.upset <- CSF.upregulated
CSF.up.upset[CSF.up.upset > 0] <- 1
CSF.up.upset <- CSF.up.upset %>% mutate(CSFupGeneSums = rowSums(CSF.up.upset))
CSF.up <- rownames_to_column(CSF.up.upset, var = "Geneid")

CSF.counts.down <- filter(res.df, log2FoldChange < -log2(1.5))
CSF.downregulated <- filter(countsSub, rownames(countsSub) %in% CSF.counts.down$gene_id)
CSF.down.upset <- CSF.downregulated
CSF.down.upset[CSF.down.upset > 0] <- 1
CSF.down.upset <- CSF.down.upset %>% mutate(CSFdownGeneSums = rowSums(CSF.down.upset))
CSF.down <- rownames_to_column(CSF.down.upset, var = "Geneid")

forupset.csf <- full_join(x=CSF.down, y=CSF.up)

#Filter plasma data
res.df <- filter(res.df, padj < 0.05)
Plasma.counts.up <- filter(res.df, log2FoldChange > log2(1.5))
Plasma.upregulated <- filter(countsSub, rownames(countsSub) %in% Plasma.counts.up$gene_id)
Plasma.up.upset <- Plasma.upregulated
Plasma.up.upset[Plasma.up.upset > 0] <- 1
Plasma.up.upset <- Plasma.up.upset %>% mutate(PlasmaupGeneSums = rowSums(Plasma.up.upset))
Plasma.up <- rownames_to_column(Plasma.up.upset, var = "Geneid")

Plasma.counts.down <- filter(res.df, log2FoldChange < -log2(1.5))
Plasma.downregulated <- filter(countsSub, rownames(countsSub) %in% Plasma.counts.down$gene_id)
Plasma.down.upset <- Plasma.downregulated
Plasma.down.upset[Plasma.down.upset > 0] <- 1
Plasma.down.upset <- Plasma.down.upset %>% mutate(PlasmadownGeneSums = rowSums(Plasma.down.upset))
Plasma.down <- rownames_to_column(Plasma.down.upset, var = "Geneid")

forupset.plasma <- full_join(x=Plasma.down, y=Plasma.up)

#join csf and plasma data together
forupset <- full_join(forupset.plasma, forupset.csf, by = c("Geneid"))

#add up genesums
forupset <- column_to_rownames(forupset, var = "Geneid")
upset.all <- transmute(forupset, PlasmaUpregulated = PlasmaupGeneSums > 1,
                                  PlasmaDownregulated = PlasmadownGeneSums > 1,
                                  CSFUpregulated = CSFupGeneSums > 1,
                                  CSFDownregulated = CSFdownGeneSums > 1)
#make upset plot
upset.all[upset.all > 0] <- 1
upset.all[is.na(upset.all)] = 0
upset(upset.all, nsets = 4, order.by = "freq")
options(ragg.max_dim = 1000000000000)
ggsave("pathtodirectory/filename.png", width = 50, height = 50, limitsize = F)
