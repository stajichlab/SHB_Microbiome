#!/usr/bin/env R

library(ape)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(plyr)
library(data.table)
library(tidyr)
library(tidyverse)
library(cowplot)
library(microViz)
# read in metadata

meta = read.table("metadata_FD_202307_3Pool.tsv",header=TRUE,row.names=1,sep="\t",stringsAsFactors=FALSE)

sampleData <- sample_data(meta)


# I had to fix this file by hand to get rid of leading '#'
otus <- read.table("results/amptk.dada2/FD_202307_3Pool.final.txt",header=T,sep="\t",row.names=1)

otumat <- as(as.matrix(otus), "matrix")
OTU = otu_table(otumat, taxa_are_rows = TRUE)

head(OTU)

taxa <- read.table("results/amptk.dada2/FD_202307_3Pool.ASVs.taxonomy.fix.txt", header=T,sep="\t",row.names=1)
taxmat <- as(as.matrix(taxa),"matrix")
TAX = tax_table(taxmat)

treefile = "results/amptk.dada2/FD_202307_3Pool.ASVs.tree.phy"
tree = read.tree(treefile)

# phyloseq object
physeq = phyloseq(OTU,TAX,sampleData,tree)


physeq.prune = prune_taxa(taxa_sums(physeq) > 1, physeq)
# some QC
readcount = data.table(as(sample_data(physeq.prune), "data.frame"),
                       TotalReads = sample_sums(physeq.prune), 
                       keep.rownames = TRUE)
setnames(readcount, "rn", "SampleID")

#For plotting, use command below.

SeqDepth = ggplot(readcount, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")

png("./Figures/Fungal_SequencingDepth.png", units="in", width = 5.8, height = 5.8, res = 300 )
SeqDepth
pdf("./Figures/Fungal_SequencingDepth.pdf", width = 5.8, height = 5.8 )
SeqDepth
dev.off()
dev.off()


pdf("./Figures/Barchart_Pamericanus.pdf")
physeq.prune.Pamericans <- ps_filter(physeq.prune,host=="persea_americana")
psTopNOTUs = names(sort(taxa_sums(physeq.prune.Pamericans), TRUE)[1:5000])
pstop.prune = prune_taxa(psTopNOTUs, physeq.prune.Pamericans)
p<-plot_bar(pstop.prune, x = "Sample", y = "Abundance", fill ="Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  theme_cowplot(12)
p

p <- plot_bar(pstop.prune, x="Class", fill = "Class", facet_grid = end_status~status) +
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")
p

p<-plot_bar(pstop.prune, x = "Sample", y = "Abundance", fill ="Class") + 
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack") +
  theme_cowplot(12)
p

physeq.prune.Pracemosa <- ps_filter(physeq.prune,host=="platanus_racemosa")
pdf("./Figures/Barchart_Pracemosa.pdf")
psTopNOTUs = names(sort(taxa_sums(physeq.prune.Pracemosa), TRUE)[1:5000])
pstop.prune = prune_taxa(psTopNOTUs, physeq.prune.Pracemosa)
p<-plot_bar(pstop.prune, x = "Sample", y = "Abundance", fill ="Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  theme_cowplot(12)
p


psTopNOTUs = names(sort(taxa_sums(physeq.prune.Pracemosa), TRUE)[1:5000])
pstop.prune = prune_taxa(psTopNOTUs, physeq.prune.Pracemosa)
p<-plot_bar(pstop.prune, x = "Sample", y = "Abundance", fill ="Class") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  theme_cowplot(12)
p

# set.seed(1)
# physeq.prune.rarefy = rarefy_even_depth(physeq.prune, sample.size = 1000, 
#                                         replace = FALSE, trimOTUs = FALSE)
# physeq.prune.rarefy
# 
# ps.dist = phyloseq::distance(physeq.prune, "bray")
# 
# physeq.prune.ord <- ordinate(physeq.prune, "PCoA", "bray")
# adonis(ps.dist ~ end_status, as(sample_data(physeq.prune),"data.frame"))
# 
# psendstatus = plot_ordination(physeq.prune, physeq.prune.ord, 
#                               type = "samples", color = "end_status") + 
#   theme_bw() + ggtitle("Fungal Beta Diversity (PCoA) in P. americana by End Status") + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   stat_ellipse(geom = "polygon", alpha = 1/12, aes(fill = end_status)) +
#   annotate("text", x = -0.3, y = 0.38, label = "PERMANOVA, p = 0.001")
# psendstatus
# 
# 
# physeq.prune.x = prune_taxa(taxa_sums(physeq.prune.rarefy) > 1, physeq.prune.rarefy)
# physeq.prune.x.ord <- ordinate(physeq.prune.x, "PCoA", "bray")
# 
# psendstatus.rarefy = plot_ordination(physeq.prune.x, physeq.prune.x.ord, 
#                               type = "samples", color = "end_status") + 
#   theme_bw() + ggtitle("Fungal Beta Diversity (PCoA) in P. americana by End Status W Rarefaction") + 
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   stat_ellipse(geom = "polygon", alpha = 1/12, aes(fill = end_status)) +
#   annotate("text", x = -0.3, y = 0.38, label = "PERMANOVA, p = 0.001")
