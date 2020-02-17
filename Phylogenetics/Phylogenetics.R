library(phyloseq)
#BiocManager::install("treeio")
BiocManager::install("ggtree")

library(treeio)
library(ggtree)
library(ggplot2)

OsloFjord_phyloseq<-readRDS("../Phyloseq and Microbiome/OsloFjord_phyloseq.rds")
OsloFjord_epa_tree<-read.jplace("../../Lectures/EPA/epa/RAxML_portableTree.EPA_out.jplace")

OsloFjord_epa_tree
ggtree(OsloFjord_epa_tree, layout='slanted') + coord_flip()
