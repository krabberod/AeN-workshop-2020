### Brief introduction to Network generation in R
# Again we make use of the phyloseq object OsloFjord_phyloseq
# 
OsloFjord_phyloseq <- readRDS("OsloFjord_phyloseq.rds")

#library(devtools)
#install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library(Matrix)
library(igraph)
library(phyloseq)



library(tidyverse)
library(microbiome)
library(devtools)
library(seqtime)
library(ProNet)
library(standardize)


### Make network with Sparcc: ####
# Relation bewteen samples:
OsloFjord.sparcc.samples<-sparcc(otu_table(OsloFjord_phyloseq), iter=20,inner_iter=100)

## Define  threshold for SparCC correlation matrix for the graph
sparcc.graph <- abs(OsloFjord.sparcc$Cor) >= 0.1
sparcc.graph
diag(sparcc.graph) <- 0
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
# Create igraph objects
ig.sparcc <- adj2igraph(sparcc.graph)
# simple plotting
plot(ig.sparcc)
plot(ig.sparcc,vertex.size=1)

#### 
#OsloFjord.sparcc.otu<-sparcc(t(otu_table(OsloFjord.core)), iter=20,inner_iter=100)


## Define  threshold for SparCC correlation matrix for the graph
sparcc.graph.otu <- abs(OsloFjord.sparcc.otu$Cor) >= 0.5
sparcc.graph.otu

diag(sparcc.graph.otu) <- 0
sparcc.graph.otu <- Matrix(sparcc.graph.otu, sparse=TRUE)

# Create igraph objects
ig.sparcc.otu <- adj2igraph(sparcc.graph.otu)

# simple plotting
plot(ig.sparcc.otu)



plot_network(ig.sparcc.otu) # update to get better
plot(ig.sparcc.otu)


OsloFjord_phyloseq.boot<-sparccboot(t(otu_table(OsloFjord_phyloseq)), R =1)

row.names(OsloFjord.sparcc.otu$Cor)<-row.names(otu_table(OsloFjord.core))

colnames(OsloFjord.sparcc.otu$Cor)<-row.names(otu_table(OsloFjord.core))
sparcc.graph.otu <- abs(OsloFjord.sparcc.otu$Cor) >= 0.5

diag(sparcc.graph.otu) <- 0
sparcc.graph.otu<- Matrix(sparcc.graph.otu, sparse=TRUE)

sparcc.graph.otu <- adj2igraph(sparcc.graph.otu, vertex.attr = list(name=row.names(otu_table(OsloFjord.core))))
am.coord <- layout.kamada.kawai(sparcc.graph.otu)
plot(sparcc.graph.otu, layout=am.coord, vertex.size=5, main="sparcc")


dd.sparcc <- degree.distribution(sparcc.graph.otu)
hist(dd.sparcc)

write.graph(sparcc.graph.otu, "sparcc.graph.otu.graphml", format = "graphml")




#### ####

pargs <- list(rep.num=50, seed=10010, ncores=4) 
OsloFjord_phyloseq <- spiec.easi(OsloFjord_phyloseq, method='mb', lambda.min.ratio=1e-3, nlambda=30,
                                 sel.criterion='bstars', pulsar.select=TRUE, pulsar.params=pargs)




rank_names(OsloFjord_phyloseq)
tax_table(OsloFjord_phyloseq)

OsloFjord.core <- core(OsloFjord_phyloseq, detection = 0, prevalence = .2)



# Zure 2017, Genomwalker
source("~/Dropbox/Projects/Scripts/R_scripts/Trimming_algorithm.R")

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
trans_red<-makeTransparent("#e31a1c", 70)
trans_green<-makeTransparent("#33a02c", 70)


#x<- as(otu_table(subset_taxa(OsloFjord_phyloseq, Supergroup == "Alveolata")), "matrix")
#y<- as(otu_table(subset_taxa(OsloFjord_phyloseq, Supergroup == "Stramenopiles")), "matrix")
#z<- as(otu_table(subset_taxa(OsloFjord_phyloseq, Supergroup == "Opisthokonta")), "matrix")


OTU<-otu_table(rbind(y,z), taxa_are_rows = TRUE)
OTU

# Adding back the taxonomic information: the phyloseq command merges the correct data
oslofjord_reduced<-phyloseq(OTU,tax_table(OsloFjord_phyloseq))

#### Spiec-easi for sparce composonality OTUs all data ####
#OsloFjord.core.mb <- spiec.easi(oslofjord_reduced, method='mb', lambda.min.ratio=1e-5, 
#                         nlambda=30, icov.select.params=list(rep.num=60))






se <- spiec.easi(OsloFjord_phyloseq, method='mb', lambda.min.ratio=1e-2, nlambda=15)
ig2.mb <- adj2igraph(getRefit(se),  vertex.attr=list(name=taxa_names(OsloFjord_phyloseq)))

plot_network(ig2.mb, amgut2.filt.phy, type='taxa', color="Rank3")

OsloFjord.core.mb.rar<-spiec.easi(oslofjord_reduced, verbose = TRUE, sel.criterion = 'mb')

seRes<-OsloFjord.core.mb

row.names(temp) = c("Lambda", "eBIC", "Sparsity")
temp = rbind(round(seRes$lambda, 4), round(seRes$ebic.score, 0), round(seRes$sparsity, 2))
row.names(temp) = c("Lambda", "eBIC", "Sparsity")
knitr::kable(temp,  col.names = as.character(1:dim(temp)[2]))


# for the rar
OsloFjord.core.mb<-spiec.easi(OsloFjord.core, method='mb', lambda.min.ratio=1e-4, 
                           nlambda=20, icov.select.params=list(rep.num=50,ncores = 8))


#### This takes time!
OsloFjord.core.glasso <- spiec.easi(OsloFjord.core, method='glasso', lambda.min.ratio=1e-2,
                           nlambda=20, icov.select.params=list(rep.num=50))
# 
# Turning it into a graph and look at the segdes

OsloFjord.core.mb.graph <- adj2igraph(OsloFjord.core.glasso$refit)#,  vertex.attr=list(name=taxa_names(OsloFjord.core)))
plot_network(OsloFjord.core.mb.graph,soil.phyloseq, label = "Class", color = "Domain", type="taxa")

betaMat<-as.matrix(symBeta(getOptBeta(OsloFjord.core.mb)))
elist.mb <- summary(symBeta(getOptBeta(OsloFjord.core.mb), mode='maxabs')) 

# Getting weights
elist.mb <- summary(symBeta(getOptBeta(OsloFjord.core.mb), mode='maxabs')) #OBS This is not the correct one! 
E(OsloFjord.core.mb.graph)$weight<-elist.mb[,3]

plot_network(OsloFjord.core.mb.graph,soil.phyloseq, label = "Genus", color = "Domain", type="taxa")#, shape = "Depth")
soil.layout <- layout.fruchterman.reingold(OsloFjord.core.mb.graph)
#vsize <- rowMeans(clr(bactOTU.select, 1))+6 #Clr is the centered log-ratio function

plot(OsloFjord.core.mb.graph, layout=soil.layout, main="OsloFjord.core.mb", vertex.size=2)
write.graph(OsloFjord.core.mb.graph,"OsloFjord.core.mb.graphml", format = "graphml")
# Pruned edges based on Zure 2017
source("~/Dropbox/Projects/Scripts/R_scripts/Network_ConnectedComponent.r")

G<-OsloFjord.core.mb.graph
E(G)$weight<-elist.mb[,3]
#G.pos<-delete.edges(G,which(E(G)$weight<0))
#G<-G.pos
G<-delete.vertices(G,which(degree(G)<1))
weights <- unique(sort(E(G)$weight, decreasing = FALSE))
out<- binary_search_filter(g = G, low = 1, high = length(weights))$graph

# stop at 0.005939436, removes about 70 pos edges... 

soil.graph.pruned<-delete.edges(OsloFjord.core.mb.graph,which(E(OsloFjord.core.mb.graph)$weight<0.02))
#plot(soil.graph.pruned,  main="spiec.easi", vertex.size=4)
#plot_network(soil.graph.pruned,soil.phyloseq)
degree(ig.sparcc)
simplify(ig.sparcc)


### FOR RARIEFIED: 
OsloFjord.core.mb.rar.graph <- adj2igraph(OsloFjord.core.mb.rar$refit,  vertex.attr=list(name=taxa_names(soil.phyloseq)))


plot_network(ig.sparcc,OsloFjord_phyloseq, label = "Class", color = "Domain", type="taxa")

# plot_network(OsloFjord.core.mb.graph,soil.phyloseq, label = "Class", color = "Domain", type="taxa")


betaMat<-as.matrix(symBeta(getOptBeta(OsloFjord.core.mb.rar)))
elist.mb <- summary(symBeta(getOptBeta(OsloFjord.core.mb.rar), mode='maxabs')) 

# Getting weights
elist.mb <- summary(symBeta(getOptBeta(OsloFjord.core.mb.rar), mode='maxabs')) #OBS This is not the correct one! 
E(OsloFjord.core.mb.rar.graph)$weight<-elist.mb[,3]

plot_network(OsloFjord.core.mb.rar.graph,soil.phyloseq, label = "Genus", color = "Domain", type="taxa")#, shape = "Depth")
soil.layout <- layout.fruchterman.reingold(OsloFjord.core.mb.rar.graph)
#vsize <- rowMeans(clr(bactOTU.select, 1))+6 #Clr is the centered log-ratio function

plot(OsloFjord.core.mb.rar.graph, layout=soil.layout, main="OsloFjord.core.mb.rar", vertex.size=2)
write.graph(OsloFjord.core.mb.rar.graph,"OsloFjord.core.mb.rar.graphml", format = "graphml")
# Pruned edges based on Zure 2017
source("~/Dropbox/Projects/Scripts/R_scripts/Trimming_algorithm.R")

G<-OsloFjord.core.mb.rar.graph
E(G)$weight<-elist.mb[,3]
G.pos<-delete.edges(G,which(E(G)$weight<0))
G<-G.pos
G<-delete.vertices(G,which(degree(G)<1))
weights <- unique(sort(E(G)$weight, decreasing = FALSE))
out<- binary_search_filter(g = G, low = 1, high = length(weights))$graph

# stop at 0.005939436, removes about 70 pos edges... 

soil.graph.pruned<-delete.edges(OsloFjord.core.mb.rar.graph,which(E(OsloFjord.core.mb.rar.graph)$weight<0.02))
#plot(soil.graph.pruned,  main="spiec.easi", vertex.size=4)
#plot_network(soil.graph.pruned,soil.phyloseq)
degree(soil.graph.pruned)
simplify(soil.graph.pruned)




#### From Faust et al online material for network course, this is not necessary for our purpose ####
(positive<-length(betaMat[betaMat>0])/2)
(negative<-length(betaMat[betaMat<0])/2)
(total<-length(betaMat[betaMat!=0])/2)

# For coloring of edges
otu.ids<-colnames(OsloFjord.core.mb$data)

edges <- E(ig.sparcc.otu)
edge.colors=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(ig.sparcc.otu,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta>0){
    edge.colors=append(edge.colors,trans_green)
  }else if(beta<=0){
    edge.colors=append(edge.colors,trans_red)
  }
}

E(OsloFjord.core.mb.graph)$color=edge.colors

spiec.graph.b=OsloFjord.core.mb.graph
nodenames=V(spiec.graph.b)$name
#V(spiec.graph.b)$name=getTaxonomy(nodenames, taxa.f, useRownames=TRUE)
E(spiec.graph.b)$arrow.size=5
V(spiec.graph.b)$color="white"
V(spiec.graph.b)$frame.color="black"
tkplot(spiec.graph.b)

set.seed(65345690)
visualization(OsloFjord.core.mb.graph,edge.color=E(OsloFjord.core.mb.graph)$color, node.size = 1,
              layout="lgl",
              node.label.size=0.5, 
              node.fill.color = "red",
              node.label.color = "black", node.label =NA)# = V(cluster2)$Combined)

otu.ids=colnames(OsloFjord.core.mb$data)
edges=E(soil.graph)
filtered.edges=c()
for(e.index in 1:length(edges)){
  adj.nodes=ends(soil.graph,edges[e.index])
  xindex=which(otu.ids==adj.nodes[1])
  yindex=which(otu.ids==adj.nodes[2])
  beta=betaMat[xindex,yindex]
  if(beta<0){
    filtered.edges=c(filtered.edges,edges[e.index])
  }
}
soil.graph.pos=delete_edges(soil.graph, filtered.edges)
plot_network(soil.graph.pos)

#write.graph(soil.graph,"soil.graph.spieceasi.ncol.txt",format="ncol")



#### Sparcc seems to give a better result ####
#sparcc.soil <- sparcc(t(otu_table(soil.phyloseq)), iter=20,inner_iter=100)

#sparcc.soil.rar <- sparcc(t(otu_table(oslofjord_reduced)), iter=20,inner_iter=100)


elist.sparcc <- summary(sparcc.graph*sparcc.soil$Cor)
E(sparcc.soil.graph)$weight<-elist.sparcc[,3]



set.seed(23948)




G<-sparcc.graph.otu
#E(G)$weight<-elist.mb[,3]
G.pos<-delete.edges(G,which(E(G)$weight<0))
G<-G.pos

G<-delete.vertices(G,which(degree(G)<1))
weights <- unique(sort(E(G)$weight, decreasing = FALSE))
out<- binary_search_filter(g = G, low = 1, high = length(weights))$graph






#### Sparcc for depths D1 ####

row.names(sparcc.soil.d1$Cor)<-row.names(otu_table(soil.phyloseq.d1))
colnames(sparcc.soil.d1$Cor)<-row.names(otu_table(soil.phyloseq.d1))
hist(abs(sparcc.soil.d1$Cor))

sparcc.graph.d1 <- abs(sparcc.soil.d1$Cor) >= 0.5
diag(sparcc.graph.d1) <- 0
sparcc.graph.d1 <- Matrix(sparcc.graph.d1, sparse=TRUE)
sparcc.soil.graph.d1 <- adj2igraph(sparcc.graph.d1, vertex.attr = list(name=row.names(otu_table(soil.phyloseq.d1))))

elist.sparcc.d1 <- summary(sparcc.graph.d1*sparcc.soil.d1$Cor)
dd.sparcc.d1 <- degree.distribution(sparcc.soil.graph.d1)
#set.seed(23948)
plot_network(sparcc.soil.graph.d1)#, 
# Prune by Zure 2017! 

E(sparcc.soil.graph.d1)$weight<-elist.sparcc.d1[,3]
G<-sparcc.soil.graph.d1
#G.pos<-delete.edges(G,which(E(G)$weight<0))
#G<-G.pos
G<-delete.vertices(G,which(degree(G)<1))
weights <- unique(sort(E(G)$weight, decreasing = FALSE))
out<- binary_search_filter(g = G, low = 1, high = length(weights))$graph
write.graph(sparcc.soil.graph.d1, "sparcc.soil.graph.d1.graphml", format = "graphml")











