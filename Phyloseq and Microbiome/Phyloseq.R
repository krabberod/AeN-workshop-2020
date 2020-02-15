### Phyloseq

oslo_fjord<-ps_dada2


sample_names(oslo_fjord)
rank_names(oslo_fjord)
sample_variables(oslo_fjord)

# If you want to subset
#oslo_fjord <- subset_samples(oslo_fjord, Select_18S_nifH =="Yes")


# Example for subsetting on taxonomic annoation. In this example we select autotrophic taxa
oslo_fjord <- subset_taxa(oslo_fjord, Division %in% c("Chlorophyta", "Dinophyta", "Cryptophyta", 
                                              "Haptophyta", "Ochrophyta", "Cercozoa"))
oslo_fjord <- subset_taxa(oslo_fjord, !(Class %in% c("Syndiniales", "Sarcomonadea")))
# Question: How could this command be used to exclude animals? 

### Normalization using median sequencing depth
total = median(sample_sums(oslo_fjord))
standf = function(x, t=total) round(t * (x / sum(x)))
oslo_fjord = transform_sample_counts(oslo_fjord, standf)

### 7.3 Bar graphs

plot_bar(oslo_fjord, fill = "Division")
plot_bar(oslo_fjord, fill = "Supergroup")

# Make the bargraph nicer by removing OTUs boundaries. This is done by adding ggplot2 modifier.
# More on ggplot later in the workshop 
plot_bar(oslo_fjord, fill = "Division") + 
  geom_bar(aes(color=Division, fill=Division), stat="identity", position="stack")

# To regroup samples that are from different fractions (could be depth/season/other metadata)
#oslo_fjord_fraction <- merge_samples(oslo_fjord, "fraction")
#plot_bar(oslo_fjord_fraction, fill = "Division") + 
#  geom_bar(aes(color=Division, fill=Division), stat="identity", position="stack")


# Plot only Chlorophyta
oslo_fjord_chloro <- subset_taxa(oslo_fjord, Division %in% c("Chlorophyta"))

plot_bar(oslo_fjord_chloro, x="Genus", fill = "Genus")+#, facet_grid = level~Class) +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")

#### 7.4 Heatmaps
plot_heatmap(oslo_fjord, method = "NMDS", distance = "bray")

#It is very very cluttered. It is better to only consider the most abundant OTUs for heatmaps. 
#For example one can only take OTUs that represent at least 20% of reads in at least one sample. 
#Remember we normalized all the sampples to median number of reads (total). 
oslo_fjord_abund <- filter_taxa(oslo_fjord, function(x) sum(x > total*0.20) > 0, TRUE)
oslo_fjord_abund

plot_heatmap(oslo_fjord_abund, method = "NMDS", distance = "bray")


# It is possible to use different distances and different multivaraite methods. 
# For example Jaccard distance and MDS and label OTUs with Class, order by Class. 
# We can also change the Palette (the default palette is a bit ugly…).

plot_heatmap(oslo_fjord_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", # 
             taxa.label = "Genus", taxa.order = "Genus", 
             trans=NULL, low="beige", high="red", na.value="beige")

# Distance measures: 

dist_methods <- unlist(distanceMethodList)
dist_methods


# 7.5 Alpha diversity 
plot_richness(oslo_fjord, measures=c("Chao1", "Shannon"))
oslo_fjord.ord <- ordinate(oslo_fjord, "NMDS", "bray")

# 7.6 Ordination
plot_ordination(oslo_fjord, oslo_fjord.ord, type="taxa", color="Class", shape= "Division", 
                title="OTUs")

# Different panels for each taxonomic division
plot_ordination(oslo_fjord, oslo_fjord.ord, type="taxa", color="Class", 
                title="OTUs", label="Class") + 
  facet_wrap(~Division, 3)

# plot_ordination(oslo_fjord, oslo_fjord.ord, type="samples", color="fraction", 
#                 shape="level", title="Samples") + geom_point(size=3)

plot_net(oslo_fjord, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.7, color="Class", point_label="Genus")
plot_net(oslo_fjord_abund, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.8, color="Class", point_label="Genus") 


