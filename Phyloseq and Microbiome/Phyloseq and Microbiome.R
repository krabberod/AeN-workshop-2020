### Intro to Phyloseq  and Microbiome ####
# Script with examples using the Phyloseq package
# Home page for Phyloseq: https://joey711.github.io/phyloseq/index.html

library("phyloseq")
library("dplyr")     # To manipulate dataframes
library("readxl")    # To read Excel files into R
library("ggplot2")   # for high quality graphics

# we can use the phyloseq object created at the end of the Dada2 pipline
# OsloFjord_phyloseq = ps_dada2 # might need to use this line
OsloFjord_phyloseq

# A detailed tutorial for importing data and creating phyloseq objects
# https://joey711.github.io/phyloseq/import-data.html

### Inspect the phyloseq object ####
# The phyloseq object contains a lot of information:
# Read abundance table, taxonomy, and metadata

sample_names(OsloFjord_phyloseq)
otu_table(OsloFjord_phyloseq)
tax_table(OsloFjord_phyloseq)
rank_names(OsloFjord_phyloseq)
sample_variables(OsloFjord_phyloseq)


# We need some metadata!
# Loading ctd data from excel

samples_df <- read_excel("cdt-data.xlsx", sheet = "CTD")
samples <- sample_data(samples_df)
sample_names(samples)<-sample_names(OsloFjord_phyloseq)
oslo_fjord <- phyloseq(otu_table(OsloFjord_phyloseq),tax_table(OsloFjord_phyloseq),samples)

### Subsetting based on taxonomy ####
# Example for subsetting on taxonomic annotation. In this example we select autotrophic taxa
oslo_fjord_autotrophs <- subset_taxa(oslo_fjord, Division %in% c("Chlorophyta", "Dinophyta", "Cryptophyta",
                                              "Haptophyta", "Ochrophyta", "Cercozoa"))
oslo_fjord_autotrophs <- subset_taxa(oslo_fjord_autotrophs, !(Class %in% c("Syndiniales", "Sarcomonadea")))


# Question: How could this command be used to exclude animals and landplants? How?

### Normalization using median sequencing depth ####
total <- median(sample_sums(oslo_fjord))
standf <- function(x, t=total) round(t * (x / sum(x)))
oslo_fjord <- transform_sample_counts(oslo_fjord, standf)

### Bar graphs ####
# Phyloseq contains wrappers for plotting functions
plot_bar(oslo_fjord, fill = "Division")
plot_bar(oslo_fjord, fill = "Supergroup")

# Make the bargraph nicer by removing OTUs boundaries.
# This is done by adding ggplot2 modifier.
# More on ggplot later in the workshop
plot_bar(oslo_fjord, fill = "Division") +
  geom_bar(aes(color=Division, fill=Division), stat="identity", position="stack")

# To regroup samples that are from different fractions (could be depth/season/other metadata)
oslo_fjord_year <- merge_samples(oslo_fjord, "Year")
plot_bar(oslo_fjord_year, fill = "Supergroup") +
  geom_bar(aes(color=Supergroup, fill=Supergroup), stat="identity", position="stack")


# Plot only Chlorophyta
oslo_fjord_chloro <- subset_taxa(oslo_fjord, Division %in% c("Chlorophyta"))

plot_bar(oslo_fjord_chloro, x="Genus", fill = "Genus")+#, facet_grid = level~Class) +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")

#### Heatmaps ####
plot_heatmap(oslo_fjord, method = "NMDS", distance = "bray")

#It is very very cluttered. It is better to only consider the most abundant OTUs for heatmaps.
#For example one can only take OTUs that represent at least 10% of reads in at least one sample.
#Remember we normalized all the sampples to median number of reads (total).
oslo_fjord_abund <- filter_taxa(oslo_fjord, function(x) sum(x > total*0.10) > 0, TRUE)
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

### Simple diversity analysis
# Example of simple analyses
# More on diversity analysis later today
plot_richness(oslo_fjord, measures=c("Chao1", "Shannon"))

# See
?plot_richness
# for more detail, you can experiment with different functions for richness

### Ordination
oslo_fjord.ord <- ordinate(oslo_fjord_abund, "NMDS", "bray")


plot_ordination(oslo_fjord, oslo_fjord.ord, type="taxa", color="Class", shape= "Division",
                title="OTUs")

# Different panels for each taxonomic division
plot_ordination(oslo_fjord, oslo_fjord.ord, type="taxa", color="Class",
                title="OTUs", label="Class") +
  facet_wrap(~Division, 3)

# ordination of samples:
plot_ordination(oslo_fjord, oslo_fjord.ord, type="samples", color="fraction",
                 shape="level", title="Samples") + geom_point(size=3)

plot_net(oslo_fjord_abund, distance = "(A+B-2*J)/(A+B)", type = "taxa",
         maxdist = 0.7, color="Class", point_label="Genus")
plot_net(oslo_fjord_abund, distance = "(A+B-2*J)/(A+B)", type = "taxa",
         maxdist = 0.8, color="Class", point_label="Genus")

#### Microbiome ####
# Another package with many interesting functions is Microbiome. It is an extension to phyloseq
# https://microbiome.github.io/tutorials/

library("microbiome")

# for instance prevalence
# For each OTU, the fraction of samples where a given OTU is detected.
# The output ca be given as a percentage.
head(prevalence(oslo_fjord,detection = 1/100, sort = TRUE))

# Relative abundance:
oslo_fjord_rel <- transform(oslo_fjord, "compositional")

# Looking for "core members"
core_members(oslo_fjord_rel, detection = 0, prevalence = 25/100)

# making a phyloseq object of the core
pseq.core <- core(oslo_fjord_rel, detection = 0, prevalence = .5)

tax_table(pseq.core)
otu_table(pseq.core)

#
library(RColorBrewer)
p <- plot_core(oslo_fjord_rel, plot.type = "heatmap",
               prevalences = prevalences,
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .2, horizontal = TRUE)
print(p)


# Example of complex aggregating of abundance on higher level and plotting with different colors
oslo_fjord_autotrophs.melt<- oslo_fjord_autotrophs %>% psmelt() %>% arrange(Kingdom) %>%
  {aggregate(.$Abundance,by=list(.$Supergroup,.$Division,.$Class,.$Order), FUN=sum)}

oslo_fjord_autotrophs.melt %>%
  ggplot(aes(x = Group.2, y = x,fill = Group.3)) +
  #facet_grid(. ~ Group.1) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = sample(colorRampPalette(brewer.pal(8,"Dark2"))(20)),name="Class")+
  theme_bw(base_size = 8)+
  theme(legend.key.size = unit(0.2, "cm")) +
  theme(axis.title.x = element_blank() ,axis.text.x  = element_text(angle=45, vjust=0.5, colour="darkgrey")) +
  ylab("Relative Abundance\n") +
  ggtitle("All groups")
