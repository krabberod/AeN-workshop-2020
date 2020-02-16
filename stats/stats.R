### STATS

### Initiate ###
# load libraries
library("data.table")
library("vegan")
library("dplyr")
library("readxl")
library("ggplot2")
library("R.utils") # nice library if you want to make and use own functions
library("dendextend")

# define directories and load files
figs_dir <- "figs/"  # figures

ps_dada2<-readRDS(file.path(dada2_dir, "OsloFjord_phyloseq.rds"))
OsloFjord_phyloseq <- ps_dada2

# check phyloseq object
sample_names(OsloFjord_phyloseq)
otu_table(OsloFjord_phyloseq)
tax_table(OsloFjord_phyloseq)
rank_names(OsloFjord_phyloseq)
sample_variables(OsloFjord_phyloseq)

# add metadata
metadata <- read_excel("cdt-data.xlsx", sheet = "CTD")
samples <- sample_data(metadata)
sample_names(samples)<-sample_names(OsloFjord_phyloseq)
oslo_fjord <- phyloseq(otu_table(OsloFjord_phyloseq),tax_table(OsloFjord_phyloseq),samples)



### Rarefaction curves ###
# Background. Rarefaction is a technique to assess species richness from the results of sampling. Rarefaction allows the calculation of species richness for a given number of individual samples, based on the construction of a so-called rarefaction curve. This curve is a plot of the number of species as a function of the number of sequences. It is created by randomly re-sampling the pool of N samples multiple times and then plotting the average number of species found in each sample (1,2, ... N). Thus, rarefaction generates the expected number of species in a small collection of n individuals (or n samples) drawn at random from the large pool of N samples. Rarefaction curves generally grow rapidly at first, as the most common species are found, but the curves plateau as only the rarest species remain to be sampled.

rarecurve(otu_table(oslo_fjord), step = 100, xlab = "Sample Size", ylab = "OTUs", label = TRUE)

# Question A: From looking at species richness estimates and rarefaction curves, what can you conclude about the sequencing effort? Do you think 100,000 sequences per sample are enough to capture the full diversity of the microbiome?
# Question B: If we were to increase the number of sequences from 10'000 to say 10'000'000, do you see any potential problem with the algorithms in R that we use? Would the analysis take 100x more time to run than it does now? How could you imagine solving this problem?

### Get phylum stats ###

# Make functions to summarize stats
# These functions take a phyloseq object and return phyloseq data (OTUs are now taxonomic ranks)

summarize_taxa <- function(counts, taxonomy) {
  if(is.matrix(taxonomy)) {
    message('multiple taxonomies')
    apply(taxonomy, 2, summarize_taxa, counts = counts, .dims = TRUE)
  } else if(is.matrix(counts)) {
    message('multiple counts')
    require('plyr')
    apply(counts, 2, summarize_taxa, taxonomy = taxonomy)
  } else {
    message('summarize')
    tapply(counts, taxonomy, sum)
  }
}

phyloseq_summarize_taxa <- function(phylo_seq_object, taxonomic_rank = rank_names(phylo_seq_object)[1], errorIfNULL = TRUE) {
  taxa <- as(tax_table(phylo_seq_object, errorIfNULL)[, taxonomic_rank], 'character')
  sum_tax_table <- summarize_taxa(as(otu_table(phylo_seq_object), 'matrix'), taxa)
  phyloseq(otu_table(sum_tax_table, taxa_are_rows = TRUE), sample_data(phylo_seq_object, FALSE))
}

# Analyze at division level
division_table <- phyloseq_summarize_taxa(oslo_fjord, 'Division')
sort(taxa_sums(division_table)/sum(taxa_sums(division_table)))
sort(taxa_sums(division_table)/sum(taxa_sums(division_table)), decreasing = TRUE)

### Make some distribution plots ###
tdt <- data.table(tax_table(oslo_fjord),
                 TotalCounts = taxa_sums(oslo_fjord),
                 OTU = taxa_names(oslo_fjord))

taxcumsum <- tdt[, .N, by = TotalCounts]
setkey(taxcumsum, TotalCounts)
taxcumsum[, CumSum := cumsum(N)]

# Cumulation curve
plot_cumsum <- ggplot(taxcumsum, aes(TotalCounts, CumSum)) +
                  geom_point() +
                  theme_bw() +
                  xlab("Abundance per ASV") +
                  ylab("Cumulative number of ASVs") +
                  ggtitle("Cumulation curve")
ggsave(plot = plot_cumsum, filename = file.path(figs_dir, "plot_cumsum.pdf"), device = "pdf", width = 15, height = 15, scale = 1, units = "cm")

# Prevalence plot

## Use own function
# sourceDirectory('/cluster/home/alexaei/ampliconanalysis_package')

sourceDirectory('stats/functions') # load functions in directory

# fast_melt is the functions we are going to use

mdt <- fast_melt(oslo_fjord)
prevdt <- mdt[, list(Prevalence = sum(count > 0),
                    TotalCounts = sum(count)),
                    by = TaxaID]

over10 <- prevdt[prevdt$Prevalence >10,]
max_prev <- prevdt[prevdt$Prevalence == max(prevdt$Prevalence),]
over10

plot_prev <- ggplot(prevdt, aes(Prevalence)) +
                geom_histogram() +
                theme_bw() +
                ylab("Cumulative number of ASVs") +
                ggtitle("Histogram of ASV Prevalence")

ggsave(plot = plot_prev, filename = file.path(figs_dir, "plot_prevalence.pdf"), device = "pdf", width = 15, height = 15, scale = 1, units = "cm")

# Question C: Looking at the cumulation curves and prevalence plots - what are your conclusions with regards to ASV distribution over time?

### alpha diversity ###

# Estimating evenness and richness
# Background. In addition to well-known diversity indices such as Shannon and Simpson, we will compute several widely-used statistical estimators of asymptotic species richness. For example, the chao1 gives the true number of species in the assemblage sampled. These estimators aim to reduce the effect of undersampling, which inevitably biases the observed species count. (more information can be found here https://palaeo-electronica.org/2011_1/238/estimate.htm)

# We need to rarefy the ASV tables first
sample_sums(oslo_fjord)
subset_oslo_fjord <- subset_samples(oslo_fjord, sample_sums(oslo_fjord) > 3000)

minimum_reads <- min(sample_sums(subset_oslo_fjord))
minimum_reads


rarefied_oslo_fjord <- rarefy_even_depth(subset_oslo_fjord, sample.size = minimum_reads,
			rngseed = FALSE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE) #Normalizing species data

# Calculation Shannon, Pilous and ACE and adding to data matrix

richness_estimates <- estimate_richness(rarefied_oslo_fjord)

# Calculate Pielou's evenness.
# Pielou's measure of species evenness, i.e. J = H'/ln(S) where H' is Shannon Weiner diversity and S is the total number of species in a sample, across all samples in dataset.

pielou <- richness_estimates$Shannon/log(specnumber(t(otu_table(rarefied_oslo_fjord))))

# combine matrices
diversity_estimates <- cbind(richness_estimates[,c(-3, -5)], pielou)


# Calculate phylogenetic diversity
# Phylogenetic diversity (“PD”) is a measure of biodiversity, based on phylogeny (the tree of life). Faith (1992) defined the phylogenetic diversity of a set of species as equal to the sum of the lengths of all those branches on the tree that span the members of the set.


# Question D: What are the differences between the various diversity indices? What properties of the community do they determine?

### beta diversity ###
# Background. Ordination techniques such as nMDS are used to study community turnover (beta-diversity). NMDS is particularly useful for visualizing changes in community composition (beta-diversity) since multidimensional data can be ordinated into (plotted in) two-dimensional space. Furthermore, we will add environmental data to our ordination for visual inspection of environmental co-variations with bacterial community composition. This will be combined with inspection of the environmental variable vectors and their significant co-variations with community composition (beta-diversity).

qqnorm(otu_table(subset_oslo_fjord))	#Check distribution of residuals. Distribution is not normal.
qqnorm(decostand(otu_table(subset_oslo_fjord), method="hellinger"))	#Check distribution of residuals. Better, but far from ideal.

#The OTU abundance data is transformed using Hellinger method. The distance matrix is then computed using a coefficient (Bray-Curtis) ignoring double zeroes since they are very numerous, uninformative and so would exaggerate similarity between samples sharing little or no OTUs.

stand_asv <- decostand(otu_table(subset_oslo_fjord), method="hellinger")

#Create the nMDS object.
asv_mds <- metaMDS(stand_asv, distance = "bray", autotransform = FALSE, try = 200)

# Goodness of representation - stress value
# As a rule of thumb literature has identified the following cut-off values for stress-level:
# Higher than 0.2 is poor (risks for false interpretation).
# 0.1 - 0.2 is fair (some distances can be misleading for interpretation).
# 0.05 - 0.1 is good (can be confident in inferences from plot).
# Less than 0.05 is excellent (this can be rare).

asv_mds$stress

gof <- goodness(asv_mds)

plot(asv_mds, display = "sites", type = "n", ylab=c(-4,4))
points(asv_mds, display = "sites", cex = 2*gof/mean(gof))

# Shepards test/goodness of fit
goodness(asv_mds) # Produces a results of test statistics for goodness of fit for each point

stressplot(asv_mds) # Produces a Shepards diagram

asv_mds_3d <- metaMDS(stand_asv, distance = "bray", autotransform = FALSE, try = 200, k = 3)

stress_vec <- numeric(10)
for(i in seq(10)){
  stress_vec[i] <- metaMDS(asv_mds, distance = "bray", autotransform = FALSE, try = 200, k=i)$stress
}

plot(seq(10),stress_vec, type = 'o', ylab = "Stress", xlab = "Number of dimensions",
     col="tomato", pch=19)
abline(h=0.2, lty=2)
abline(h=0.05, lty=3)

# Plot the nMDS object. Is the figure easily readable?
plot(asv_mds)

# Creating the same plot displaying samples only.
ordiplot(asv_mds, display="sites", type="t")


### Fit environmental data ###
# Fit environmental variables to the ordination. Environmental variable values must be scaled and centered since they vary in their physical dimensions. The function envfit does it by default to columns of numeric values present in the dataframe passed to the function.

env_data <- sample_data(subset_oslo_fjord)

env_fit <- envfit(asv_mds, env_data, permu=999, na.rm=T)

# Check significance of environmental vector fitting
env_fit[[1]]

# Plot environmental variable vectors onto the ordination.
plot(env_fit, p.max = 0.05, cex=0.4)

# Check for autocorrelations

### Community analysis - hierarchical clustering ###

# Compute distance matrix for community composition.
asv_dist <- vegdist(stand_asv, method="bray")

# Cluster communities using unweighted pair group method with arithmetic mean (UPGMA) method.
asv_clust <- hclust(asv_dist, method="average")

env_scaled <- as.matrix(scale(env_data, scale=T, center=T))
rownames(env_scaled)= env_data[,1]

# Compute distance matrix for environmental metadata.
env_dist <- vegdist(env_scaled, method="euclidean")

# Cluster environments using UPGMA method.
env_clust <- hclust(env_dist, method="average")

# Convert clustering objects into dendrograms
asv_dendro = as.dendrogram(asv_clust)
env_dendro = as.dendrogram(env_clust)

# Compare dendrograms with a tanglegram.
dl = dendlist(OTU.dendro, env.dendro)

pdf(file = file.path(figs_dir, "dendograms.pdf")
tanglegram(dl, sort = F, common_subtrees_color_lines = F, highlight_distinct_edges = F, highlight_branches_lwd = F, main_left = "Community composition", main_right = "Environmental conditions", common_subtrees_color_branches = FALSE, margin_inner = 5, margin_outer = 5, axes=F)
dev.off()

# Check dendrogram comparison plot in your "M:/pc/Dokumenter/" folder.
# Question D: Are dendrograms more similar if you use the environmental variables that fit most significantly your community composition ordination? Hint: retry the above with env_table variables you consider more relevant.
# Question E: Did you notice anything interesting or striking from the time series? Did you observe distinct groups?

### Unifrac distance ###


