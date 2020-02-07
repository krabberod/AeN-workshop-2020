### Phyloseq

oslo_fjord<-ps_dada2


sample_names(oslo_fjord)
rank_names(oslo_fjord)
sample_variables(oslo_fjord)

### Normalization using median sequencing depth
total = median(sample_sums(oslo_fjord))
standf = function(x, t=total) round(t * (x / sum(x)))
oslo_fjord = transform_sample_counts(oslo_fjord, standf)

### Rarefaction curves
rarecurve(OTU.table, step = 100, xlab = "Sample Size", ylab = "OTUs", label = TRUE)

### alpha diversity
alpha.div = diversity(OTU.table, index = "shannon", MARGIN = 1, base = exp(1))
alpha.div
#Calculate Fisher's alpha diversity.
alpha.div$fisher = fisher.alpha(OTU.table, MARGIN = 1, se = FALSE)
#Calculate Pielou's evenness.
alpha.div$pielou = alpha.div/log(specnumber(OTU.table))
