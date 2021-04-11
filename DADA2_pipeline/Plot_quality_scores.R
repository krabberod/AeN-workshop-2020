

plot_list = list()
for (i in 1:length(fns)) {
  plot_list[[i]]<- plotQualityProfile(fns[i])
}
ml<- marrangeGrob(plot_list, nrow=4, ncol=4) # skriver til skjerm 4x4 plot

ml # for å printe på skjermen
ggsave("multipage.pdf", ml , scale = 1, width = 200,  height = 200, units = "mm") # skriver til pdf

plotQualityProfile(fns_R1, aggregate = T)
plotQualityProfile(fns_R2, aggregate = T)
plotQualityProfile(filt_R1, aggregate = T)
plotQualityProfile(filt_R2, aggregate = T)
