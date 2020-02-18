### ggplot
# References and flWLlates for this exercise can be found here http://r-statistics.co/Complete-Ggplot2-Tutorial-Part1-With-R-Code.html and http://r-statistics.co/ggplot2-Tutorial-With-R.html


### Understanding the Ggplot Syntax ###
# The syntax for constructing ggplots could be puzzling if you are a beginner or work primarily with base graphics. The main difference is that, unlike base graphics, ggplot works with dataframes and not individual vectors. All the data needed to make the plot is typically be contained within the dataframe supplied to the ggplot() itself or can be supplied to respective geoms. More on that later.

# The second noticeable feature is that you can keep enhancing the plot by adding more layers (and themes) to an existing plot created using the ggplot() function.
# Let’s initialize a basic ggplot based on the midwest dataset.

# set up
library("ggplot2")
library("phyloseq")

stats_dir <- "stats/"
figs_dir <- "figs/"
dir.create(stats_dir)
dir.create(figs_dir)

CTD_data <- read.table("CTD_data.tsv", header=T, check.names=F)	#flWL is fluorescence wet lab
oslo_fjord <- readRDS(str_c(dada2_dir, "OsloFjord_phyloseq.rds"))
metadata <- sample_data(oslo_fjord)
metadata$flWL <- CTD_data$flWL
metadata$flWL <- CTD_data$flWL
metadata$DoY <- CTD_data$DoY
metadata$Year <- CTD_data$Year
colnames(metadata)

# Init ggplot
ggplot(metadata, aes(x = Date, y = flWL))  # Date and flWL are columns in 'metadata'

# A blank ggplot is drawn. Even though the x and y are specified, there are no points or lines in it. This is because, ggplot doesn’t assume that you meant a scatterplot or a line chart to be drawn. I have only told ggplot what dataset to use and what columns should be used for X and Y axis. I haven’t explicitly asked it to draw any points.

# Also note that aes() function is used to specify the X and Y axes. That’s because, any information that is part of the source dataframe has to be specified inside the aes() function.

### How to Make a Simple Scatterplot ###
# Let’s make a scatterplot on top of the blank ggplot by adding points using a geom layer called geom_point.

ggplot(metadata, aes(x = Date, y = flWL)) + geom_point()

# We got a basic scatterplot, where each point represents a county. However, it lacks some basic components such as the plot title, meaningful axis labels etc. Moreover most of the points are concentrated on the bottom portion of the plot, which is not so nice. You will see how to rectify these in upcoming steps.

# Like geom_point(), there are many such geom layers which we will see in a subsequent part in this tutorial series. For now, let’s just add a smoothing layer using geom_smooth(method='lm'). Since the method is set as lm (short for linear model), it draws the line of best fit.

g <- ggplot(metadata, aes(x = Date, y = flWL)) + geom_point() + geom_smooth(method="lm")  # set se=FALSE to turnoff confidence bands

plot(g)

# The line of best fit is in blue.
# Question A: Can you find out what other method options are available for geom_smooth? (note: see ?geom_smooth). You might have noticed that majority of points lie in the bottom of the chart which doesn’t really look nice. So, let’s change the Y-axis limits to focus on the lower half.

### Adjusting the X and Y axis limits ###
# The X and Y axis limits can be controlled in 2 ways.

## Method 1: By deleting the points outside the range
# This will change the lines of best fit or smoothing lines as compared to the original data.

# This can be done by xlim() and ylim(). You can pass a numeric vector of length 2 (with max and min values) or just the max and min values itself.

# Delete the points outside the limits
#g + xlim(c(42400, 43100)) + ylim(c(0, 15))   # deletes points

# In this case, the chart was not built from scratch but rather was built on top of g. This is because, the previous plot was stored as g, a ggplot object, which when called will reproduce the original plot. Using ggplot, you can add more layers, themes and other settings on top of this plot.

# Did you notice that the line of best fit became more horizontal compared to the original plot? This is because, when using xlim() and ylim(), the points outside the specified range are deleted and will not be considered while drawing the line of best fit (using geom_smooth(method='lm')). This feature might come in handy when you wish to know how the line of best fit would change when some extreme values (or outliers) are removed.

## Method 2: Zooming In
# The other method is to change the X and Y axis limits by zooming in to the region of interest without deleting the points. This is done using coord_cartesian().

# Let’s store this plot as g1.
g1 <- g + coord_cartesian(xlim=c(42400, 43100), ylim=c(0, 15))  # zooms in
plot(g1)

# Since all points were considered, the line of best fit did not change.

### How to Change the Title and Axis Labels ###
# Let’s add the plot title and labels for X and Y axis. This can be done in one go using the labs() function with title, x and y arguments. Another option is to use the ggtitle(), xlab() and ylab().

# Add Title and Labels
g1 + labs(title="Water fluorescence (wet lab) in Oslo Fjord", subtitle="in 2016 and 2017", y="Fluorescence (wet lab)", x="Day of the Year", caption="CTD")

# or

g1 + ggtitle("Water fluorescence (wet lab) in Oslo Fjord", subtitle="in 2016 and 2017") + xlab("Date") + ylab("Fluorescence (wet lab)")

### How to Change the Color and Size of Points ###
## How to Change the Color and Size To Static?
# We can change the aesthetics of a geom layer by modifying the respective geoms. Let’s change the color of the points and the line to a static value.

ggplot(metadata, aes(x = DoY, y = flWL)) +
  geom_point(col="steelblue", size=3) +   # Set static color and size for points
  geom_smooth(method="lm", col="firebrick") +  # change the color of line
  coord_cartesian(xlim=c(0, 365), ylim=c(0, 1)) +
  labs(title="Water fluorescence (wet lab) in Oslo Fjord", subtitle="in 2016 and 2017", y="Fluorescence (wet lab)", x="Day of the Year", caption="CTD")

## How to Change the Color To Reflect Categories in Another Column?
# Suppose if we want the color to change based on another column in the source dataset (metadata), it must be specified inside the aes() function.

gg <- ggplot(metadata, aes(x = DoY, y = flWL)) +
  geom_point(aes(col=Year), size=3) +  # Set color to vary based on state categories.
  geom_smooth(method="lm", col="firebrick", size=2) +
  coord_cartesian(xlim=c(0, 365), ylim=c(0, 1)) +
  labs(title="Water fluorescence (wet lab) in Oslo Fjord", subtitle="in 2016 and 2017", y="Fluorescence (wet lab)", x="Day of the Year", caption="CTD")
plot(gg)

# Now each point is colored based on the state it belongs because of aes(col=state). Not just color, but size, shape, stroke (thickness of boundary) and fill (fill color) can be used to discriminate groupings.

# As an added benefit, the legend is added automatically. If needed, it can be removed by setting the legend.position to None from within a theme() function.

gg + theme(legend.position="None")  # remove legend

# Also, You can change the color palette entirely.

#gg + scale_colour_brewer(palette = "Set1")  # change color palette

# More of such palettes can be found in the RColorBrewer package

library(RColorBrewer)
head(brewer.pal.info, 10)

### How to Change the X Axis Texts and Ticks Location ###
## How to Change the X and Y Axis Text and its Location?
# Alright, now let’s see how to change the X and Y axis text and its location. This involves two aspects: breaks and labels.

# Step 1: Set the breaks
# The breaks should be of the same scale as the X axis variable. Note that I am using scale_x_continuous because, the X axis variable is a continuous variable. Had it been a date variable, scale_x_date could be used. Like scale_x_continuous() an equivalent scale_y_continuous() is available for Y axis.

# Change breaks
gg + scale_x_continuous(breaks=seq(0, 365, 30))

# Step 2: Change the labels You can optionally change the labels at the axis ticks. labels take a vector of the same length as breaks.

gg + scale_x_continuous(breaks=seq(0, 365, 30), labels = letters[1:13])

# If you need to reverse the scale, use scale_x_reverse().
# Reverse X Axis Scale
gg + scale_x_reverse()

## How to Customize the Entire Theme in One Shot using Pre-Built Themes?
# Finally, instead of changing the theme components individually (which I discuss in detail in part 2), we can change the entire theme itself using pre-built themes. The help page ?theme_bw shows all the available built-in themes.

# This again is commonly done in couple of ways. * Use the theme_set() to set the theme before drawing the ggplot. Note that this setting will affect all future plots. * Draw the ggplot and then add the overall theme setting (eg. theme_bw())

gg <- gg + scale_x_continuous(breaks=seq(0, 365, 30))

# method 1: Using theme_set()
theme_set(theme_classic())  # not run
gg

# method 2: Adding theme Layer itself.
gg + theme_bw() + labs(subtitle="BW Theme")
gg + theme_classic() + labs(subtitle="Classic Theme")

# For more customized and fancy themes have a look at the ggthemes package and the ggthemr package.

### Bar plots - plotting taxonomy ###
# plot divisions
color_palette = c("#D1BBD7", "#AE76A3", "#882E72", "#1965B0",
                        "#5289C7", "#7BAFDE", "#4EB265", "#90C987",
                        "#CAE0AB", "#F7F056", "#F6C141", "#F1932D",
                        "#E8601C", "#DC050C", "#72190E")

taxonomy = "Division"

classGlom = tax_glom(oslo_fjord, taxrank = taxonomy)

taxon_table = otu_table(classGlom)
tax_matrix = as(tax_table(classGlom), 'matrix')
rownames(taxon_table) = tax_matrix[,taxonomy]
tax_table = prop.table(taxon_table, margin = 2)*100
tax_table <- tax_table[order(rowSums(-taxon_table)),]
rowSums(taxon_table)
dim(taxon_table)
tax_table <- tax_table[1:15,]
rownames(tax_table)

# calculate proportions based on phylum data
library(reshape)
melted_tax_table = melt(tax_table)
colnames(melted_tax_table) <- c("Var1", "Var2", "value")
melted_tax_table <- arrange(melted_tax_table, Var1, desc(value))
melted_tax_table$Var2 <- factor(melted_tax_table$Var2 , levels = unique(melted_tax_table$Var2))

division_colors <- setNames(color_palette, levels(melted_tax_table$Var1))

pdf(file.path(figs_dir, "Barplot_divisions.pdf"), width=120/25.4, height=160/25.4, pointsize = 6)
  ggplot(melted_tax_table, aes(x = Var2, y = value, fill = Var1))+
     geom_bar(stat = "identity", position = "stack", alpha = .5) +
     guides(fill = guide_legend(title = taxonomy)) +
     coord_flip() +
     theme(axis.text = element_text(size=6),
       axis.title = element_text(size=10, face="bold"),
       legend.text = element_text(size=8),
       plot.background = element_rect(fill = "white"),
       panel.background = element_rect(fill = "white"),
       axis.line.x = element_line(color = "grey")) +
     xlab("Sample") +
     ylab("Standardized number of reads") +
     scale_fill_manual(values = division_colors) +
     scale_x_discrete(limits = rev(levels(melted_tax_table$Var2))
     )
dev.off()

### Making a nice RDA plot ###

metadata1 = data.frame(sample_data(classGlom))
metadata1 = metadata1[,c()]

library(vegan)
rda_div_meta <- rda(decostand(t(tax_table), method="hellinger"), as.matrix(metadata1))

rda_div_meta
summary_rda_div_meta <- summary(rda_div_meta)
RsquareAdj(rda_div_meta)

library(ggfortify)
#p_rda <- autoplot(rda_div_meta, layers = c("species", "biplot"), legend.position = "none")
p_rda <- biplot(rda_div_meta, scaling = "species", display = c("sites", "species"))

rda_scores_env = vegan::scores(rda_div_meta, display = "sites")
rda_scores_species = vegan::scores(rda_div_meta, display = "species")

names = rownames(rda_scores_species)
names
melted_tax_table$Var1

rda_plot_species <- ggplot(data.frame(rda_scores_species), aes(x = PC1, y = PC2, color = names)) +
                geom_point(size = 1, alpha = .5) +
                scale_color_manual(values = division_colors) +
				xlab("RDA1") +
				ylab("RDA2")

mult = 0.2

rda_biplot_division <- rda_plot_species +
  geom_segment(data = data.frame(rda_scores_env), aes(x = 0, xend = mult * PC1,
                    y = 0, yend = mult * PC2),
                arrow = arrow(length = unit(0.125, "cm")), colour = "red", alpha = .4, size = 0.1) +
  geom_text(data = data.frame(rda_scores_env),
            aes(x = mult * PC1, y = mult * PC2, label = rownames(rda_scores_env),
                hjust = 0.5 * (1-sign(PC1)), vjust = 0.5 * (1-sign(PC2))),
                color = "red", size = 1.5, alpha = .4) +
  coord_cartesian(xlim = c(-0.25, 0.52), ylim = c(-0.22, 0.2)) +
  xlab("RDA1") +
  ylab("RDA2") +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(color = "grey"),
		legend.position = "bottom", 
		legend.text=element_text(size=4, 
		margin = margin(t = 0)), 
		legend.spacing.x = unit(0, 'cm'), 
		legend.spacing.y = unit(0, 'cm'), 
		legend.key.height=unit(0.1, 'lines'),
		panel.grid.minor = element_blank(),
		panel.grid.major = element_blank(),
		legend.title = element_blank())

pdf(file.path(figs_dir, "rda_divisions.pdf"), width=80/25.4, height=80/25.4, pointsize = 6)
rda_biplot_division
dev.off()
