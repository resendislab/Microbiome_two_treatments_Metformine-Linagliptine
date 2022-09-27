# Diversity analysis with phyloseq

rm(list=ls())

library("ape")
library("rlang")
library("dplyr")
library("tidyverse")
library("microbiome")
library("knitr")
library("phyloseq")
library("ggsignif")

#####

# Import data from OTU table
otu_out <- as.matrix(read.csv("otu_table.txt", header = T, row.names = NULL, sep = "\t", check.names=F))
pre_otu <- t(otu_out)
colnames(pre_otu) <- pre_otu[1,]
pre_otu <- pre_otu[-1,]
otu <- matrix(as.numeric(pre_otu), ncol = ncol(pre_otu)) 
rownames(otu) <- rownames(pre_otu)
colnames(otu) <- colnames(pre_otu)

# Import data from taxonomy table
taxa_out <- as.matrix(read.table("taxa_table.txt", header = T, row.names = NULL, sep = "\t", check.names=F))
rownames(taxa_out) <- taxa_out[,1]
taxa <- taxa_out[,-1]

# Add sample data
sample_out <- as.data.frame(read.table("sample_data.txt", header = T, row.names = 1, sep = "\t", check.names=F))

# Built the phyloseq object
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(taxa)
samplingdata <- sample_data(sample_out)
physeq = phyloseq(OTU, TAX, samplingdata)

#####

# Normalize number of reads in each sample using median sequencing depth
total = median(sample_sums(physeq))
standf = function(x, t=total) round(t * (x / sum(x)))
physeq2 = transform_sample_counts(physeq, standf)
physeq2 = filter_taxa(physeq2, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

# Alpha diversity
p <- plot_richness(physeq2, x="status")
p <- p + geom_boxplot() 
p <- p + geom_signif(comparisons = list(c("sin tratamiento", "metformina"), c("sin tratamiento", "linagliptina/metformina"), c("metformina", "linagliptina/metformina")), 
                       map_signif_level=TRUE, textsize = 2, step_increase = .1) 

# Save image
ggsave("alpha_index_Fisher_status.png",
       plot = p,
       device = "png",
       width = 6,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

#####

q <- plot_richness(physeq2, x="groups")
q <- q + geom_boxplot() 
q <- q + geom_signif(comparisons = list(c("sin tratamiento", "metformina_6"), c("sin tratamiento", "metformina_12"), c("sin tratamiento", "linagliptina/metformina_6"), c("Sin tratamiento", "linagliptina/metformina_12"), c("metformina_6", "linagliptina/metformina_6"), c("metformina_12", "linagliptina/metformina_12"), c("metformina_6", "metformina_12"), c("linagliptina/metformina_6", "linagliptina/metformina_12")), 
                       map_signif_level=TRUE, textsize = 2, step_increase = .1) 

#save image
ggsave("alpha_index_Fisher_groups.png",
       plot = q,
       device = "png",
       width = 6,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

#####

physeq3 <- microbiome::transform(physeq, "log10p")
physeq3 = filter_taxa(physeq, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

GP=physeq3
GP.ord <- ordinate(GP, "NMDS", "jaccard")

r = plot_ordination(GP, GP.ord, type="samples", color="groups") 
r = r + stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  theme_bw()

# Save image
ggsave("beta_index_NMDS_groups.png",
       plot = r,
       device = "png",
       width = 6,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)