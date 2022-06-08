#diversity analysis with phyloseq
#lets try with your data
rm(list=ls()) #Eliminar todas las listas

#ubicar directorio
setwd("C:/Users/Intel Core i5/Desktop/TXS Y MICROBIOTA/AN?LISIS BIOINFORMATICO/diversity phyloseq/Phyloseq")
getwd()

#import data from otu table
otu.out <- as.matrix(read.csv("otu_table.txt", header = T, row.names = NULL, sep = "\t", check.names=F))
str(otu.out)
#generar la transpuesta
otu1 <- t(otu.out)
class(otu1)
#assign colnames
colnames(otu1) <-otu1[1,]
#remove duplicate rows
otu2<- otu1[-1,]
class(otu2)

#convertir a numerica
ncol(otu2)
otu5 <- matrix(
  as.numeric(otu2), ncol = ncol(otu2)) 

#reassign rownames and colnames
rownames(otu5) <- rownames(otu2)
colnames(otu5) <- colnames(otu2)

class(otu5)
str(otu5)

#import data from taxonomy
taxa.out <- as.matrix(read.table("taxa_table.txt", header = T, row.names = NULL, sep = "\t", check.names=F))
str(taxa.out)

#assign rownames
rownames(taxa.out) <- taxa.out[,1]

#remove duplicate columns
taxa<- taxa.out[,-1]
class(taxa)

#prepare to make phyloseq object
OTU = otu_table(otu5, taxa_are_rows = TRUE)
TAX = tax_table(taxa)
OTU
TAX

physeq = phyloseq(OTU, TAX)
physeq

#test_physeq object
plot_bar(physeq, fill = "Family")

#add sample data
#import data from taxonomy
sample.out <- as.data.frame(read.table("sample_data.txt", header = T, row.names = 1, sep = "\t", check.names=F))
str(sample.out)

#prepare a phyloseq_object
samplingdata <- sample_data(sample.out)

#Agregar random_tree

library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)

#physeq1

physeq1 = merge_phyloseq(physeq, samplingdata, random_tree)
physeq1

#prepare new physeq object with sample data

##physeq2 = phyloseq(OTU, TAX, samplingdata)
##physeq2

physeq2 = phyloseq(OTU, TAX, samplingdata, random_tree)
physeq2

#Intall rlang.

install.packages("rlang")
install.packages("dplyr")
library(rlang)
library(dplyr)

#based on physeq2 object
library("tidyverse")
library("microbiome")
library("knitr")
library("phyloseq")

#data transformations 
physeq4 <- microbiome::transform(physeq2, "log10p")
physeq4 = filter_taxa(physeq4, function(x) sum(x > 3) > (0.2*length(x)), TRUE)


#normalize
#Normalize number of reads in each sample using median sequencing depth.

total = median(sample_sums(physeq2))
standf = function(x, t=total) round(t * (x / sum(x)))
physeq3 = transform_sample_counts(physeq2, standf)

physeq3 = filter_taxa(physeq3, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

GP = physeq4

ordu = ordinate(GP, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(GP, ordu, color="Treatment")

#devtools::install_github("Artjom-Metro/ggsignif")
#https://stackoverflow.com/questions/29263046/how-to-draw-the-boxplot-with-significant-level
#https://github.com/const-ae/ggsignif/issues/8
#https://rdrr.io/cran/ggsignif/



#alpha diversity
library(ggsignif)
p1 <- plot_richness(physeq3, x="Treatment", measures = "Fisher")

p1 <- plot_richness(physeq3, x="Treatment")

p2 <- p1 + geom_boxplot() 

p3 <- p2 + geom_signif(comparisons = list(c("Basal", "Metformin"), c("Treatment", "Linagliptin/metformin"), c("Metformin", "Linagliptin/metformin")), 
              map_signif_level=TRUE, textsize = 2, step_increase = .1) 

#save image
ggsave("alpha_index_Fisher_status1.png",
       plot = p3,
       device = "png",
       width = 6,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)


#new transformation
physeq5 <- microbiome::transform(ps.gen, "log10p")
physeq5 = filter_taxa(ps.gen, function(x) sum(x > 3) > (0.2*length(x)), TRUE)


GP=physeq5

GP.ord <- ordinate(GP, "NMDS", "jaccard")
v1=sample_data(GP)$status
q1 = plot_ordination(GP, GP.ord, type="taxa", color="Phylum", title="taxa")
print(q1)

f1 <- q1 + stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  theme_bw()

#save image
ggsave("beta_index_NMDS_phylum_Treatment.png",
       plot = f1,
       device = "png",
       width = 6,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)


sample_data(GP)$status


f1 + facet_wrap(~Phylum, 5)


q2 = plot_ordination(GP, GP.ord, type="samples", color="Treatment") 
f3 = q2 + stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  theme_bw()

#save image
ggsave("beta_index_NMDS_Treatment1.png",
       plot = f3,
       device = "png",
       width = 6,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

f4 = q2 + facet_wrap(~status, 5)

#q2 + geom_polygon(aes(fill=Treatment)) + geom_point(size=5) + ggtitle("samples")


q4 = plot_ordination(GP, GP.ord, type="split", color="Phylum", label="Treatment", title="split") 
q4


GP = physeq3
wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 5), A=0.5*nsamples(GP))
GP1 = prune_taxa(wh0, GP)
GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))
phylum.sum = tapply(taxa_sums(GP1), tax_table(GP1)[, "Family"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
GP1 = prune_taxa((tax_table(GP1)[, "Family"] %in% top5phyla), GP1)


#human = get_variable(GP1, sam_data) %in% c("Type")
#sample_data(GP1)$human <- factor(human)


#esta funcion obtiene la grafica de nmbds puede acomodarse ya sea phylum o lo que sea

GP.ord <- ordinate(GP, "NMDS", "bray")
p4 = plot_ordination(GP, GP.ord, type="taxa", color="Family", title="taxa")
print(p4)

#p1 +  stat_ellipse(type = "norm", linetype = 2) +
#stat_ellipse(type = "t") +
#theme_bw()



p5 = plot_ordination(GP1, GP.ord, type="samples", color="Treatment")#,shape="human") 

#p1 + facet_wrap(~ta2, 3)


nmds_c_vs_c <- p2 +  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  theme_bw()  + scale_color_manual(name = "Tipo", labels = c("ST", "M", "LM"), values = c("blue", "red","orange"))

#p2 + geom_polygon(aes(fill=Type)) + geom_point(size=5) + ggtitle("samples")


ggsave(plot = nmds_c_vs_c,filename =  "nmds_status1.pdf", width=10, height=5)

#########################################
##### for Groups##########################
#########################################

#physeq1

physeq1 = merge_phyloseq(physeq, samplingdata, random_tree)
physeq1

#prepare new physeq object with sample data

##physeq2 = phyloseq(OTU, TAX, samplingdata)
##physeq2

physeq2 = phyloseq(OTU, TAX, samplingdata, random_tree)
physeq2

#Intall rlang.

install.packages("rlang")
install.packages("dplyr")
library(rlang)
library(dplyr)

#based on physeq2 object
library("tidyverse")
library("microbiome")
library("knitr")
library("phyloseq")

#data transformations 
physeq4 <- microbiome::transform(physeq2, "log10p")
physeq4 = filter_taxa(physeq4, function(x) sum(x > 3) > (0.2*length(x)), TRUE)


#normalize
#Normalize number of reads in each sample using median sequencing depth.

total = median(sample_sums(physeq2))
standf = function(x, t=total) round(t * (x / sum(x)))
physeq3 = transform_sample_counts(physeq2, standf)

physeq3 = filter_taxa(physeq3, function(x) sum(x > 3) > (0.2*length(x)), TRUE)

GP = physeq4

ordu = ordinate(GP, "PCoA", "unifrac", weighted=TRUE)
plot_ordination(GP, ordu, color="groups")

#devtools::install_github("Artjom-Metro/ggsignif")
#https://stackoverflow.com/questions/29263046/how-to-draw-the-boxplot-with-significant-level
#https://github.com/const-ae/ggsignif/issues/8
#https://rdrr.io/cran/ggsignif/



#alpha diversity
library(ggsignif)
p1 <- plot_richness(physeq3, x="Groups", measures = "Fisher")

p1 <- plot_richness(physeq3, x="Groups")

p2 <- p1 + geom_boxplot() 

p3 <- p2 + geom_signif(comparisons = list(c("Basal", "Metformin_6"), c("Basal", "Metformin_12"), c("Basal", "Linagliptin/metformin_6"), c("Basal", "Linagliptin/metformin_12"), c("Metformin_6", "Linagliptin/metformin_6"), c("Metformin_12", "Linagliptin/metformin_12"), c("Metformin_6", "Metformin_12"), c("Linagliptin/metformin_6", "Linagliptin/metformin_12")), 
                       map_signif_level=TRUE, textsize = 2, step_increase = .1) 

#save image
ggsave("alpha_index_Fisher_groups.png",
       plot = p3,
       device = "png",
       width = 6,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)


#new transformation
physeq5 <- microbiome::transform(physeq2, "log10p")
physeq5 = filter_taxa(physeq2, function(x) sum(x > 3) > (0.2*length(x)), TRUE)


GP=physeq5

GP.ord <- ordinate(GP, "NMDS", "jaccard")
v1=sample_data(GP)$status
q1 = plot_ordination(GP, GP.ord, type="taxa", color="Phylum", title="taxa")
print(q1)

f1 <- q1 + stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  theme_bw()

#save image
ggsave("beta_index_NMDS_phylum_groups.png",
       plot = f1,
       device = "png",
       width = 6,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)


sample_data(GP)$status


f1 + facet_wrap(~Phylum, 5)


q2 = plot_ordination(GP, GP.ord, type="samples", color="Groups") 
f3 = q2 + stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  theme_bw()

#save image
ggsave("beta_index_NMDS_groups.png",
       plot = f3,
       device = "png",
       width = 6,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)

f4 = q2 + facet_wrap(~groups, 5)

#q2 + geom_polygon(aes(fill=status)) + geom_point(size=5) + ggtitle("samples")


q4 = plot_ordination(GP, GP.ord, type="split", color="Phylum", label="Groups", title="split") 
q4


GP = physeq3
wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 5), A=0.5*nsamples(GP))
GP1 = prune_taxa(wh0, GP)
GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))
phylum.sum = tapply(taxa_sums(GP1), tax_table(GP1)[, "Family"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
GP1 = prune_taxa((tax_table(GP1)[, "Family"] %in% top5phyla), GP1)


#human = get_variable(GP1, sam_data) %in% c("Type")
#sample_data(GP1)$human <- factor(human)


#esta funcion obtiene la grafica de nmbds puede acomodarse ya sea phylum o lo que sea

GP.ord <- ordinate(GP, "NMDS", "bray")
p4 = plot_ordination(GP, GP.ord, type="taxa", color="Family", title="taxa")
print(p4)

#p1 +  stat_ellipse(type = "norm", linetype = 2) +
#stat_ellipse(type = "t") +
#theme_bw()



p5 = plot_ordination(GP1, GP.ord, type="samples", color="Groups")#,shape="human") 

#p1 + facet_wrap(~ta2, 3)


nmds_c_vs_c <- p2 +  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  theme_bw()  + scale_color_manual(name = "Tipo", labels = c("Basal", "Metformin_6", "Metformin_12", "Linagliptin/metformin_6", "Linagliptin/metformin_12"), values = c("blue", "red","orange"))

#p2 + geom_polygon(aes(fill=Type)) + geom_point(size=5) + ggtitle("samples")


ggsave(plot = nmds_c_vs_c,filename =  "nmds_groups.pdf", width=10, height=5)

#alpha diversity
library(ggsignif)
p1 <- plot_richness(physeq3, x="Month", measures = "Fisher")

p1 <- plot_richness(physeq3, x="Month")

p2 <- p1 + geom_boxplot() 

p3 <- p2 + geom_signif(comparisons = list(c("0_month", "6_month"), c("0_month", "12_month"), c("6_month", "12_month")), 
                       map_signif_level=TRUE, textsize = 2, step_increase = .1) 

#save image
ggsave("alpha_index_Fisher_month.png",
       plot = p3,
       device = "png",
       width = 6,
       height = 5,
       units = c("in", "cm", "mm", "px"),
       dpi = 300)
