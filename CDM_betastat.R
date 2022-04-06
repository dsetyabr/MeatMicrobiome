###load library
# install.packages(c("vegan","ggplot2","ggpubr","data.table","phyloseq","qiime2R","tidyr","naniar","raster"))
# install.packages("qiime2R")
library(vegan) 
library(ggplot2)
library(ggpubr)
library(data.table)
library(phyloseq)
library(qiime2R)
library(tidyr)
library(naniar)
library(raster)
library(dplyr)
library(tidyverse)
library(fantaxtic)
library(tibble)
library(microbiomeMarker)

###set working directory
#mac
setwd("~/Desktop/microbiomeCCDA")

###set treatment position
positions <- c("INI","WA", "DA", "DWA","UDA")
metadata <- metadata %>% mutate(Treatment=parse_factor(Treatment, levels=positions))

###load data
metadata <- read.csv("trial/alpha/CDM_alpha.csv", na.strings = c("","NA"), header=TRUE)
colnames(metadata)[colnames(metadata) == "sample.id"] = "SampleID"
Wunifrac <- read_qza("trial/beta/weighted_unifrac_pcoa_results.qza")
BC <- read_qza("trial/beta/bray_curtis_pcoa_results.qza")

#weighted unifrac plot
Wunifrac$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color=`Treatment`, shape=`Treatment`)) +
  ggtitle("Weighted Unifrac")+
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_q2r()+
  stat_ellipse()
ggsave("CDM_weightuni.tiff", width = 5, height = 5, path = "Trial/image") 
#Bray Curtis plot
BC$data$Vectors %>%
  dplyr::select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color=`Treatment`, shape=`Treatment`)) +
  ggtitle("Bray-Curtis")+
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_q2r()+
  stat_ellipse()
ggsave("CDM_BC.tiff", width = 5, height = 5, path = "Trial/image")

###Phyloseq 
#table
physeq<-qza_to_phyloseq(
  features="trial/rarefied_table.qza",
  tree="trial/rooted-tree.qza",
  "trial/taxonomy.qza",
  metadata = "trial/metadata.tsv"
)

#plot (using fantaxtic)

physeq1 <- get_top_taxa(physeq_obj = physeq, n = 15, relative = TRUE,
                       discard_other = FALSE, other_label = "Other")
pshyseq1 <- name_taxa(physeq1, label = "Unkown", species = F, other_label = "Other")

fantaxtic_bar(physeq1, color_by = "Genus", label_by = "Genus", other_label = "Other",order_alg = "alph",facet_by = "Treatment")+
  theme(axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size=20))+
  theme(plot.title =  element_text(size=25, face="bold"))+
  theme(legend.title=element_text(size=15),legend.text = element_text(size=15))
ggsave("CDM_Genus1.tiff", width = 15, height = 10, path = "Trial/image")

fantaxtic_bar(physeq1, color_by = "Order", label_by = "Order", other_label = "Other",order_alg = "alph",facet_by = "Treatment")+
  theme(axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size=20))+
  theme(plot.title =  element_text(size=25, face="bold"))+
  theme(legend.title=element_text(size=15),legend.text = element_text(size=15))
ggsave("CDM_Order1.tiff", width = 15, height = 10, path = "Trial/image")

###LEFSE (using microbiomeMarker)
#LEfSe (Linear discriminant analysis Effect Size) determines the features 
#(organisms, clades, operational taxonomic units, genes, or functions) 
#most likely to explain differences between classes by coupling standard 
#tests for statistical significance with additional tests encoding 
#biological consistency and effect relevance.

lefse <-run_lefse(physeq, group="Treatment", wilcoxon_cutoff = 0.05,taxa_rank = "Genus",kw_cutoff = 0.05,
                  multigrp_strat = TRUE,
                  lda_cutoff = 4)
head(marker_table(lefse))
marker <- (marker_table(lefse))

write.csv(marker,"~/Desktop/microbiomeCCDA/trial/beta/CDM_lefse.csv")

###ANCOM (using microbiomeMarker)
#Determine taxa whose absolute abundances, per unit volume, of the ecosystem (e.g. gut) 
#are significantly different with changes in the covariate of interest (e.g. the group effect).

ancom <-run_ancom(physeq,group = "Treatment",p_adjust = "fdr", taxa_rank = "Genus")
markerancom <- (marker_table(ancom))

write.csv(markerancom,"~/Desktop/microbiomeCCDA/trial/beta/CDM_ancom.csv")

###ANCOMBC (using microbiomeMarker)
ancombc <-run_ancombc(physeq,group = "Treatment",formula = "Treatment", taxa_rank = "Genus", p_adjust = "fdr")
markerancombc <- (marker_table(ancombc))

write.csv(markerancombc,"~/Desktop/microbiomeCCDA/trial/beta/CDM_ancombc.csv")

###PERMANOVA
#BC
dist.bray <- phyloseq::distance(physeq, method = "bray")
dist.bray <- as.dist(dist.bray)

permanovabc <- adonis(dist.bray ~ metadata$Treatment, permutations = 999)
permanovabc

#checking pairwise permanova
OTU1 = as(otu_table(physeq), "matrix")
if(taxa_are_rows(physeq1)){OTU1 <- t(OTU1)}
OTUdf = as.data.frame(OTU1)

PairBC <- pairwise.adonis(OTUdf,metadata$Treatment)

ps.disper <- betadisper(dist.bray, metadata$Treatment)
permutest(ps.disper, pairwise = TRUE)


#WeightUNI
OTU_Unifrac <- UniFrac(physeq, weighted = TRUE)

permanovaWuni <- adonis(OTU_Unifrac~ metadata$Treatment, permutations = 999)
permanovaWuniphy

BRD_Wu <- betadisper(OTU_Unifrac, type = c("centroid"), group = metadata$Treatment)
BRD_Wu
boxplot(BRD_Wu)
TukeyHSD(BRD_Wu) ## calculate the distance from the centroids to one group to another
plot(BRD_Wu)

ps.disper <- betadisper(OTU_Unifrac, metadata$Treatment)
permutest(ps.disper, pairwise = TRUE)

###PAIRWISE FUNCTION
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 

