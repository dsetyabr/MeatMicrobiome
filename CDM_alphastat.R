###load library
#install.packages(c("jtools","afex","emmeans","lubridate","jtools","sjstats"))
library(afex)
library(lme4)
library(emmeans)
library(lubridate)
library(ggplot2)
library("cowplot")
theme_set(theme_grey())
library(jtools)
library(ggpubr)
library(sjstats)
rm(list = ls ())
library(dplyr)
library(tidyverse)

###set working directory
#mac
setwd("~/Desktop/microbiomeCCDA")

###load data
metadata <- read.csv("trial/alpha/CDM_alpha.csv", na.strings = c("","NA"), header=TRUE)

#set treatment position
positions <- c("INI","WA", "DA", "DWA","UDA")
metadata <- metadata %>% mutate(Treatment=parse_factor(Treatment, levels=positions))

###STAT ANALYSIS
#1. Observed OTUs
#2. Chao1 (measures richness of the environment)
#3. Pielou_e (measures evenness)
#4. Faith_pd (phylogenetic diversity)

set_sum_contrasts() # important for afex

###INTERACTION WITH INI
M1<- mixed(observed_features ~ Treatment*Source + (1|Animal), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1)
stat<-as.data.frame(M1[[1]])
rownames(stat)[rownames(stat) == c("Treatment","Source","Treatment:Source")] = 
  c("ObsOtu-Treatment", "ObsOtu-Source","ObsOtu-Treatment:Source")
allstat <- stat

ls<-lsmeans(M1, pairwise~Treatment*Source, adjust="tukey")
lsdf<-as.data.frame(ls[[1]])
lsdf$exp <- "ObsOtu"
allls <- lsdf

M2 <- mixed(pielou_evenness ~ Treatment*Source + (1|Animal), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2)
stat<-as.data.frame(M2[[1]])
rownames(stat)[rownames(stat) == c("Treatment","Source","Treatment:Source")] = 
  c("Pielou-Treatment", "Pielou-Source","Pielou-Treatment:Source")
allstat <- rbind(allstat,stat)

ls<-lsmeans(M1, pairwise~Treatment*Source, adjust="tukey")
lsdf<-as.data.frame(ls[[1]])
lsdf$exp <- "pielou"
allls <- rbind(allls,lsdf)

M3 <- mixed(chao1 ~ Treatment*Source + (1|Animal), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3)
stat<-as.data.frame(M3[[1]])
rownames(stat)[rownames(stat) == c("Treatment","Source","Treatment:Source")] = 
  c("Chao-Treatment", "Chao-Source","Chao-Treatment:Source")
allstat <- rbind(allstat,stat)

ls<-lsmeans(M1, pairwise~Treatment*Source, adjust="tukey")
lsdf<-as.data.frame(ls[[1]])
lsdf$exp <- "chao"
allls <- rbind(allls,lsdf)

M4 <- mixed(faith_pd ~ Treatment*Source + (1|Animal), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4)
stat<-as.data.frame(M4[[1]])
rownames(stat)[rownames(stat) == c("Treatment","Source","Treatment:Source")] = 
  c("Faith-Treatment", "Faith-Source","Faith-Treatment:Source")
allstat <- rbind(allstat,stat)

ls<-lsmeans(M1, pairwise~Treatment*Source, adjust="tukey")
lsdf<-as.data.frame(ls[[1]])
lsdf$exp <- "faith"
allls <- rbind(allls,lsdf)

#save ls and stat
write.csv(allstat,"~/Desktop/microbiomeCCDA/trial/alpha/CDM_alpha_stat.csv")
write.csv(allls,"~/Desktop/microbiomeCCDA/trial/alpha/CDM_alpha_lsmeans.csv")

#checking assumptions
plot(M1$full_model)
plot(M2$full_model)
plot(M3$full_model)
plot(M4$full_model)

#check normality
qqnorm(residuals(M1$full_model))
qqnorm(residuals(M2$full_model))
qqnorm(residuals(M3$full_model))
qqnorm(residuals(M4$full_model))

###PLOT
a <- afex_plot(M1, x = "Treatment", trace="Source", id = "Animal", dodge = 0.4, point_arg = list(size = 4),  line_arg = list(linetype = 0), mapping="color") +
  theme_bw() + theme(legend.position="bottom") +
  labs(y = "Observed ASVs", x = "Treatment") +
  theme(legend.position="none") 

c <- afex_plot(M2, x = "Treatment",trace = "Source", id = "Animal", dodge = 0.4, point_arg = list(size = 4),line_arg = list(linetype = 0), mapping=c("shape", "color")) +
  theme_bw() + theme(legend.position="bottom") +
  labs(y = "Evenness (Pielou)", x = "Treatment") +
  theme(legend.position="none")

b <- afex_plot(M3, x = "Treatment", trace = "Source",id = "Animal", dodge = 0.4, point_arg = list(size = 4),line_arg = list(linetype = 0), mapping="color") +
  theme_bw() + theme(legend.position="bottom") +
  labs(y = "Chao 1", x = "Treatment") +
  theme(legend.position="none")

d <- afex_plot(M4, x = "Treatment", trace = "Source",id = "Animal", dodge = 0.4, point_arg = list(size = 4),line_arg = list(linetype = 0), mapping="color") +
  theme_bw() + theme(legend.position="bottom") +
  labs(y = "Faith_pd", x = "Treatment") +
  theme(legend.position="none")

legend_b <- get_legend(a+ theme(legend.position = "bottom"))

e<-plot_grid(a, b, c , d, labels = "AUTO")
final<-plot_grid(e,legend_b, ncol=1,rel_heights = c(1, .1))
final


###PLOTTING OTHER

# D <- ggplot(data = metadata, aes(x = Age, y = faith_pd)) + 
#   geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
#   theme_classic() +  ylab("Faith_pd") +xlab ("Animal Age (month)") +
#   theme(axis.title.x = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14)) + 
#   theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10)) 
# theme(text=element_text(family="Times New Roman"))

# e <- ggplot(data = metadata, aes(x = Treatment, y = observed_features)) + 
#   geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
#   theme_classic() +  ylab("Observed ASVs") +xlab ("Days relative to the end of the study") +
#   theme(axis.title.x = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14)) + 
#   theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10)) 
# theme(text=element_text(family="Times New Roman"))

# l <- ggplot(data = metadata, aes(x = Ave.temp, y = faith_pd)) + 
#   geom_point() + geom_smooth(method = "lm", se=TRUE) + scale_color_brewer(palette = "Dark2") + 
#   theme_classic() +  ylab("Faith") +xlab ("Average daily temperature (°F)") +
#   theme(axis.title.x = element_text(color="black", size=14), axis.title.y = element_text(color="black", size=14)) + 
#   theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10)) 
# theme(text=element_text(family="Times New Roman"))

# ggarrange(i,j,k,l, labels = c("a", "b", "c", "d"),
#           ncol = 2, nrow=2, font.label = list(size = 20))