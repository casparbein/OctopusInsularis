## load packages
library(stringr)
library(ggplot2)

## FIT-Values from vcftools
## 69inds_40MD
fis <- read.table("~/LMU/octopus/Fis/69inds_40MD.het", header = T)
View(fis)
fis[,1]

## 64inds_20MD
fis_64 <- read.table("~/LMU/octopus/Fis/64inds_20MD.het", header = T)
View(fis_64)
fis_64[,1]


## Missing data per individual, from vcftools, 69inds_40MD
miss <- read.table("~/LMU/octopus/Fis/69inds_40MD_miss.imiss", header = T)
View(miss)

miss$original_reads <- read_number_rem$reads_raw
View(miss)

ggplot(miss, aes(F_MISS, original_reads)) +  
  geom_point(aes(color = population, size = 5)) + 
  colScale

miss$population <- read_number_rem$population
plot(miss$N_MISS, miss$original_reads, color=miss$INDV)

## Missing data per individual, from vcftools, 64inds_20MD
miss_64 <- read.table("~/LMU/octopus/Fis/64inds_20MD_miss.imiss", header = T)
View(miss)

## assign locations
## this is a convenient script for cluster assignment including the function "pops_creator"
source("~/LMU/octopus/R scripts/stats_ggplot_clean.R")
fis <- pops_creator(aframe = fis, cat = fis[,1])
fis_64 <- pops_creator(aframe = fis_64, cat = fis_64[,1])

## assign clusters, 69inds_40MD, names as in the manuscript
fis$population[fis$population == "AR" | fis$population == "CE" | fis$population == "FN" |
                 fis$population == "RN"] <- "N-Coastal"

fis$population[fis$population == "AL" | fis$population == "BA" ] <- "S-Coastal"

fis$population[fis$population == "TM"] <- "S-Oceanic"

fis$population[fis$population == "ASC" | fis$population == "STH" ] <- "N-Atlantic"

fis$population[fis$population == "SPS" ] <- "N-SPS"

fis$population[fis$population == "OIC" ] <- "N-Carribean"

## assign clusters, 64inds_20MD
fis_64$population[fis_64$population == "AR" | fis_64$population == "CE" | fis_64$population == "FN" |
                 fis_64$population == "RN"] <- "N-Coastal"

fis_64$population[fis_64$population == "TM"] <- "S-Oceanic"

fis_64$population[fis_64$population == "AL" | fis_64$population == "BA" ] <- "S-Coastal"

fis_64$population[fis_64$population == "ASC" ] <- "N-Atlantic"

fis_64$population[fis_64$population == "SPS" ] <- "N-SPS"

fis_64$population[fis_64$population == "OIC" ] <- "N-Carribean"


## add a missing data variable
fis$miss <- miss$F_MISS
fis_64$miss <- miss_64$F_MISS

##piSNP values from DnaSP
pis <- read.csv("~/LMU/octopus/Fis/pairwise_n_69_40_all.csv", sep = ";", header = T)
View(pis)

## assign locations

source("C:/Users/aquan/Documents/LMU/octopus/R scripts/stats_ggplot_clean.R")
pis<- pops_creator(aframe = pis, cat = pis[,1])

## assign clusters

pis$population2 <- c(rep("S-Coastal", 4), rep("N-Coastal", 10), rep("N-Atlantic", 2), rep("S-Coastal",9), 
                     rep("N-Coastal", 10), rep("N-Coastal", 6), rep("N-Carribean",3), rep("N-Coastal",9), rep("N-SPS",3), 
                     rep("N-Atlantic",1), rep("S-Oceanic", 12))

pis$miss <- miss$F_MISS 

## AMOVA of Fit values
View(fis)
anova_fis <- aov(F ~ population, data = fis)
summary(anova_fis)

## pairwise t-test
## Fit
pwt <- pairwise.t.test(fis$F, fis$population, alternative="two.sided", p.adjust.method="BH", pool.sd = F)
pwt$p.value

## piSNP
pwt_pi <- pairwise.t.test(pis$X69inds, pis$population2, alternative="two.sided", p.adjust.method="BH", pool.sd = F)
pwt_pi$p.value

## plot missing data against Fit
missing_plot <- ggplot(fis, aes(population, miss)) + 
  geom_boxplot(fill= c("#A945FF", "#ae5d5d" , "#1D72F5","#77CE61",
                       "#FF9326",
                       "#DF0101"),  outlier.shape = NA) + geom_jitter()

missing_plot

## plot Fit stratified by cluster assignment
fit_plot <- ggplot(fis, aes(population, F)) + 
  geom_boxplot(fill= c("#A945FF", "#ae5d5d" , "#1D72F5","#77CE61",
                       "#FF9326",
                        "#DF0101"),  outlier.shape = NA) +
  geom_jitter(color="black", size=1, alpha=0.9) + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank())+
  ylab(expression(F[I][T])) + 
  xlab("cluster")+
  #ggtitle(label = expression(paste("Boxplots of ", F[I][T],"values of individuals",sep=" "), 
  #                   subtitle = "grouped by cluster ID"))
  ggtitle(label= title, subtitle = "grouped by cluster ID")
  
title = expression(paste(F[I][T]," - values of individuals",sep=" "))
fit_plot


fit_plot_64 <- ggplot(fis_64, aes(population, F)) + 
  geom_boxplot(fill= c("#A945FF", "#ae5d5d" , "#1D72F5","#77CE61",
                       "#FF9326",
                       "#DF0101"),  outlier.shape = NA) +
  geom_jitter(color="black", size=1, alpha=0.9) + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x=element_blank())+
  ylab(expression(F[I][T])) + 
  xlab("cluster")+
  #ggtitle(label = expression(paste("Boxplots of ", F[I][T],"values of individuals",sep=" "), 
  #                   subtitle = "grouped by cluster ID"))
  ggtitle(label= title, subtitle = "grouped by cluster ID")

title = expression(paste(F[I][T]," - values of individuals",sep=" "))
fit_plot_64


## Averages, 69inds_40MD
## calculate an average Fit value for each cluster
tapply(fis$F, fis$population, mean)

## calculate average missingness for each cluster
tapply(fis$miss, fis$population, mean)

## average number of sites
tapply(fis$N_SITES, fis$population, mean)

## average exp homozygotic sites
tapply(fis$E.HOM., fis$population, mean)

## average obs homozygotic sites
tapply(fis$O.HOM., fis$population, mean)

## Averages, 64inds_20MD
## calculate an average Fit value for each cluster
tapply(fis_64$F, fis_64$population, mean)

## calculate average missingness for each cluster
tapply(fis_64$miss, fis_64$population, mean)

## average number of sites
tapply(fis_64$N_SITES, fis_64$population, mean)

## average exp homozygotic sites
tapply(fis_64$E.HOM., fis_64$population, mean)

## average obs homozygotic sites
tapply(fis_64$O.HOM., fis_64$population, mean)

## average obs heterozygosity
tapply(fispis_64$obs_het, fispis_64$cluster, mean)
 
## average exp heterozygosity
tapply(fispis_64$exp_het, fispis_64$cluster, mean)

## calculate obs/exp heterozygosity, 69inds_40MD

obs_het <- 1 - ((fis$O.HOM.)/(fis$N_SITES)) 
obs_het

exp_het <- 1 - ((fis$E.HOM.)/(fis$N_SITES)) 
exp_het

1-(obs_het/exp_het)

## 64inds_20MD
obs_het_64 <- 1 - ((fis_64$O.HOM.)/(fis_64$N_SITES)) 
obs_het_64

exp_het_64 <- 1 - ((fis_64$E.HOM.)/(fis_64$N_SITES)) 
exp_het_64

1-(obs_het_64/exp_het_64)


## piSNP plots
pair_pi_plot_69 <- ggplot(pis, aes(population2, X69inds)) + 
  geom_boxplot(fill = c("#A945FF", "#ae5d5d" , "#1D72F5","#77CE61",
                        "#FF9326",
                        "#DF0101"), outlier.shape = NA) +
  geom_jitter(color="black", size=1, alpha=0.9) + 
  theme_bw()+  
  ylab("pairwise nucleotide difference (SNP)") + 
  xlab("cluster")+
  #ggtitle(label = expression(paste("Boxplots of ", F[I][T],"values of individuals",sep=" "), 
  #                   subtitle = "grouped by cluster ID"))
  ggtitle(label= "Nucleotide diversity of samples \nbased on SNPs", subtitle = "grouped by cluster ID")
pair_pi_plot_69

pair_pi_plot_64 <- ggplot(pis, aes(population2, X64inds)) + 
  geom_boxplot(fill = c("#A945FF", "#ae5d5d" , "#1D72F5","#77CE61",
                        "#FF9326",
                        "#DF0101"), outlier.shape = NA) +
  geom_jitter(color="black", size=1, alpha=0.9) + 
  theme_bw()+  
  ylab("pairwise nucleotide difference (SNP)") + 
  xlab("cluster")+
  #ggtitle(label = expression(paste("Boxplots of ", F[I][T],"values of individuals",sep=" "), 
  #                   subtitle = "grouped by cluster ID"))
  ggtitle(label= "Nucleotide diversity of samples \nbased on SNPs", subtitle = "grouped by cluster ID")
pair_pi_plot_64

## combine the plots
library(cowplot)
ggdraw() + 
  draw_plot(fit_plot  + theme(axis.title.x = element_blank()#, axis.title.y = element_blank()
), x = 0, y = 0.5, 
height = 0.45, width = 0.5) + 
  draw_plot(fit_plot_64  + theme(axis.title.x = element_blank()#, axis.title.y = element_blank()
  ), x = 0.5, y = 0.5, 
  height = 0.45, width = 0.5) + 
  draw_plot(pair_pi_plot_69  + theme(axis.title.x = element_blank(), #axis.title.y = element_blank()
  ), x = 0, y = 0.0, 
  height = 0.5, width = 0.5) +
  draw_plot(pair_pi_plot_64  + theme(axis.title.x = element_blank(), #axis.title.y = element_blank()
  ), x = 0.5, y = 0.0, 
  height = 0.5, width = 0.5) +
  draw_plot_label(
    c("A", "B", "Dataset: 69inds_40MD", "Dataset: 64inds_20MD"),
    c(0.0, 0.5, 0.1, 0.6),
    c(0.99, 0.99, 0.99, 0.99),
    size = c(11, 11, 14, 14)
  )

## regression analysis: 
## Test if there is a significant relationsship between individual missingness
## and Fit/piSNP

pis

library(lme4)
normal_model_Fit <- lm(F ~ miss, data = fis)
summary(normal_model_Fit)

normal_model_Pis <- lm(X69inds ~ miss, data = pis)
summary(normal_model_Pis)

## plot regression
fis_cols <- assign_cols(fis, 11, guide = "True", fill = "False")

ggplot(fis, aes(miss, F)) + 
  geom_point(aes(color = population)) + 
  geom_smooth(method='lm', se = F) + 
  fis_cols

ggplot(pis, aes(miss, X69inds)) + 
  geom_point(aes(color = population), size = 6) + 
  geom_smooth(method='lm', se = F) + 
  fis_cols

## create tables for the manuscript
pis_64 <- pis[complete.cases(pis), -2] 
pis_69 <- pis[,-3]

pis_64$miss <- miss_64$F_MISS 
View(pis_64)
View(fis_64)

library(tidyverse)
fispis <- fis %>% 
  inner_join(pis_69, by = c("population"="population2", "miss"="miss", "INDV"="X")) %>%
  rename("locality"="population.y",
         "pi_SNP" = "X69inds",
         "fit"="F",
         "ind"="INDV",
         "cluster"="population",
         "num_sites"="N_SITES",
         "exp_hom"="E.HOM.",
         "obs_hom"="O.HOM."
         ) %>%
  relocate(locality, .after=ind) %>%
  relocate(cluster, .after=locality) %>%
  relocate(num_sites, .after=cluster) %>%
  relocate(pi_SNP, .after=fit) 

fispis$obs_het <- obs_het
fispis$exp_het <- exp_het

View(fispis)


fispis_64 <- fis_64 %>% 
  inner_join(pis_64, by = c("population"="population2", "miss"="miss", "INDV"="X")) %>%
  rename("locality"="population.y",
         "pi_SNP" = "X64inds",
         "fit"="F",
         "ind"="INDV",
         "cluster"="population",
         "num_sites"="N_SITES",
         "exp_hom"="E.HOM.",
         "obs_hom"="O.HOM."
  ) %>%
  relocate(locality, .after=ind) %>%
  relocate(cluster, .after=locality) %>%
  relocate(num_sites, .after=cluster) %>%
  relocate(pi_SNP, .after=fit) 

fispis_64$obs_het <- obs_het_64
fispis_64$exp_het <- exp_het_64

View(fispis_64)


fispis <- fispis%>% 
  mutate(exp_hom=NULL,
        obs_hom=NULL) %>%
  relocate(obs_het, .after=num_sites) %>%
  relocate(exp_het, .after=obs_het)
View(fispis)

fispis_64 <- fispis_64%>% 
  mutate(exp_hom=NULL,
         obs_hom=NULL) %>%
  relocate(obs_het, .after=num_sites) %>%
  relocate(exp_het, .after=obs_het)
View(fispis_64)

## export
setwd("~/LMU/octopus/tables")
write.table(fispis_64, "het_etc_64.csv", sep = ";", row.names = F,
            dec = ".")

## export p-values
write.table(pwt_pi$p.value, "pisnp_pval.csv", sep = ";", row.names = T,
            dec = ".")
write.table(pwt$p.value, "fit_pval.csv", sep = ";", row.names = T,
            dec = ".")


