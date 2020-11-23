# Data analysis of the project
# NANOSED - Adsorption of metallic nanoparticles to estuarine sediments: what implication for denitrification?


# ----------------------------------------
# Denitrification rates response to Cu NPs
# ----------------------------------------


# Libraries
library(tidyverse)     # data and plots    
library(vegan)         # ecological multivariate analyses
library(dabestr)       # test significance no p-value based
library(envalysis)     # environmental analysis
library(scales)        # map data to aesthetics
library(cowplot)       # arranging plots

set.seed(57)

# Import data
dnf_data <- read.csv("Copper_nanoparticle_estuary_denitrification.csv")

dnf_data$Conc <- factor(dnf_data$Conc)

p14<-ggplot(dnf_data, aes(x=Conc, y=DNF, fill= Conc)) + 
  geom_dotplot(binaxis="y", stackdir="center") +
  geom_boxplot() +
  ylab(expression(nmol~N~~g^{-1}~h^{-1})) +
  xlab(expression(Cu~NPs~(µg~g^{-1}~sediment))) +
  scale_fill_manual(values = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF")) +
  theme_publish() +
  theme(
    axis.title.x = element_text(size=16),
    axis.text.x  = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.y  = element_text(size=16),
    legend.title = element_text(size=16),
    legend.text = element_text(size=16),
    legend.position ="none")

p15<-ggplot(dnf_data, aes(x=Conc, y=NO2sed, fill= Conc)) + 
  geom_dotplot(binaxis="y", stackdir="center") +
  geom_boxplot() +
  ylab(expression(nmol~N-NO[2]^{"-"}~~g^{-1}~h^{-1})) +
  xlab(expression(Cu~NPs~(µg~g^{-1}~sediment))) +
  scale_fill_manual(values = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF")) +
  theme_publish() +
  theme(
    axis.title.x = element_text(size=16),
    axis.text.x  = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.y  = element_text(size=16),
    legend.title = element_text(size=16),
    legend.text = element_text(size=16),
    legend.position ="none")

plot_grid6 <- plot_grid(p14,p15,
                        ncol=2, nrow=1)
plot_grid6

# Are these differences significant?

kruskal.test(DNF ~ Conc, data = dnf_data) 
kruskal.test(NO2sed ~ Conc, data = dnf_data) 
