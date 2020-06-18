# Data analysis of the project
# NANOSED - Adsorption of metallic nanoparticles to estuarine sediments: what implication for denitrification?


# ----------------------------------------
# Denitrification genes response to Cu NPs
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
data <- read.csv("Copper_nanoparticle_estuary_microcosms.csv")

# ---------------------------------------
# Can 16S rRNA be used as reference gene?
# ---------------------------------------


# Import data 
data_reference <- data %>% 
  select(EXP, Time, Foldch16S)  %>%
  filter(EXP == "EXP50" | EXP == "EXP150")

ggplotRegression <- function (fit) {
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(color = "indianred") +
    ggtitle("Linear regression of fold change in 16S expression over time") +
    theme_publish() +
    ylab("Log2 Fold change") +
    xlab("Time (h)") +
    stat_smooth(method = "lm", color = "black") +
    labs(subtitle = paste("R2 = ", signif(summary(fit)$r.squared, 3),
                          "Intercept =", signif(fit$coef[[1]], 3),
                          "Slope =", signif(fit$coef[[2]], 3),
                          "P =", signif(summary(fit)$coef[2,4], 3)))
}

ggplotRegression(lm(Foldch16S ~ Time, data = data_reference))


# --------------------------
# Changes in gene expression
# --------------------------


data_gene <- data %>% 
  select(EXP, Time, Conc, Size, nosZ, nirS, nirK) %>%
  filter(EXP == "EXP50" | EXP == "EXP150") %>%
  filter(Conc != "0")

data_gene$Conc <- factor(data_gene$Conc)
data_gene$Time <- factor(data_gene$Time, levels = rev(unique(data_gene$Time)))

conc.labs <- c("0 µg/g","0.01 µg/g", "0.1 µg/g", "1 µg/g")
names(conc.labs) <- c("0","0.01", "0.1", "1")

plot_nirS50 <- ggplot(data_gene %>% filter(EXP=="EXP50"), aes(x = Time, y = nirS, fill = nirS>0)) + 
  geom_bar(stat = "identity", position = "dodge") +  
  facet_grid(~Conc, labeller = labeller(Conc = conc.labs)) +
  scale_fill_manual(values = c("#FDDBC7","#E64B35FF")) +
  ylab(expression(paste("Log"[2], " fold-change in ", italic("nirS"), " transcripts "))) + 
  xlab("Time (h)") + 
  geom_hline(yintercept = 0) + 
  coord_flip() + 
  scale_y_continuous(breaks = seq(-8, 14, by = 4)) +
  guides(fill=F) +
  theme_publish()

plot_nirS150 <- ggplot(data_gene %>% filter(EXP=="EXP150"), aes(x = Time, y = nirS, fill = nirS>0)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_grid(~Conc, labeller = labeller(Conc = conc.labs)) +
  scale_fill_manual(values = c("#FDDBC7","#E64B35FF")) +
  ylab(expression(paste("Log"[2], " (fold-change) in ", italic("nirS"), " transcripts "))) + 
  xlab("Time (h)") + 
  geom_hline(yintercept = 0) + 
  coord_flip() + 
  scale_y_continuous(breaks = seq(-6, 6, by = 2)) +
  guides(fill=F) +
  theme_publish() 

plot_nosZ50 <- ggplot(data_gene %>% filter(EXP=="EXP50"), aes(x = Time, y = nosZ, fill = nosZ>0)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_grid(~Conc, labeller = labeller(Conc = conc.labs)) +
  scale_fill_manual(values = c("#D1E5F0","#4DBBD5FF")) +
  ylab(expression(paste("Log"[2], " (fold-change) in ", italic("nosZ"), " transcripts "))) + 
  xlab("Time (h)") + 
  geom_hline(yintercept = 0) + 
  coord_flip() + 
  scale_y_continuous(breaks = seq(-8, 14, by = 4)) +
  guides(fill=F) +
  theme_publish()

plot_nosZ150 <- ggplot(data_gene %>% filter(EXP=="EXP150"), aes(x = Time, y = nosZ, fill = nosZ>0)) + 
  geom_bar(stat = "identity", position = "dodge") +  
  facet_grid(~Conc, labeller = labeller(Conc = conc.labs)) +
  scale_fill_manual(values = c("#D1E5F0","#4DBBD5FF")) +
  ylab(expression(paste("Log"[2], " (fold-change) in ", italic("nosZ"), " transcripts "))) + 
  xlab("Time (h)") + 
  geom_hline(yintercept = 0) + 
  coord_flip() + 
  scale_y_continuous(breaks = seq(-6, 6, by = 2)) +
  guides(fill=F) +
  theme_publish()

plot_nirK50 <- ggplot(data_gene %>% filter(EXP=="EXP50"), aes(x = Time, y = nirK, fill = nirK>0)) + 
  geom_bar(stat = "identity", position = "dodge") +  
  facet_grid(~Conc, labeller = labeller(Conc = conc.labs)) +
  scale_fill_manual(values = c("#AEFCEE","#00A087FF")) +
  ylab(expression(paste("Log"[2], " (fold-change) in ", italic("nirK"), " transcripts "))) + 
  xlab("Time (h)") + 
  geom_hline(yintercept = 0) + 
  coord_flip() + 
  scale_y_continuous(breaks = seq(-6, 6, by = 2)) +
  guides(fill=F) +
  theme_publish()

plot_grid2 <- plot_grid(plot_nirS50, plot_nirS150,
                        plot_nosZ50, plot_nosZ150, 
                        plot_nirK50,
                        ncol=2, nrow=3)
plot_grid2

# Differentially expressed genes

sum_fun <- function(data_gene, nosZ, nirS, nirK){
  nosZ <- enquo(nosZ)
  nirS <- enquo(nirS)
  nirK <- enquo(nirK)
  
  data_gene %>%
    summarize(avgnosZ = mean(!!nosZ, na.rm = T), sdnosZ = sd(!!nosZ, na.rm = T),
              avgnirS = mean(!!nirS, na.rm = T), sdnirS = sd(!!nirS, na.rm = T),
              avgnirK = mean(!!nirK, na.rm = T), sdnirK = sd(!!nirK, na.rm = T))
}

sum_fun(data_gene, nosZ, nirS, nirK)


# ---------------------------
# Changes in gene copy number
# ---------------------------


gene_copy <- data %>% 
  select(EXP, Time, Conc, copy16, copynos, copynir) %>%
  filter(EXP == "EXP50")

gene_copy$Conc <- factor(gene_copy$Conc)
gene_copy$Time <- factor(gene_copy$Time)

library(wesanderson)
pal <- wes_palette("Zissou1", 10, type = "continuous")

nir50_heatmap <- ggplot(data = gene_copy, 
                        mapping = aes(x = Time, y = Conc, fill = copynir)) +
  geom_tile() +
  scale_fill_gradientn(colours = pal, name = "Copy number / \ng sediment", label=scientific) + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  coord_equal() +
  ylab(expression(Cu~NPs~(µg~g^{-1}))) +
  xlab("Time (h)") +
  theme_publish() +
  theme(
    axis.title.x = element_text(size=16),
    axis.text.x  = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.y  = element_text(size=16),
    legend.title = element_text(size=14),
    legend.text = element_text(size=14),
    legend.position = "right") +
  ggtitle("nirS Cu NPs < 50 nm")

nos50_heatmap <- ggplot(data = gene_copy, 
                        mapping = aes(x = Time, y = Conc, fill = copynos)) +
  geom_tile() +
  scale_fill_gradientn(colours = pal, name = "Copy number / \ng sediment", label=scientific) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  coord_equal() +
  ylab(expression(Cu~NPs~(µg~g^{-1}))) +
  xlab("Time (h)") +
  theme_publish() +
  theme(
    axis.title.x = element_text(size=16),
    axis.text.x  = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.y  = element_text(size=16),
    legend.title = element_text(size=14),
    legend.text = element_text(size=14),
    legend.position = "right") +
  ggtitle("nosZ Cu NPs < 50 nm")

ref16S50_heatmap<- ggplot(data = gene_copy, 
                          mapping = aes(x = Time, y = Conc, fill = copy16)) +
  geom_tile() +
  scale_fill_gradientn(colours = pal, name = "Copy number / \ng sediment", label=scientific) + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  coord_equal() +
  ylab(expression(Cu~NPs~(µg~g^{-1}))) +
  xlab("Time (h)") +
  theme_publish() +
  theme(
    axis.title.x = element_text(size=16),
    axis.text.x  = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.y  = element_text(size=16),
    legend.title = element_text(size=14),
    legend.text = element_text(size=14),
    legend.position = "right") +
  ggtitle("16S Cu NPs < 50 nm")

plot_grid3 <- plot_grid(nir50_heatmap, nos50_heatmap,ref16S50_heatmap,
                        ncol=1, nrow=3)
plot_grid3

# Create data set for Log copy genes divided by Log 16S

gene_copy_nr<- log10(gene_copy %>% select(copy16, copynos, copynir))
gene_copy_nr<- merge(gene_copy %>% select(Time, Conc), gene_copy_nr, by ="row.names", sort = F)

gene_copy_nr <- transform(gene_copy_nr, lognos = copynos / copy16)
gene_copy_nr <- transform(gene_copy_nr, lognir = copynir / copy16)

plot_gene_nos <-ggplot(gene_copy_nr, aes(x = Conc, y = lognos, fill= Conc)) +
  geom_dotplot(binaxis="y", stackdir="center") +
  geom_boxplot() +
  ylab(expression(Log10~italic(nosZ)~"/"~Log10~italic("16S"))) +
  xlab(expression(Cu~NPs~(µg~g^{-1}~sediment))) +
  scale_fill_manual(values = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF")) +
  theme_publish() +
  theme(
    axis.title.x = element_text(size=16),
    axis.text.x  = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.y  = element_text(size=16),
    legend.position ="none") +
  ggtitle("nosZ / 16S Cu NPs < 50 nm")


plot_gene_nir <-ggplot(gene_copy_nr, aes(x = Conc, y = lognir, fill= Conc)) +
  geom_dotplot(binaxis="y", stackdir="center") +
  geom_boxplot() +
  ylab(expression(Log10~italic(nirS)~"/"~Log10~italic("16S"))) +
  xlab(expression(Cu~NPs~(µg~g^{-1}~sediment))) +
  scale_fill_manual(values = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF")) +
  theme_publish() +
  theme(
    axis.title.x = element_text(size=16),
    axis.text.x  = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.y  = element_text(size=16),
    legend.position ="none") +
  ggtitle("nirS / 16S Cu NPs < 50 nm")

plot_grid4 <- plot_grid(plot_gene_nir, plot_gene_nos,
                        ncol=2, nrow=1)
plot_grid4 

# Are these differences significant?

kruskal.test(lognir ~ Conc, data = gene_copy_nr) 
pairwise.wilcox.test(gene_copy_nr$lognir, gene_copy_nr$Conc, p.adjust.method = "BH")

kruskal.test(lognos ~ Conc, data = gene_copy_nr)
pairwise.wilcox.test(gene_copy_nr$lognos, gene_copy_nr$Conc, p.adjust.method = "BH")


# ----------------------------
# Gene copy number correlation
# ----------------------------


ggplotRegression2 <- function (fit) {
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(color = "indianred") +
    stat_smooth(method = "lm", color = "indianred") +
    labs(subtitle = paste("R2 =", signif(summary(fit)$r.squared, 3),
                          "P =", signif(summary(fit)$coef[2,4], 3))) +
    ylab(expression(paste("Log"[10], italic("nirS"), " abundance "))) + 
    xlab(expression(paste("Log"[10], italic("nosZ"), " abundance "))) 
  
}

p10<- ggplotRegression2(lm(copynir ~ copynos, data = gene_copy_nr[gene_copy_nr$Conc=="0",])) +
  ggtitle("0 µg/g") + theme_publish()
p11<- ggplotRegression2(lm(copynir ~ copynos, data = gene_copy_nr[gene_copy_nr$Conc=="0.01",])) +
  ggtitle("0.1 µg/g") + theme_publish()
p12<- ggplotRegression2(lm(copynir ~ copynos, data = gene_copy_nr[gene_copy_nr$Conc=="0.1",])) +
  ggtitle("0.1 µg/g") + theme_publish()
p13<- ggplotRegression2(lm(copynir ~ copynos, data = gene_copy_nr[gene_copy_nr$Conc=="1",])) +
  ggtitle("1 µg/g") + theme_publish()

plot_grid5 <- plot_grid(p10, p11, p12, p13,
                        ncol=2, nrow=2)
plot_grid5
