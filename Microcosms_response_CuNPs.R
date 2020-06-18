# Data analysis of the project
# NANOSED - Adsorption of metallic nanoparticles to estuarine sediments: what implication for denitrification?


# -----------------------------
# Microcosms response to Cu NPs
# -----------------------------


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


# --------------------------------------
# Changes in microcosms due to salinity
# --------------------------------------


# Select data for salinity comparison
data_sal <- data %>% 
  select(EXP, Time, Conc, Size, Salinity, O2, pH, OM, 
         NO3wat, NO2wat, NH4wat, NO3sed, NO2sed, NH4sed,
         Cuenv, Feenv) %>% 
  filter(EXP == "EXP1" | EXP == "EXP2")

# Transform data with z-score
data_zscores<- scale(data_sal[,c(6:16)],center = T, scale = T) %>% as.data.frame()

# Create factors in metadata
Metadata <- data_sal %>% 
  select(Time, Conc, Salinity)
Metadata$Time <- factor(Metadata$Time)
Metadata$Conc <- factor(Metadata$Conc)

# Calculate distance - Euclidean 
dis.eucli <- vegdist(data_zscores, method="euclidean", na.rm=TRUE)

# PCoA
pcoa.ord <- cmdscale(dis.eucli, eig = T)

# Create a dataframe to plot the PCoA 
PCoA.ord <-data.frame(PCoA1 = pcoa.ord$points[,1], PCoA2 = pcoa.ord$points[,2])
PCoA.ord.df<- merge(PCoA.ord , Metadata, by="row.names", sort = F)
PCoA.ord.df <- PCoA.ord.df[,-1]

# How much % of the variance is explained by the axes?
round(pcoa.ord$eig*100/sum(pcoa.ord$eig),1)

# Plot PCoA
ggplot(data = PCoA.ord.df, aes(PCoA1, PCoA2)) + 
  geom_point(aes(color = Time, shape = Conc, size = 1)) +
  labs( y="PCoA2 (19.6%)", x="PCoA1 (22.5%)") +
  facet_wrap(Salinity~.) +
  scale_color_manual(name="Time (h)",
                     labels = c("1","4","7","24","48","144"),
                     values =c ("#E64B35FF","#4DBBD5FF","#3C5488FF","#00A087FF","#F39B7FFF","#91D1C2FF")) +
  scale_shape_manual(name="Conc (µg/g)",
                     labels = c("0","0.01","0.1","1"),
                     values = c(16,17,15,3)) +
  guides(size = F, 
         shape = guide_legend(order=1, override.aes = list(size=3)), 
         color = guide_legend(order=2, override.aes = list(size=3))) +
  theme_publish() +
  theme(
    axis.title.x = element_text(size=16),
    axis.text.x  = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.y  = element_text(size=16),
    legend.title = element_text(size=16),
    legend.text = element_text(size=16),
    legend.position = "right",
    legend.box = "vertical") +
  ggtitle("Microcosms with Cu NPs < 50 nm - effect of salinity")

# Check for differences between Groups - PERMANOVA
adonis2(dis.eucli ~ Salinity*Time, data = Metadata, by ="margin")
adonis2(dis.eucli ~ Salinity*Conc, data = Metadata, by ="margin")

# Are differences in PERMANOVA due to the dispersion of the results?
beta_dispersion_data <- betadisper(dis.eucli, Metadata$Time)
permutest(beta_dispersion_data)

# Changes in microcosms due to salinity over time

# Create a data frame for each EXP
data_lm_EXP1 <- data_sal %>% filter(EXP == "EXP1") 
data_lm_EXP2 <- data_sal %>% filter(EXP == "EXP2")

# EXP1
x1 <- data_lm_EXP1 %>% select(Time)
y1 <- data_lm_EXP1  %>%
  select(O2, pH,  OM, NO3wat, NO2wat, NH4wat, NO3sed, NO2sed, NH4sed, Cuenv, Feenv)

out_lm1 <- data.frame(NULL)              
for (i in 1:length(y1)) {
  m1 <- summary(lm(y1[,i] ~ x1[,1]))       
  out_lm1[i, 1] <- names(y1)[i]           
  out_lm1[i, 2] <- m1$coefficients[1,1]  
  out_lm1[i, 3] <- m1$coefficients[2,1]   
  out_lm1[i, 4] <- m1$coefficients[2,4]   
}
names(out_lm1) <- c("y.variable", "intercept", "slope", "p.value")
out_lm1

# EXP2
x2 <- data_lm_EXP2 %>% select(Time)
y2 <- data_lm_EXP2  %>%
  select(O2, pH, OM, NO3wat, NO2wat, NH4wat, NO3sed, NO2sed, NH4sed, Cuenv, Feenv)

out_lm2 <- data.frame(NULL)              
for (i in 1:length(y2)) {
  m2 <- summary(lm(y2[,i] ~ x2[,1]))       
  out_lm2[i, 1] <- names(y2)[i]           
  out_lm2[i, 2] <- m2$coefficients[1,1]  
  out_lm2[i, 3] <- m2$coefficients[2,1]   
  out_lm2[i, 4] <- m2$coefficients[2,4]   
}
names(out_lm2) <- c("y.variable", "intercept", "slope", "p.value")
out_lm2

# T-test to compare slopes of two regressions 

#####
fit1_ph <- lm(pH~Time, data_lm_EXP1)
s1_ph <- summary(fit1_ph)$coefficients
fit2_ph <- lm(pH~Time, data_lm_EXP2)
s2_ph <- summary(fit2_ph)$coefficients
#
db_ph <- (s2_ph[2,1]-s1_ph[2,1])
sd_ph <- sqrt(s2_ph[2,2]^2+s1_ph[2,2]^2)
df_ph <- (fit1_ph$df.residual+fit2_ph$df.residual)
td_ph <- db_ph/sd_ph
t_test_ph <- 2*pt(-abs(td_ph), df_ph)
#####
fit1_no3 <- lm(NO3wat~Time, data_lm_EXP1)
s1_no3 <- summary(fit1_no3)$coefficients
fit2_no3 <- lm(NO3wat~Time, data_lm_EXP2)
s2_no3 <- summary(fit2_no3)$coefficients
#
db_no3 <- (s2_no3[2,1]-s1_no3[2,1])
sd_no3 <- sqrt(s2_no3[2,2]^2+s1_no3[2,2]^2)
df_no3 <- (fit1_no3$df.residual+fit2_no3$df.residual)
td_no3 <- db_no3/sd_no3
t_test_no3 <- 2*pt(-abs(td_no3), df_no3)

#####
out_fit <- data.frame(NULL)  
out_fit[1, 1] <- t_test_ph           
out_fit[1, 2] <- t_test_no3 
names(out_fit) <- c("pH", "NO3 water")
out_fit


# -----------------------------------------------
# Changes in microcosms with different NPs sizes 
# -----------------------------------------------


data_np <- data %>% 
  select(EXP, Time, Conc, Size, Salinity, O2, pH, OM, 
         NO3wat, NO2wat, NH4wat, NO3sed, NO2sed, NH4sed,
         Cuenv, Feenv, Cuexc, Feexc, DAPI) %>%
  filter(EXP == "EXP50" | EXP == "EXP150")

# Normalize data with z-score
data_np_zscores <- scale(data_np[,c(6:18)], center = T, scale = T) %>% as.data.frame()

# Create factors in metadata
Metadata_np <- data_np %>% 
  select(Time, Conc, Size)
Metadata_np$Time <- factor(Metadata_np$Time)
Metadata_np$Conc <- factor(Metadata_np$Conc)

# Calculate distances - Euclidean 
dis.eucli.np <- vegdist(data_np_zscores, method="euclidean", na.rm=TRUE)

# PCoA 
pcoa.ord.np <- cmdscale(dis.eucli.np, eig = T)

# Create a dataframe to plot the PCoA 

# Create a dataframe to plot the PCoA 
PCoA.ord.np <-data.frame(PCoA1 = pcoa.ord.np$points[,1], PCoA2 = pcoa.ord.np$points[,2])
PCoA.ord.np.df<- merge(PCoA.ord.np, Metadata_np, by="row.names", sort = F)
PCoA.ord.np.df <- PCoA.ord.np.df[,-1]

# How much % of the variance is explained by the axes?
round(pcoa.ord.np$eig*100/sum(pcoa.ord.np$eig),1)

# Order the NPs sizes
PCoA.ord.np.df$Size = factor(PCoA.ord.np.df$Size, levels=c("50nm","150nm"))

# Plot PCoA
ggplot(data = PCoA.ord.np.df, aes(PCoA1, PCoA2)) + 
  geom_point(aes(color = Size, shape = Time, size = 1)) +
  stat_ellipse(aes(colour = Size, group = Size)) +
  labs( y="PCoA2 (16.1%)", x="PCoA1 (28.4%)") +
  scale_color_manual(name="Size",
                     labels = c("50nm","150nm"),
                     values = c("#FBB4AE","#B3CDE3")) +
  scale_shape_manual(name="Time (h)",
                     labels = c("1","4","7","24","48","144"),
                     values =c (15,16,17,3,4,5)) +
  guides(size = F, 
         shape = guide_legend(order=2, override.aes = list(size=3)), 
         color = guide_legend(order=1, override.aes = list(size=2))) +
  theme_publish() +
  theme(
    axis.title.x = element_text(size=16),
    axis.text.x  = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.y  = element_text(size=16),
    legend.title = element_text(size=16),
    legend.text = element_text(size=16),
    legend.position = "right",
    legend.box = "vertical") +
  ggtitle("Microcosms with Cu NPs < 50 nm and < 150 nm")

# Check for differences between Groups - PERMANOVA
adonis2(dis.eucli.np ~ Size * Conc, data = Metadata_np)
adonis2(dis.eucli.np ~ Size * Time, data = Metadata_np)

# Are differences in PERMANOVA due to the dispersion of the results?
beta_dispersion_np_size <- betadisper(dis.eucli.np, Metadata_np$Size)
permutest(beta_dispersion_np_size)

beta_dispersion_np_time <- betadisper(dis.eucli.np, Metadata_np$Time)
permutest(beta_dispersion_np_time)

# Changes in microcosms with different sizes over time

# Create a data frame for each EXP
data_lm_EXP50 <- data_np %>% filter(EXP == "EXP50") 
data_lm_EXP150 <- data_np %>% filter(EXP == "EXP150")

# EXP50
x1 <- data_lm_EXP50 %>% select(Time)
y1 <- data_lm_EXP50  %>%
  select(O2, pH, OM, NO3wat, NO2wat, NH4wat, NO3sed, NO2sed, NH4sed, Cuenv, Feenv, Cuexc, Feexc, DAPI)

out_lm50 <- data.frame(NULL)              
for (i in 1:length(y1)) {
  m1 <- summary(lm(y1[,i] ~ x1[,1]))       
  out_lm50[i, 1] <- names(y1)[i]           
  out_lm50[i, 2] <- m1$coefficients[1,1]  
  out_lm50[i, 3] <- m1$coefficients[2,1]   
  out_lm50[i, 4] <- m1$coefficients[2,4]   
}
names(out_lm50) <- c("y.variable", "intercept", "slope", "p.value")
out_lm50

# EXP150
x2 <- data_lm_EXP150 %>% select(Time)
y2 <- data_lm_EXP150 %>%
  select(O2, pH, OM, NO3wat, NO2wat, NH4wat, NO3sed, NO2sed, NH4sed, Cuenv, Feenv, Cuexc, Feexc, DAPI)

out_lm150 <- data.frame(NULL)              
for (i in 1:length(y2)) {
  m2 <- summary(lm(y2[,i] ~ x2[,1]))       
  out_lm150[i, 1] <- names(y2)[i]           
  out_lm150[i, 2] <- m2$coefficients[1,1]  
  out_lm150[i, 3] <- m2$coefficients[2,1]   
  out_lm150[i, 4] <- m2$coefficients[2,4]   
}
names(out_lm150) <- c("y.variable", "intercept", "slope", "p.value")
out_lm150

# T-test to compare slopes of two regressions 

#####
fit1_ph <- lm(pH~Time, data_lm_EXP50)
s1_ph <- summary(fit1_ph)$coefficients
fit2_ph <- lm(pH~Time, data_lm_EXP150)
s2_ph <- summary(fit2_ph)$coefficients
#
db_ph <- (s2_ph[2,1]-s1_ph[2,1])
sd_ph <- sqrt(s2_ph[2,2]^2+s1_ph[2,2]^2)
df_ph <- (fit1_ph$df.residual+fit2_ph$df.residual)
td_ph <- db_ph/sd_ph
t_test_ph <- 2*pt(-abs(td_ph), df_ph)
#####
fit1_no3 <- lm(NO3wat~Time, data_lm_EXP50)
s1_no3 <- summary(fit1_no3)$coefficients
fit2_no3 <- lm(NO3wat~Time, data_lm_EXP150)
s2_no3 <- summary(fit2_no3)$coefficients
#
db_no3 <- (s2_no3[2,1]-s1_no3[2,1])
sd_no3 <- sqrt(s2_no3[2,2]^2+s1_no3[2,2]^2)
df_no3 <- (fit1_no3$df.residual+fit2_no3$df.residual)
td_no3 <- db_no3/sd_no3
t_test_no3 <- 2*pt(-abs(td_no3), df_no3)
#####
fit1_no2 <- lm(NO2sed~Time, data_lm_EXP50)
s1_no2 <- summary(fit1_no2)$coefficients
fit2_no2 <- lm(NO2sed~Time, data_lm_EXP150)
s2_no2 <- summary(fit2_no2)$coefficients
#
db_no2 <- (s2_no2[2,1]-s1_no2[2,1])
sd_no2 <- sqrt(s2_no2[2,2]^2+s1_no2[2,2]^2)
df_no2 <- (fit1_no2$df.residual+fit2_no2$df.residual)
td_no2 <- db_no2/sd_no2
t_test_no2 <- 2*pt(-abs(td_no2), df_no2)
#####
fit1_fe <- lm(Feexc~Time, data_lm_EXP50)
s1_fe <- summary(fit1_fe)$coefficients
fit2_fe <- lm(Feexc~Time, data_lm_EXP150)
s2_fe <- summary(fit2_fe)$coefficients
#
db_fe <- (s2_fe[2,1]-s1_fe[2,1])
sd_fe <- sqrt(s2_fe[2,2]^2+s1_fe[2,2]^2)
df_fe <- (fit1_fe$df.residual+fit2_fe$df.residual)
td_fe <- db_fe/sd_fe
t_test_fe <- 2*pt(-abs(td_fe), df_fe)

#####
out_fit <- data.frame(NULL)  
out_fit[1, 1] <- t_test_ph           
out_fit[1, 2] <- t_test_no3 
out_fit[1, 3] <- t_test_no2
out_fit[1, 4] <- t_test_fe
names(out_fit) <- c("pH", "NO3 water", "NO2 sed","Fe exchangeable")
out_fit


# -------------------------------
# Changes in the number of cells
# -------------------------------


# Create factors in data_np
data_np$Size <- factor(data_np$Size, levels=c("50nm","150nm"))
data_np$Conc <- factor(data_np$Conc)

# Evaluate increases in cell number

# Test for data normality with Anova & Shapiro-Wilk test (if not normal use Kruskal-wallis)
dapi.aov <- aov(DAPI ~ Conc, data = data_np)
dapi.aov_residuals <- residuals(object = dapi.aov)
shapiro.test(x = dapi.aov_residuals)

# Kruskal-wallis
kruskal.test(DAPI ~ Conc, data = data_np)

# If significant compare control without NPs with microcosms with NPs - Wilcoxon signed-rank test
pairwise.wilcox.test(data_np$DAPI, data_np$Conc, p.adjust.method = "BH")

# Estimation plot

DAPI_best <- dabest(data_np, Conc, DAPI,
                    idx = c("0", "0.01", "0.1", "1"),
                    paired = FALSE)

plot(DAPI_best, color.column = Size,
     rawplot.markersize = 3,
     rawplot.ylabel = "DAPI-stained cells \n (number/g sediment)",
     effsize.ylabel = "Unpaired \n mean difference",
     axes.title.fontsize = 10,
     palette = "Pastel1")


# ------------------------
# Changes in metal content
# ------------------------


# Evaluate increases in metal concentration

# Test for data normality with Anova & Shapiro-Wilk test
cuenv.aov <- aov(Cuenv ~ Conc, data = data_np)
cuenv.aov_residuals <- residuals(object = cuenv.aov)
shapiro.test(x = cuenv.aov_residuals)

# Kruskal-wallis test
kruskal.test(Cuenv ~ Conc, data = data_np)

# Test for data normality with Anova & Shapiro-Wilk test
cuexc.aov <- aov(Cuexc ~ Conc, data = data_np)
cuexc.aov_residuals <- residuals(object = cuexc.aov)
shapiro.test(x = cuexc.aov_residuals)

# Kruskal-wallis test
kruskal.test(Cuexc ~ Conc, data = data_np)
pairwise.wilcox.test(data_np$Cuexc, data_np$Conc, p.adjust.method = "BH")

# Test for data normality with Anova & Shapiro-Wilk test
feenv.aov <- aov(Feenv ~ Conc, data = data_np)
feenv.aov_residuals <- residuals(object = feenv.aov)
shapiro.test(x = feenv.aov_residuals)

# Kruskal-wallis test (data is not normal)
kruskal.test(Feenv ~ Conc, data = data_np)

# Test for data normality with Anova & Shapiro-Wilk test
feexc.aov <- aov(Feexc ~ Conc, data = data_np)
feexc.aov_residuals <- residuals(object = feexc.aov)
shapiro.test(x = feexc.aov_residuals)

# Kruskal-wallis test
kruskal.test(Feexc ~ Conc, data = data_np)

# Estimation plots

Cu_env <- dabest(data_np, Conc, Cuenv,
                 idx = c("0", "0.01", "0.1", "1"),
                 paired = FALSE)
Cu_exc <- dabest(data_np, Conc, Cuexc,
                 idx = c("0", "0.01", "0.1", "1"),
                 paired = FALSE)
Fe_env <- dabest(data_np, Conc, Feenv,
                 idx = c("0", "0.01", "0.1", "1"),
                 paired = FALSE)
Fe_exc <- dabest(data_np, Conc, Feexc,
                 idx = c("0", "0.01", "0.1", "1"),
                 paired = FALSE)

p1<- plot(Cu_env, color.column = Size,
          rawplot.markersize = 2,
          rawplot.ylabel = "Cu environmental \n (µg/g sediment)",
          effsize.ylabel = "Unpaired \n mean difference",
          axes.title.fontsize = 10,
          palette = "Pastel1")
p2<- plot(Cu_exc, color.column = Size,
          rawplot.markersize = 2,
          rawplot.ylabel = "Cu exchangeable \n (µg/g sediment)",
          effsize.ylabel = "Unpaired \n mean difference",
          axes.title.fontsize = 10,
          palette = "Pastel1") 
p3<- plot(Fe_env, color.column = Size,
          rawplot.markersize = 2,
          rawplot.ylabel = "Fe environmental \n (mg/g sediment)",
          effsize.ylabel = "Unpaired \n mean difference",
          axes.title.fontsize = 10,
          palette = "Pastel1")
p4<- plot(Fe_exc, color.column = Size,
          rawplot.markersize = 2,
          rawplot.ylabel = "Fe exchangeable \n (mg/g sediment)",
          effsize.ylabel = "Unpaired \n mean difference",
          axes.title.fontsize = 10,
          palette = "Pastel1")

plot_grid1 <- plot_grid(p1, p2, p3, p4, 
                        ncol=2, nrow=2)
plot_grid1

# Normalise Cu concentration

data_np <- transform(data_np, ratiocuenv = Cuenv / Feenv)

ratio_Cu <- dabest(data_np, Conc, ratiocuenv,
                   idx = c("0", "0.01", "0.1", "1"),
                   paired = FALSE)
plot(ratio_Cu , color.column = Size,
     rawplot.markersize = 3,
     rawplot.ylabel = "Cu/Fe environmental \n fraction",
     effsize.ylabel = "Unpaired \n mean difference",
     axes.title.fontsize = 10,
     palette = "Pastel1")
