#path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Public-Repo"
path_to_repo = "C:/Users/cmptrsn2/OneDrive - University of Illinois - Urbana/Ch1/Public-Repo"
setwd(path_to_repo)

library(ggplot2)
library(ggfortify)
library(patchwork)
library(ggrepel)
library(tidyr)
library(dplyr)

# treatment names
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

## read in all soil data
soil.dat = read.csv("Soil_Analysis/Clean_Data/All_Soil_Data_2023.csv")

## read in biomass data
bm.dat = read.csv("Understory_Analysis/Clean_Data/Biomass_By_Quadrat_and_Year.csv")
bm.dat.2022 = bm.dat[which(bm.dat$year == 2022),-which(colnames(bm.dat) == "year")]

# soil data labels
ag.cols = c("fPOM","g475mm","g2mm","g250um","g53um","l53um","MWD")
meq.cols = c("k.meq","ca.meq","mg.meq","h.meq")
text.cols = c("Sand","Silt","Clay")

# combine biomass and soil dta
all.dat = left_join(soil.dat[,-which(colnames(soil.dat) %in% c(meq.cols,"volumetric.moisture","ton","pon"))], bm.dat.2022, 
                    by=c("trt","trt.full","num","quad"))

# remove reference
df = all.dat[-which(all.dat$trt == "R"),]

# update column names
colnames(df)[6:40] = c("Temp","Moisture","BD","TN","TC","CN","TOC",
                              text.cols,"pH","P","K","Ca","Mg","SOM",
                              "NO3","NH4","CEC","POC",ag.cols,"TIC",
                              "MAOC","HB","HL","FWD",
                              "MB","PAB","HJB")

# run PCA
PCA = prcomp(~ TIC+POC+MAOC+CN+MWD+CEC+Sand+Clay+pH+P+K+Ca+Mg+TN+NO3+NH4+BD+Moisture+Temp+HL+FWD+MB+PAB+HJB,
                 data=df, scale=T, center=T)

# plot PCA with autoplot
df$trt.full = factor(df$trt.full, levels=trt.names)
df$num = factor(df$num)
autoplot(PCA, x=1, y=2, data=df, 
         loadings=T, loadings.label=T, size=4, loadings.label.colour='black',
         col="trt.full", shape="num", loadings.label.size=4, 
         loadings.label.vjust = -0.15, loadings.label.hjust = 0.1) +
         labs(color="Treatment", shape="Plot")

autoplot(PCA, x=1, y=3, data=df, 
         loadings=T, loadings.label=T, size=4, loadings.label.colour='black',
         col="trt.full", shape="num", loadings.label.size=4, 
         loadings.label.vjust = -0.15, loadings.label.hjust = 0.1) +
  labs(color="Treatment", shape="Plot")

autoplot(PCA, x=2, y=3, data=df, 
         loadings=T, loadings.label=T, size=4, loadings.label.colour='black',
         col="trt.full", shape="num", loadings.label.size=4, 
         loadings.label.vjust = -0.15, loadings.label.hjust = 0.1) +
         labs(color="Treatment", shape="Plot")

# look at rotation vectors
summary(PCA)
PCA$rotation[,1:5]

sort(PCA$rotation[,1]) # Silt/BD - SOM/CEC
sort(PCA$rotation[,2]) # Sand - Silt - Clay/P
sort(PCA$rotation[,3]) # Clay/sand - CEC/Silt

# Percentage of variance explained
eigenvalues <- PCA$sdev^2
plot(eigenvalues/sum(eigenvalues), type = "b",
     xlab = "Principal Component",
     ylab = "Percentage of Variance Explained", ylim=c(0,0.2))
round(eigenvalues/sum(eigenvalues)*100,2)
round(cumsum(eigenvalues/sum(eigenvalues)*100),2)

## plot PCA with ggplot/ggrepel

#PCA loadings
PCAloadings <- data.frame(Variables = rownames(PCA$rotation), PCA$rotation)

# add PCA scores to the dataset
df[, c('PC1','PC2','PC3')] = PCA$x[, 1:3]

# save variable loadings in a separate dataset
rot = as.data.frame(PCA$rotation[, 1:3])
rot$var = rownames(PCA$rotation)

# rescale the loadings to fit nicely within the scatterplot of our data
mult = max(abs(df[, c('PC1','PC2','PC3')])) / max(abs(rot[, 1:3])) / 2
rot[, 1:3] = rot[, 1:3] * mult

# make plot
p1 = ggplot(data=rot, 
       aes(x=0, y=0, xend=PC1, yend=PC2, label=var)) +
       geom_point(data=df, aes(x=PC1, y=PC2, color=trt.full, shape=num), 
                  inherit.aes=FALSE, size=3) +
       geom_segment(color='red', arrow=arrow(length=unit(0.03,"npc"))) +
       geom_label_repel(aes(PC1 * 1.001, PC2 * 1.001)) +
       theme_bw() + labs(color="Treatment", shape="Plot") +
       scale_y_continuous(limits=c(-4,4)) + 
       scale_x_continuous(limits=c(-4,4)) + 
       labs(x="PC1 (16.8%)",y="PC2 (12.6%)") +
       theme(text = element_text(size=14),
             panel.grid = element_blank(),legend.position='none')
p2 = ggplot(data=rot, 
            aes(x=0, y=0, xend=PC1, yend=PC3, label=var)) +
      geom_point(data=df, aes(x=PC1, y=PC3, color=trt.full, shape=num), 
                 inherit.aes=FALSE, size=3) +
      geom_segment(color='red', arrow=arrow(length=unit(0.03,"npc"))) +
      geom_label_repel(aes(PC1 * 1.001, PC3 * 1.001)) +
      theme_bw() + labs(color="Treatment", shape="Plot") +
      scale_y_continuous(limits=c(-4,4)) + 
      scale_x_continuous(limits=c(-4,4)) + 
      labs(x="PC1 (16.8%)",y="PC3 (10.7%)") +
      theme(text = element_text(size=14),
            panel.grid = element_blank())
p3 = p1 + p2

# write plot
ggsave("Figures/Figure5_Soil_Vegetation_PCA.jpeg", 
       plot = p3, width = 34, height = 14, units="cm")



