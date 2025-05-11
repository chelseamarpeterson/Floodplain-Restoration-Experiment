setwd("C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/TreeData")

library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)

# treatment names
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

# read in csvs
ht.dat = read.csv("ZonalSt_2020_Veg_Height.csv", header=T)
dens.dat = read.csv("ZonalSt_2020_Veg_Density.csv", header=T)
int.dat = read.csv("ZonalSt_2020_Intensity.csv", header=T)
colnames(ht.dat)[which(colnames(ht.dat) %in% c("Treatment","Plot_Num"))] = c("trt","num")
colnames(dens.dat)[which(colnames(dens.dat) %in% c("Treatment","Plot_Num"))] = c("trt","num")
colnames(int.dat)[which(colnames(int.dat) %in% c("Treatment","Plot_Num"))] = c("trt","num")

# add treatment names column
ht.dat$trt.full = rep(0, nrow(ht.dat))
dens.dat$trt.full = rep(0, nrow(dens.dat))
int.dat$trt.full = rep(0, nrow(int.dat))
for (i in 1:6) {ht.dat$trt.full[which(ht.dat$trt == trts[i])] = trt.names[i]}
for (i in 1:6) {dens.dat$trt.full[which(dens.dat$trt == trts[i])] = trt.names[i]}
for (i in 1:6) {int.dat$trt.full[which(int.dat$trt == trts[i])] = trt.names[i]}

# melt df
stats = c("MIN","MAX","RANGE","MEAN","STD")
ht.melt = melt(ht.dat[,c("trt.full","num",stats)], id.vars=c("trt.full","num"))
dens.melt = melt(dens.dat[,c("trt.full","num",stats)], id.vars=c("trt.full","num"))
int.melt = melt(int.dat[,c("trt.full","num",stats)], id.vars=c("trt.full","num"))

# combine dfs
#ht.melt$var2 = "Height"
#dens.melt$var2 = "Density"
#int.melt$var2 = "Intensity"
#all.melt = rbind(ht.melt,dens.melt,int.melt)
#colnames(all.melt)

# make ggplot
p1 = ggplot(ht.melt, aes(y=factor(trt.full, levels=trt.names), x=value)) + 
       geom_boxplot() + 
       facet_wrap(~variable, ncol=5, scales="free_x") + 
       labs(y="",x="Veg height (m)")

p2 = ggplot(dens.melt, aes(y=factor(trt.full, levels=trt.names), x=value)) + 
       geom_boxplot() + 
       facet_wrap(.~variable, ncol=5, scales="free_x") + 
       labs(y="",x="Veg density")

p3 = ggplot(int.melt, aes(y=factor(trt.full, levels=trt.names), x=value)) + 
      geom_boxplot() + 
      facet_wrap(.~variable, ncol=5, scales="free_x") + 
      labs(y="",x="Intensity")
p1/p2/p3

## plot biomass variables v. lidar data

# read in biomass data
bm.dat = read.csv("All_Biomass_Variables.csv", header=T)
colnames(bm.dat)[4:19] = c("FWD count","CWD area","Snag area","Dead stems","Live woody C stock",
                           "Live stems","SOM","FWD C stock","FWD C:N ratio","Litter C stock", 
                           "Litter C:N ratio","P. arun C stock","H. jap C stock","Mixed biomass C stock",
                           "Total herbaceous C stock","Herb. biomass C:N ratio")

# join lidar and bm data
ht.bm.dat = right_join(bm.dat, ht.dat, by=c("trt","trt.full","num"))
dens.bm.dat = right_join(bm.dat, dens.dat, by=c("trt","trt.full","num"))
int.bm.dat = right_join(bm.dat, int.dat, by=c("trt","trt.full","num"))
                         
# melt dataframes
vars = c("Live woody C stock","Total herbaceous C stock","Litter C stock","FWD C stock","SOM")
df.melt.ht = melt(ht.bm.dat[,c("trt.full","num", vars, stats)], id.vars=c("trt.full","num",vars))
df.melt.dens = melt(dens.bm.dat[,c("trt.full","num", vars, stats)], id.vars=c("trt.full","num",vars))
df.melt.int = melt(int.bm.dat[,c("trt.full","num", vars, stats)], id.vars=c("trt.full","num",vars))
df.melt2.ht = melt(df.melt.ht, id.vars=c("trt.full","num","variable","value"))
df.melt2.dens = melt(df.melt.dens, id.vars=c("trt.full","num","variable","value"))
df.melt2.int = melt(df.melt.int, id.vars=c("trt.full","num","variable","value"))

# plot lidar vars v. biomass vars
colnames(df.melt2.ht) = c("trt.full","num","var1","value1","var2","value2" )
colnames(df.melt2.dens) = c("trt.full","num","var1","value1","var2","value2" )
colnames(df.melt2.int) = c("trt.full","num","var1","value1","var2","value2" )
ggplot(df.melt2.ht, aes(x=value1, y=value2, color=factor(trt.full, levels=trt.names))) + 
        geom_point() + theme(legend.title = element_blank()) +
        labs(x="Veg height (m)", y="", title="C stock v. height stats") +
        facet_grid(var2~var1, scales="free")

ggplot(df.melt2.dens, aes(x=value1, y=value2, color=factor(trt.full, levels=trt.names))) + 
        geom_point() + theme(legend.title = element_blank()) +
        labs(x="Veg density", y="", title="C stock v. density stats") +
        facet_grid(var2~var1, scales="free")

ggplot(df.melt2.int, aes(x=value1, y=value2, color=factor(trt.full, levels=trt.names))) + 
        geom_point() + theme(legend.title = element_blank()) +
        labs(x="Intensity", y="", title="C stock v. intensity stats") +
        facet_grid(var2~var1, scales="free")


## make correlation plots
library(corrplot)

vars2 = c("Live woody C stock","Snag area","CWD area","Live stems","Dead stems",
          "Total herbaceous C stock","P. arun C stock","H. jap C stock",
          "Mixed biomass C stock","Herb. biomass C:N ratio",
          "Litter C stock","Litter C:N ratio","FWD C stock","FWD C:N ratio")
ht.cor = cor(ht.bm.dat[,c(vars2,stats)])
dens.cor = cor(dens.bm.dat[,c(vars2,stats)])
int.cor = cor(int.bm.dat[,c(vars2,stats)])
corrplot(ht.cor, type="upper", title="Vegetation data v. height stats")
corrplot(dens.cor, type="upper", title="Vegetation data v. density stats")
corrplot(int.cor, type="upper", title="Vegetation data v. intensity stats")

## quick pca
library(ggfortify)
library(ggrepel)  
library(ggalt)
library(plyr)

lid.dat = data.frame(trt=ht.dat$trt, trt.full=ht.dat$trt.full, num=ht.dat$num,
                     mean.ht = ht.dat$MEAN, mean.den = dens.dat$MEAN, mean.int = int.dat$MEAN,
                     sd.ht = ht.dat$STD, sd.den = dens.dat$STD, sd.int = int.dat$STD)
colnames(bm.dat)[4:19] = c("FWD count","CWD area","Snag area","Dead stems","Live woody C stock",
                          "Live stems","SOM","FWD C stock","FWD C:N ratio","Litter C stock", 
                          "Litter C:N ratio","P. arun C stock","H. jap C stock","Mixed biomass C stock",
                          "Total herbaceous C stock","Herb. biomass C:N ratio")
bm.df = right_join(bm.dat, lid.dat, by=c("trt","trt.full","num"))
pca.biomass = prcomp(bm.df[,-which(colnames(bm.df) %in% c("trt","trt.full","num","Total herbaceous C stock"))], 
                     center=T, scale=T)
summary(pca.biomass)

#calculate total variance explained by each principal component
var_explained = pca.biomass$sdev^2 / sum(pca.biomass$sdev^2)

#create scree plot
plot(seq(1:length(var_explained)), var_explained, type="b")

# get loadings matrix
CAloadings <- data.frame(Variables = rownames(pca.biomass$rotation), pca.biomass$rotation)

# add PCA scores to the dataset
bm.df[, c('PC1', 'PC2')] = pca.biomass$x[, 1:2]

# save variable loadings in a separate dataset
rot = as.data.frame(pca.biomass$rotation[, 1:2])
rot$var = rownames(pca.biomass$rotation)

# rescale the loadings to fit nicely within the scatterplot of our data
mult = max(abs(bm.df[, c('PC1', 'PC2')])) / max(abs(rot[, 1:2])) / 2
rot[, 1:2] = rot[, 1:2] * mult

# if there are many variables to plot, you can play with ggrepel 
var1 = paste("PC1", paste("(", 100*round(var_explained[1], 4), "%)", sep=""), sep=" ")
var2 = paste("PC2", paste("(", 100*round(var_explained[2], 4), "%)", sep=""), sep=" ")
ggplot(data = rot, aes(x=0, y=0, xend=PC1, yend=PC2, label=var)) +
        geom_point(data = bm.df, aes(x=PC1, y=PC2, color=factor(trt.full, levels=trt.names)), inherit.aes=FALSE, size=4) +
        geom_segment(color = 'red', arrow = arrow(length = unit(0.03, "npc"))) +
        geom_text_repel(aes(PC1 * 1, PC2 * 1),color = 'red') +
        labs(x = var1, y = var2, color="") + 
        theme(legend.title = element_blank(), legend.text = element_text(size=12), 
              axis.title = element_text(size=14), axis.text = element_text(size=10))
