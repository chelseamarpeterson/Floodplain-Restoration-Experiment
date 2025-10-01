path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo"
setwd(path_to_repo)

library(ggplot2)
library(ggfortify)
library(patchwork)
library(ggrepel)
library(tidyr)
library(dplyr)
library(vegan)

## treatment names
trt.letters = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

## read in all soil data
soil.df = read.csv("Soil_Analysis/Clean_Data/Soil_Data_by_Quadrat_June2023.csv")

# soil data labels
ag.cols = c("fpom","g475mm","g2mm","g250um","g53um","l53um","mwd")
cat.cols = c("k","ca","mg")
meq.cols = c("k.meq","ca.meq","mg.meq")
text.cols = c("Sand","Silt","Clay")

# variables to omit
leave.out = c("ton.percent","pon.percent","texture.class","n.release",
              "volumetric.moisture","coarse.material",
              "k.sat","ca.sat","mg.sat","h.sat","h.meq",
              ag.cols[1:6],meq.cols)
soil.df = soil.df[,-which(colnames(soil.df) %in% leave.out)]

## read in biomass data
bm.df = read.csv("Understory_Analysis/Clean_Data/Biomass_By_Quadrat_and_Year.csv")
bm.df.2022 = bm.df[which(bm.df$year == 2022),-which(colnames(bm.df) == "year")]

# combine biomass and soil data
soil.bm.df = left_join(soil.df, bm.df.2022, 
                       by=c("full.treatment.name","treatment","plot","quadrat"))

# remove reference
soil.bm.df = soil.bm.df[-which(soil.bm.df$treatment == "R"),]

# update column names
colnames(soil.bm.df)[5:36] = c("Temp","Moist","BD","TN","TC","CN","SOC",text.cols,
                               "pH","P","K","Ca","Mg","SOM","NO3","NH4","CEC","POC",
                               "EV","HT","MWD","TIC","MAOC","POC_MAOC",
                               "HB","HL","FWD","MB","PAB","HJB")
 
# select explanatory and response variables
y.cols = c("SOC","POC","MAOC","POC_MAOC")
x.cols = c("Temp","Moist","BD","TN","CN","Sand","Clay",
           "pH","P","NO3","NH4","CEC","EV","HT",
           "MWD","HL","FWD","MB","PAB","HJB")

# scale the data
df.scale = data.frame(scale(soil.bm.df[,c(y.cols,x.cols)], center=T, scale=T))
y.data = df.scale[,y.cols]
x.data = df.scale[,x.cols]

# run the RDA
rda_model <- rda(y.data ~., data = x.data)

# RDA summary data
summary(rda_model)
smry <- summary(rda_model)
df1  <- data.frame(smry$sites[,1:3])       # PC1 and PC2
df1$full.treatment.name = soil.bm.df$full.treatment.name
df2  <- data.frame(smry$species[,1:3])     # loadings for PC1 and PC2
df3  <- data.frame(smry$biplot[,1:3])

# variance explained by each RDA axis
var.rda1 = round(smry$cont$importance["Proportion Explained","RDA1"]*100,1)
var.rda2 = round(smry$cont$importance["Proportion Explained","RDA2"]*100,1)
var.rda1.label = paste("RDA1 (",var.rda1,"%)",sep="")
var.rda2.label = paste("RDA2 (",var.rda2,"%)",sep="")

# make biplot
rda.plot <- ggplot(data=df1) + 
                   geom_point(aes(x=RDA1, y=RDA2, 
                                  color=factor(full.treatment.name, levels=trt.names),
                                  shape=factor(full.treatment.name, levels=trt.names)),
                              size=2) + 
                   guides(shape=guide_legend(title="Treatment"),
                          color=guide_legend(title="Treatment")) +
                   geom_hline(yintercept=0, linetype="dotted") +
                   geom_vline(xintercept=0, linetype="dotted") +
                   coord_fixed() + 
                   geom_segment(data=df3, 
                                aes(x=0, xend=RDA1, y=0, yend=RDA2), 
                                color="blue", 
                                arrow=arrow(length=unit(0.01,"npc"))) +
                   geom_segment(data=df2, 
                                aes(x=0, xend=RDA1, y=0, yend=RDA2), 
                                color="red", 
                                arrow=arrow(length=unit(0.01,"npc"))) +
                   geom_label_repel(data=df3, 
                                    aes(x=RDA1, y=RDA2, label=rownames(df3)),
                                    color="blue",fill="transparent",
                                    alpha=0.7,size=2.5) +
                   geom_label_repel(data=df2,
                                    aes(x=RDA1, y=RDA2, label=rownames(df2)),
                                    color="red", fill="transparent",
                                    alpha=0.7,size=3) +
                   labs(x=var.rda1.label,y=var.rda2.label) + 
                   theme(text=element_text(size=12))
rda.plot

#df3$rows = rownames(df3)
#var.sort1 = sort(df3[,"RDA1"], decreasing=T, index.return=T)
#df3[var.sort1$ix, c("RDA1","rows")]

#var.sort2 = sort(df3[,"RDA2"], decreasing=T, index.return=T)
#df3[var.sort2$ix, c("RDA2","rows")]

# write plot
ggsave("Figures/Figure5_Soil_Vegetation_RDA_HTElev.jpeg", 
       plot=rda.plot, width=26, height=20, units="cm", dpi=600)

# compute Pearsons r^2 to confirm inference about correlations
# among SOC, MAOC, POC, and POC:MAOC ratio
#cor(df.scale$SOC, df.scale$MAOC, method = "pearson")
#cor(df.scale$SOC, df.scale$POC, method = "pearson")
#cor(df.scale$SOC, df.scale$POC_MAOC, method = "pearson")

