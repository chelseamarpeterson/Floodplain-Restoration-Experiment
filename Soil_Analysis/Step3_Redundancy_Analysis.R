path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo"
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

# combine biomass and soil data
all.dat = left_join(soil.dat[,-which(colnames(soil.dat) %in% c(meq.cols,"volumetric.moisture","ton","pon"))], bm.dat.2022, 
                    by=c("trt","trt.full","num","quad"))

# remove reference
df = all.dat[-which(all.dat$trt == "R"),]

# update column names
colnames(df)[6:42] = c("Temp","Moisture","BD","CM","TN","TC","CN","SOC",text.cols,
                       "pH","P","K","Ca","Mg","SOM","NO3","NH4","CEC","POC",
                       ag.cols,"TIC","MAOC","POC_MAOC","HB","HL","FWD","MB","PAB","HJB")
################################################################################
# redundancy analysis
library(vegan)

quad.df = read.csv("Soil_Analysis/Raw_Data/Quadrat_GIS_Data_PolygonMean_WGS1984Aux.csv", header=T)
colnames(quad.df) = c("trt","num","quad","trt_num_quad","EV","HT","shape_length","shape_area")
quad.df = quad.df[,c("trt","num","quad","EV","HT")]
df$num = as.integer(df$num)
quad.df$num = as.integer(quad.df$num)
quad.df$quad = paste("Q", quad.df$quad, sep="")
df = left_join(df, quad.df, by=c("trt","num"))
y.cols = c("SOC","POC","MAOC","POC_MAOC")
x.cols = c("TN","CN","MWD","CEC","Sand","Clay","pH","P","NO3",
           "NH4","BD","Moisture","Temp","HL","FWD","MB","PAB","HJB",
           "EV","HT")
df.scale = data.frame(scale(df[,c(y.cols,x.cols)], center=T, scale=T))
y.data = df.scale[,y.cols]
x.data = df.scale[,x.cols]
rda_model <- rda(y.data ~., data = x.data)
#plot(rda_model)
summary(rda_model)

smry <- summary(rda_model)
df1  <- data.frame(smry$sites[,1:3])       # PC1 and PC2
df1$trt.full = df$trt.full
df2  <- data.frame(smry$species[,1:3])     # loadings for PC1 and PC2
df3  <- data.frame(smry$biplot[,1:3])

rda.plot1 <- ggplot(data=df1) + 
                   geom_point(aes(x=RDA1, y=RDA2, 
                                  color=factor(trt.full, levels=trt.names),
                                  shape=factor(trt.full, levels=trt.names)),
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
                   labs(x="RDA1 (30.4%)",y="RDA2 (17.9%)") + 
                   theme(text=element_text(size=12))
rda.plot1

df3$rows = rownames(df3)
var.sort1 = sort(df3[,"RDA1"], decreasing=T, index.return=T)
df3[var.sort1$ix, c("RDA1","rows")]

var.sort2 = sort(df3[,"RDA2"], decreasing=T, index.return=T)
df3[var.sort2$ix, c("RDA2","rows")]

# write plot
ggsave("Figures/Figure5_Soil_Vegetation_RDA_HTElev.jpeg", 
       plot=rda.plot1, width=26, height=20, units="cm", dpi=600)

# compute Pearsons r^2 to confirm inference about correlations
# among SOC, MAOC, POC, and POC:MAOC ratio
cor(df.scale$SOC, df.scale$MAOC, method = "pearson")
cor(df.scale$SOC, df.scale$POC, method = "pearson")
cor(df.scale$SOC, df.scale$POC_MAOC, method = "pearson")

########################################################################
# Determine which variables are most related to each SOC variable
library(rethinking)

# make data lists
out = c("plot","trt","trt.full","num","quad","volumetric.moisture",
        "k.meq","ca.meq","mg.meq","h.meq","pon","fpom",
        "som","poc","maoc","tic")
c.vars = c("toc")
x.vars = colnames(all.dat)[-which(colnames(all.dat) %in% c(c.vars,out))]
n.c = length(c.vars)
n.x = length(x.vars)
c.lists = list()
v.means = list()
for (i in 1:n.c) {
  c = c.vars[i]
  v.means[[c]] = mean(all.dat[,c])
  c.lists[[c]] = list()
  for (j in 1:n.x) {
    x = x.vars[j]
    v.means[[x]] = mean(all.dat[,x])
    c.lists[[c]][[x]] = list(y=all.dat[,c]/v.means[[c]],
                             x=all.dat[,x]/v.means[[x]])
  }
}

# run models
c.models = list()
for (i in 1:n.c) {
  c = c.vars[i]
  c.models[[c]] = list()
  for (j in 1:n.x) {
    x = x.vars[j]
    c.models[[c]][[x]] = ulam(alist(y ~ normal(mu, sigma),
                               log(mu) <- a*x,
                               a ~ dnorm(0, 1),
                               sigma ~ dexp(1)),
                         data=c.lists[[c]][[x]], chains=1, log_lik=T)
  }
}

# compare waics across models
waic.df = data.frame(matrix(nrow=0,ncol=8))
colnames(waic.df) = c("WAIC","SE","dWAIC","dSE","pWAIC","weight","y.var","x.var")
for (i in 1:n.c) {
  c = c.vars[i]
  waic.c = compare(c.models[[c]][["temperature"]],
                   c.models[[c]][["gravimetric.moisture"]],
                   c.models[[c]][["bulk.density"]],
                   c.models[[c]][["bulk.n"]],
                   c.models[[c]][["bulk.cn"]],
                   c.models[[c]][["sand"]],
                   c.models[[c]][["silt"]],
                   c.models[[c]][["clay"]],
                   c.models[[c]][["ph"]],
                   c.models[[c]][["p"]],
                   c.models[[c]][["no3"]],
                   c.models[[c]][["nh4"]],
                   c.models[[c]][["cec"]],
                   c.models[[c]][["poc_maoc"]],
                   c.models[[c]][["mwd"]],
                   c.models[[c]][["Herbaceous.litter"]],
                   c.models[[c]][["Fine.woody.debris"]],
                   func=WAIC)
  waic.c$y.var = c
  waic.c$x.var = rownames(waic.c)
  waic.df = rbind(waic.df, waic.c)
}

# rename model column
waic.df$model = rep(0,nrow(waic.df))
waic.df$coef = rep(0,nrow(waic.df))
for (i in 1:n.x) {
  x = x.vars[i]
  x.long.label = paste("c.models[[c]][[\"", x, "\"]]",sep="")
  x.short.label = x
  waic.df$model[which(waic.df$x.var == x.long.label)] = x.short.label
  #sample.x = samples(c.models[["toc"]][[x]])
  #waic.df$coef[which(waic.df$x.var == x.long.label)] = c.models[["toc"]][[x]]
}

model.sort = sort(waic.df[waic.df$y.var == "toc","WAIC"],index.return=T)
model.order = waic.df[waic.df$y.var == "toc","model"][model.sort$ix]
ggplot(waic.df, 
       aes(x=WAIC, y=factor(model, levels=model.order))) + 
       geom_point() + 
       facet_wrap(.~y.var, scales="free_x")

write.csv(waic.df,"Soil_Analysis/Posteriors/Univariate_Linear_Model_WAIC_Comparison.csv",
          row.names=F)
