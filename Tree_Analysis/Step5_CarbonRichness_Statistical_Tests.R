path_to_tree_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Public-Repo/Tree_Analysis"
path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Public-Repo"
setwd(path_to_tree_folder)

library(readxl)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(rethinking)
library(allodb)

# treatment names
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")
n.trt = length(trts)

################################################################################
## statistical analyses on live and dead C stock components by plot ----

# read in biomass C stock estimates
veg.soil.df = read.csv("Clean_Data/Clean_Veg_Soil_Cstocks_Richness.csv")

# estimate aggregate C stocks in different pools
veg.soil.df$total.live.carbon =  veg.soil.df$live.woody.c.stock + veg.soil.df$live.stem.carbon + veg.soil.df$total.C.Mg.ha
veg.soil.df$total.dead.carbon =  veg.soil.df$snag.carbonmin + veg.soil.df$dead.stem.carbonmin + veg.soil.df$cwd.carbon + veg.soil.df$int.fwd.carbon + veg.soil.df$FWD.C.Mg.ha
veg.soil.df$total.veg.carbon = veg.soil.df$total.live.carbon + veg.soil.df$total.dead.carbon
veg.soil.df$total.ecosystem.carbon = veg.soil.df$TC + veg.soil.df$total.veg.carbon
veg.soil.df$n.total = veg.soil.df$N.herb + veg.soil.df$N.tree

# variable names
vars = colnames(veg.soil.df)[4:30]
n.v = length(vars)

# convert treatments into indices
veg.soil.df$tid = as.numeric(factor(veg.soil.df$trt))

# make simplified data lists
all.df = list()
all.mean = list()
all.sd = list()
for (i in 1:n.v) {
  v = vars[i]
  all.df[[v]] = list(tid = veg.soil.df$tid, y = veg.soil.df[,v]/mean(veg.soil.df[,v]))
  all.mean[[v]] = mean(veg.soil.df[,v])
  all.sd[[v]] = sd(veg.soil.df[,v])
}

# look at histograms
#par(mfrow=c(3,4))
#for (i in 1:n.v) {
#  v = vars[i]
  #hist(all.df[[v]]$y)
#  hist(log(all.df[[v]]$y), main="")
#}

# fit simple models
set.seed(12)
m.all = list()
for (i in 1:n.v) {
  v = vars[i]
  m.all[[v]] = ulam(alist(y ~ dnorm(mu, sigma),
                          log(mu) <- a[tid],
                          a[tid] ~ dnorm(0, 1),
                          sigma ~ dexp(1)),
                    data = all.df[[v]], chains=5)
}

# MCMC checks
for (i in 1:n.v) { traceplot(m.all[[vars[i]]]) }
for (i in 1:n.v) { trankplot(m.all[[vars[i]]]) }

# get posterior distributions
post.all = list()
for (i in 1:n.v) {
  v = vars[[i]]
  post.all[[v]] = extract.samples(m.all[[v]])
}

# variable in new order with labels
var.order = c("n.total","N.herb","N.tree",
              "MAOM.C","POM.C","SOC","TIC","TC",
              "int.fwd.carbon","cwd.carbon",
              "snag.carbonmin","dead.stem.carbonmin",
              "live.woody.c.stock","live.stem.carbon",
              "FWD.C.Mg.ha","HL.C.Mg.ha","total.C.Mg.ha",
              "total.live.carbon","total.dead.carbon",
              "total.veg.carbon","total.ecosystem.carbon")
var.labels = c("Total richness","Herbaceous species","Tree species",
               "MAOM-C (< 53 \U00B5m)","POM-C (\u2265 53 \U00B5m)",
               "TOC","Carbonates","TC",
               "Fine woody debris (2.5-7.5 cm)","Coarse woody debris (\u2265 7.6 cm)",
               "Standing dead trees (\u2265 2.5 cm)","Dead stems (< 2.5 cm)",
               "Live trees (\u2265 2.5 cm)","Live stems (< 2.5 cm)",
               "Fine woody debris (< 2.5 cm)","Herbaceous litter","Herbaceous biomass",
               "Total live vegetation","Total dead vegetation",
               "Total vegetation","Total ecosystem")

# extract posterior samples and create dataframes
all.samples = list()
all.a = list()
a.df = data.frame(matrix(nrow=0, ncol=7))
colnames(a.df) = c("var",trts)
n.v.short = length(var.order)
for (i in 1:n.v.short) {
  v1 = var.order[i]
  v2 = var.labels[i]
  
  # store posterior samples in lists
  all.samples[[v2]] = extract.samples(m.all[[v1]])
  all.a[[v2]] = exp(all.samples[[v2]]$a) * all.mean[[v1]]
  
  # add a's to dataframe
  a.row = data.frame(all.a[[v2]])
  colnames(a.row) = trts
  a.row$var = var.labels[i]
  a.df = rbind(a.df, a.row[,c("var",trts)])
}

# create df with means and HPDIs by treatment and species
post.df = data.frame(matrix(nrow=n.v.short*n.trt, ncol=8))
colnames(post.df) = c("tid","trt","var","mean","5","95","25","75")
k = 1
for (j in 1:n.v.short) {
  for (i in 1:n.trt) {
    v = var.labels[j]
    a.v = all.a[[v]]
    a.mean = mean(a.v[,i], 2)
    a.hpdi.90 = HPDI(a.v[,i], prob=0.90)
    a.hpdi.50 = HPDI(a.v[,i], prob=0.50)
    post.df[k,"tid"] = i
    post.df[k,"trt"] = trt.names[i]
    post.df[k,"var"] = var.labels[j]
    post.df[k,"mean"] = a.mean
    post.df[k,"5"] = a.hpdi.90[1]
    post.df[k,"95"] = a.hpdi.90[2]
    post.df[k,"25"] = a.hpdi.50[1]
    post.df[k,"75"] = a.hpdi.50[2]
    k = k + 1
  }
}

# write posterior HPDIs and means to file
write.csv(post.df,"Posteriors/CarbonRichness_Means_Intervals_5Chains.csv",
          row.names=F)


# melt posterior samples into dataframe
new.df = all.a[["Total richness"]]
colnames(new.df) = trt.names
all.melt = melt(new.df, id.vars=trt.names)
all.melt$var = "Total richness"
all.melt = all.melt[,2:4]
colnames(all.melt) = c("trt","value","var")
for (v in var.labels[2:n.v.short]) {
  new.df = all.a[[v]]
  colnames(new.df) = trt.names
  melt.df = melt(new.df, id.vars=trt.names) 
  melt.df$var = v
  melt.df = melt.df[,2:4]
  colnames(melt.df) = c("trt","value","var")
  all.melt = rbind(all.melt, melt.df)
}
all.melt$var = factor(all.melt$var, levels=var.labels)

# plot posterior distributions
ggplot(all.melt, aes(x=value, fill=trt)) + 
  geom_density(alpha=0.5, linewidth=0.2) + 
  facet_wrap(.~var, scales = "free") +
  labs(y="Density",x="") + 
  theme(legend.title=element_blank())
