path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo"
setwd(path_to_repo)

library(readxl)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(rethinking)
library(allodb)

# treatment names
trt.letters = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")
n.t = length(trt.letters)

################################################################################
## statistical analyses on live and dead C stock components by plot

# read in biomass C stock estimates
veg.soil.df = read.csv("Tree_Analysis/Clean_Data/Clean_Veg_Soil_Cstocks_Richness.csv")

# variable names
vars = colnames(veg.soil.df)[4:40]
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

# fit simple models
set.seed
m.all = list()
for (i in 1:n.v) {
  v = vars[i]
  print(v)
  m.all[[v]] = ulam(alist(y ~ dnorm(mu, sigma),
                          log(mu) <- a[tid],
                          a[tid] ~ dnorm(0, 1),
                          sigma ~ dexp(1)),
                    data = all.df[[v]], chains=5)
}

# MCMC checks
for (i in 1:n.v) { traceplot(m.all[[vars[i]]]) }
for (i in 1:n.v) { trankplot(m.all[[vars[i]]]) }

# variable in new order with labels
var.labels = c("Fine woody debris (2.5-7.5 cm)",
               "Coarse woody debris (\u2265 7.6 cm)",
               "Standing dead trees (\u2265 2.5 cm)",
               "Live F. pennsylvanica (\u2265 2.5 cm)",
               "Standing dead F. pennsylvanica (\u2265 2.5 cm)",
               "Dead stems (< 2.5 cm)",
               "Aboveground live stems (< 2.5 cm)",
               "Belowground live stems (< 2.5 cm)",
               "Live stems (< 2.5 cm)",
               "Aboveground live trees (\u2265 2.5 cm)",
               "Belowground live trees (\u2265 2.5 cm)",
               "Live trees (\u2265 2.5 cm)",
               "Herbaceous biomass",
               "Mixed biomass",
               "P. arundinacea biomass",
               "H. japonicus biomass",
               "Herbaceous litter",
               "Fine woody debris (< 2.5 cm)",
               "Fine woody debris C:N ratio",
               "Herbaceous litter C:N ratio",
               "Herbaceous biomass C:N ratio",
               "MAOC",
               "POC",
               "SOC",
               "TIC",
               "TC",
               "Aboveground woody biomass",
               "Belowground woody biomass",
               "Total live vegetation",
               "Total dead vegetation",
               "Total vegetation",
               "Total ecosystem",
               "Live trees w/ F. pennsylvanica",
               "Standing dead trees w/o F. pennsylvanica",
               "Herbaceous species",
               "Tree species",
               "Total richness")
length(vars)
length(var.labels)

# extract posterior samples and create dataframes
all.samples = list()
all.a = list()
a.df = data.frame(matrix(nrow=0, ncol=7))
colnames(a.df) = c("var",trt.letters)
for (i in 1:n.v) {
  v1 = vars[i]
  v2 = var.labels[i]
  
  # store posterior samples in lists
  all.samples[[v2]] = extract.samples(m.all[[v1]])
  all.a[[v2]] = exp(all.samples[[v2]]$a) * all.mean[[v1]]
  
  # add a's to dataframe
  a.row = data.frame(all.a[[v2]])
  colnames(a.row) = trt.letters
  a.row$var = var.labels[i]
  a.df = rbind(a.df, a.row[,c("var",trt.letters)])
}

# create df with means and HPDIs by treatment and species
post.df = data.frame(matrix(nrow=n.v*n.t, ncol=8))
colnames(post.df) = c("tid","trt","var","mean","5","95","25","75")
k = 1
for (j in 1:length(vars)) {
  for (i in 1:n.t) {
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
write.csv(post.df,
          "Tree_Analysis/Posteriors/CarbonRichness_Means_Intervals_5Chains.csv",
          row.names=F)

# plot continuous distibutions
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

################################################################################
### statistical analysis on carbon stock data by species

## coarse woody debris C stocks

# read in CWD data
cwd.melt = read.csv("Tree_Analysis/Clean_Data/CWD_Carbon_Stocks_By_Species.csv", header=T)

# convert treatment & species into indices
cwd.melt$tid = as.numeric(factor(cwd.melt$trt))
cwd.melt$spid = as.numeric(factor(cwd.melt$species))

# make simplified data list
mean.cwd = mean(cwd.melt$cwd.carbon)
d.cwd.sp = list(tid = cwd.melt$tid, 
                spid = cwd.melt$spid,
                y = cwd.melt$cwd.carbon/mean.cwd) 

# fit simple model
set.seed(12)
m.cwd.sp = ulam(alist(y ~ dnorm(mu, sigma),
                      log(mu) <- a[tid, spid],
                      matrix[tid, spid]: a ~ dnorm(0, 1),
                      sigma ~ dexp(1)),
                data = d.cwd.sp, chains=5)

# check model
#traceplot(m.cwd.sp)
#trankplot(m.cwd.sp)

# get posterior distributions
post.cwd.sp = extract.samples(m.cwd.sp)
a.cwd = exp(post.cwd.sp$a)*mean.cwd
spp.cwd = levels(factor(cwd.melt$species))
n.spp.cwd = length(spp.cwd)

# create df with means and HPDIs by treatment and species
post.df.cwd = data.frame(matrix(nrow=n.t*n.spp.cwd, ncol=9))
colnames(post.df.cwd) = c("tid","spid","trt","spp","mean","5","95","25","75")
k = 1
for (i in 1:n.t) {
  for (j in 1:n.spp.cwd) {
    a.mean = mean(a.cwd[,i,j], 2)
    a.hpdi.90 = HPDI(a.cwd[,i,j], prob=0.90)
    a.hpdi.50 = HPDI(a.cwd[,i,j], prob=0.50)
    post.df.cwd[k,"tid"] = i
    post.df.cwd[k,"spid"] = j
    post.df.cwd[k,"trt"] = trt.names[i]
    post.df.cwd[k,"spp"] = spp.cwd[j]
    post.df.cwd[k,"mean"] = a.mean
    post.df.cwd[k,"5"] = a.hpdi.90[1]
    post.df.cwd[k,"95"] = a.hpdi.90[2]
    post.df.cwd[k,"25"] = a.hpdi.50[1]
    post.df.cwd[k,"75"] = a.hpdi.50[2]
    k = k + 1
  }
}

## snag C stocks

# read in snag data
snag.melt = read.csv("Tree_Analysis/Clean_Data/Snag_Carbon_Stocks_By_Species.csv", header=T)

# convert treatment and species into id
snag.melt$tid = as.numeric(factor(snag.melt$trt))
snag.melt$spid = as.numeric(factor(snag.melt$species))

# make simplified data list
snag.mean = mean(snag.melt$snag.carbonmin)
d.snag.sp = list(tid = snag.melt$tid, 
                 spid = snag.melt$spid,
                 y = snag.melt$snag.carbonmin/snag.mean) 

# fit simple model
set.seed(12)
m.snag.sp = ulam(alist(y ~ dnorm(mu, sigma),
                       log(mu) <- a[tid, spid],
                       matrix[tid, spid]:a ~ dnorm(0, 1),
                       sigma ~ dexp(1)),
                 data = d.snag.sp, chains=5)

# check model
#traceplot(m.snag.sp)
#trankplot(m.snag.sp)

# get posterior distributions
post.snag.sp = extract.samples(m.snag.sp)
a.snag = exp(post.snag.sp$a)*snag.mean
spp.snag = levels(factor(snag.melt$species))
n.spp.snag = length(spp.snag)

# create df with means and HPDIs by treatment and species
post.df.snag = data.frame(matrix(nrow=n.t*n.spp.snag, ncol=9))
colnames(post.df.snag) = c("tid","spid","trt","spp","mean","5","95","25","75")
k = 1
for (i in 1:n.t) {
  for (j in 1:n.spp.snag) {
    a.mean = mean(a.snag[,i,j], 2)
    a.hpdi.90 = HPDI(a.snag[,i,j], prob=0.90)
    a.hpdi.50 = HPDI(a.snag[,i,j], prob=0.50)
    post.df.snag[k,"tid"] = i
    post.df.snag[k,"spid"] = j
    post.df.snag[k,"trt"] = trt.names[i]
    post.df.snag[k,"spp"] = spp.snag[j]
    post.df.snag[k,"mean"] = a.mean
    post.df.snag[k,"5"] = a.hpdi.90[1]
    post.df.snag[k,"95"] = a.hpdi.90[2]
    post.df.snag[k,"25"] = a.hpdi.50[1]
    post.df.snag[k,"75"] = a.hpdi.50[2]
    k = k + 1
  }
}

## large woody biomass C stocks

# prep dataframe
treeC.dat = read.csv("Tree_Analysis/Clean_Data/WoodyBiomass_C_Stocks_By_Species.csv")
treeC.dat.2022 = treeC.dat[which(treeC.dat$year == 2022), 
                             c("trt","num","genus","species","bmC_ha3")]

# make full species name column
treeC.dat.2022 = unite(treeC.dat.2022, col='spp', c('genus','species'), sep=' ')

# make wide df
treeC.wide = pivot_wider(treeC.dat.2022, 
                          id_cols = c(trt, num),
                          names_from = spp, 
                          values_from = bmC_ha3)
treeC.wide[is.na(treeC.wide)] = 0

# reshape dataframe
treeC.melt = melt(treeC.wide, id.vars=c("trt","num"), variable.name="spp")

# convert treatments and species into indices
treeC.melt$tid = as.numeric(factor(treeC.melt$trt))
treeC.melt$spid = as.numeric(factor(treeC.melt$spp))

# make simplified data list
mean.treeC = mean(treeC.melt$value)
d.treeC = list(tid = treeC.melt$tid, 
                spid = treeC.melt$spid,
                y = treeC.melt$value/mean.treeC)

# fit simple model
set.seed(12)
m.treeC = ulam(alist(y ~ dnorm(mu, sigma),
                      log(mu) <- a[tid, spid],
                      matrix[tid, spid]: a ~ dnorm(0, 1),
                      sigma ~ dexp(1)),
                data=d.treeC, chains=5)

# model checks
#traceplot(m.bmC)
#trankplot(m.bmC)

# get posterior distributions
post.treeC = extract.samples(m.treeC)
a.scale = exp(post.treeC$a) * mean.treeC
spp.treeC = levels(factor(sp.C.melt$spp))
n.spp.treeC = length(spp.treeC)

# create df with means and HPDIs by treatment and species
post.df.treeC = data.frame(matrix(nrow=n.t*n.spp.treeC, ncol=9))
colnames(post.df.treeC) = c("tid","spid","trt","spp","mean","5","95","25","75")
k = 1
for (i in 1:n.t) {
  for (j in 1:n.spp.bmC) {
    a.mean = mean(a.scale[,i,j], 2)
    a.sd = sd(a.scale[,i,j])
    a.hpdi.90 = HPDI(a.scale[,i,j], prob=0.90)
    a.hpdi.50 = HPDI(a.scale[,i,j], prob=0.50)
    post.df.treeC[k,"tid"] = i
    post.df.treeC[k,"spid"] = j
    post.df.treeC[k,"trt"] = trt.names[i]
    post.df.treeC[k,"spp"] = spp.treeC[j]
    post.df.treeC[k,"mean"] = a.mean
    post.df.treeC[k,"5"] = a.hpdi.90[1]
    post.df.treeC[k,"95"] = a.hpdi.90[2]
    post.df.treeC[k,"25"] = a.hpdi.50[1]
    post.df.treeC[k,"75"] = a.hpdi.50[2]
    k = k + 1
  }
}



## herbaceous biomass C stocks

# read in biomass C stock estimates
herbC.dat = read.csv("Tree_Analysis/Clean_Data/Clean_Veg_Soil_Cstocks_Richness.csv")

# isolate biomass data for each species
herbC.vars = c("PAB.C.Mg.ha","HJB.C.Mg.ha","MB.C.Mg.ha")
herbC.spp = c("Phalaris arundinacea","Humulus japonicus","Mixed community")
herbC.spp.df = herbC.dat[,c("trt","trt.full","num",herbC.vars)]
colnames(herbC.spp.df)[4:6] = herbC.spp

# melt dataframe
herbC.spp.melt = melt(herbC.spp.df, 
                      id.vars=c("trt","trt.full","num"),
                      variable.name="species")

# convert treatmentinto indices
herbC.spp.melt$tid = as.numeric(factor(herbC.spp.melt$trt))
herbC.spp.melt$spid = as.numeric(factor(herbC.spp.melt$species))

# make simple dataframe
herbC.mean = mean(herbC.spp.melt$value)
d.herbC = list(tid = as.integer(herbC.spp.melt$tid),
               spid = as.integer(herbC.spp.melt$spid),
               y = herbC.spp.melt$value/herbC.mean)

# fit model
set.seed(12)
m.herbC.sp = ulam(alist(y ~ dnorm(mu, sigma),
                     log(mu) <- a[tid, spid],
                     matrix[tid, spid]:a ~ dnorm(0, 1),
                     sigma ~ dexp(1)),
               data=d.herbC, chains=5)

# check model
#traceplot(m.herbC.sp)
#trankplot(m.herbC.sp)

# get posterior distributions
post.herbC.sp = extract.samples(m.herbC.sp)
a.herbC.sp = exp(post.herbC.sp$a) * herbC.mean
spp.herbC = levels(factor(herbC.spp.melt$species))
n.spp.herbC = length(spp.herbC)

# create df with means and HPDIs by treatment and species
post.df.herbC = data.frame(matrix(nrow=n.t*n.spp.herbC, ncol=9))
colnames(post.df.herbC) = c("tid","spid","trt","spp","mean","5","95","25","75")
k = 1
for (i in 1:n.t) {
  for (j in 1:n.spp.herbC) {
    a.mean = mean(a.herbC.sp[,i,j], 2)
    a.hpdi.90 = HPDI(a.herbC.sp[,i,j], prob=0.90)
    a.hpdi.50 = HPDI(a.herbC.sp[,i,j], prob=0.50)
    post.df.herbC[k,"tid"] = i
    post.df.herbC[k,"spid"] = j
    post.df.herbC[k,"trt"] = trt.names[i]
    post.df.herbC[k,"spp"] = spp.herbC[j]
    post.df.herbC[k,"mean"] = a.mean
    post.df.herbC[k,"5"] = a.hpdi.90[1]
    post.df.herbC[k,"95"] = a.hpdi.90[2]
    post.df.herbC[k,"25"] = a.hpdi.50[1]
    post.df.herbC[k,"75"] = a.hpdi.50[2]
    k = k + 1
  }
}


# combine CWD, snag, woody C, and herbaocues stock posteriors. then write to file
post.df.snag$type = "dead.trees"
post.df.treeC$type = "live.trees"
post.df.cwd$type = "coarse.woody.debris"
post.df.herbC$type = "herbaceous.biomass"
post.df.all = rbind(post.df.snag, post.df.treeC, post.df.cwd,post.df.herbC)
write.csv(post.df.all,
          "Tree_Analysis/Posteriors/VegetationCarbonStocks_Species_Means_Intervals_5Chains.csv",
          row.names=F)


