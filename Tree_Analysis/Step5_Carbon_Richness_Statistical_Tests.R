path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch2_Floodplain_Experiment/Floodplain-Experiment-Repo"
setwd(path_to_repo)

library(readxl)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(rethinking)
library(allodb)
library(brms)

# treatment names
trt.letters = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")
n.t = length(trt.letters)

################################################################################
## statistical analyses on live and dead C stock components by plot

# read in biomass C stock estimates
stock.richness.df = read.csv("Tree_Analysis/Clean_Data_By_Plot/Clean_Veg_Soil_C_Stocks_Richness_by_Plot.csv")

# variable names
vars = read.csv("Tree_Analysis/Clean_Data_By_Plot/carbon.richness.metadata.csv")[,"variable"]
var.labels = read.csv("Tree_Analysis/Clean_Data_By_Plot/carbon.richness.metadata.csv", stringsAsFactors=F)[,"label"]
n.v = length(vars)

# make simplified data lists
all.lists = list()
all.mean = list()
all.sd = list()
for (i in 1:n.v) {
  v = vars[i]
  all.mean[[v]] = mean(stock.richness.df[,v])
  all.sd[[v]] = sd(stock.richness.df[,v])
  all.lists[[v]] = list(treatment = factor(stock.richness.df$treatment), 
                        strip = factor(stock.richness.df$strip), 
                        plot = factor(stock.richness.df$plot),
                        y = stock.richness.df[,v]/all.mean[[v]])
  hist(all.lists[[v]]$y)
}

# fit simple models
seed = 3141; n.iter=5000; n.chain=5
stock.model.list = list()
models = c("simple","plot.random","strip.plot.random")
model.labels = c("Fixed effects only","Plot random effects","Strip and plot random effects")
n.m = length(models)
for (i in 2:2) {
  v.i = vars[i]
  stock.model.list[[v.i]] = list()
  for (j in 2:2) {
    m.j = models[j]
    if (m.j == "simple") {
      model.fit.ij = brm(data = all.lists[[v.i]], 
                         y ~ treatment,
                         family=gaussian(link="log"),
                         prior=prior(normal(0,1), class=b),
                         chains=10, seed=seed, iter=10000,
                         control = list(adapt_delta = 0.99,
                                        max_treedepth = 12))
    } else if (m.j == "plot.random") {
      model.fit.ij = brm(data = all.lists[[v.i]], 
                         y ~ treatment + (1|treatment:plot),
                         family=gaussian(link="log"),
                         prior=prior(normal(0,2), class=b),
                         chains=10, seed=seed, iter=10000,
                         control = list(adapt_delta = 0.99,
                                        max_treedepth = 12))
    l;6} else if (m.j == "strip.plot.random") {
      model.fit.ij = brm(data = all.lists[[v.i]], 
                         y ~ treatment + (1|treatment:strip) + (1|treatment:strip:plot),
                         family=gaussian(link="log"),
                         prior=prior(normal(0,2), class=b),
                         chains=10, seed=seed, iter=10000,
                         control = list(adapt_delta = 0.99,
                                        max_treedepth = 12))
    }
    model.fit.ij = add_criterion(model.fit.ij,"loo")
    model.fit.ij = add_criterion(model.fit.ij,"waic")
    stock.model.list[[v.i]][[m.j]] = model.fit.ij
  }
}

summary(model.fit.ij)
rhat(model.fit.ij)


# save Rhat statistic for each variable and model
#rhat.df = data.frame(matrix(nrow=0,ncol=7))
colnames(rhat.df) = c("variable","variable.label",
                      "model","model.label",
                      "treatment","full.treatment.name","Rhat")
#for (i in 1:n.v) {
for (i in 1:1) {
  v.i = vars[i]
  for (j in 1:n.m) {
    m.j = models[j]
    model.sum.ij = summary(stock.model.list[[v.i]][[m.j]])
    model.df.ij = data.frame(matrix(nrow=6,ncol=7))
    colnames(model.df.ij) = c("variable","variable.label",
                              "model","model.label",
                              "treatment","full.treatment.name","Rhat")
    model.df.ij$variable = v.i
    model.df.ij$variable.label = var.labels[i]
    model.df.ij$model = m.j
    model.df.ij$model.label = model.labels[j]
    model.df.ij$treatment = trt.letters
    model.df.ij$full.treatment.name = trt.names
    model.df.ij$Rhat = model.sum.ij$fixed$Rhat
    rhat.df = rbind(rhat.df, model.df.ij)
  }
}
rhat.df$model.label = rep(0, nrow(rhat.df))
for (i in 1:n.m) { rhat.df$model.label[which(rhat.df$model == models[i])] = model.labels[i] }
#write.csv(rhat.df, "Tree_Analysis/Posteriors/Stock_Rhat_Statistic.csv", row.names=F)

# compare models for each variale with WAIC, LOO, and kfold
criteria = c("waic","loo")
criteria.labels = c("WAIC","LOO")
n.c = length(criteria)
comp.df = data.frame(matrix(nrow=0, ncol=13))
colnames(comp.df) = c("elpd_diff","se_diff","elpd","se_elpd","p","se_p","ic","se_ic",
                      "variable","variable.label","criterion","model","best.model")
for (i in 1:n.v) {
  v.i = vars[i]
  for (j in 1:n.c) {
    m1 = stock.model.list[[v.i]][[models[1]]]
    m2 = stock.model.list[[v.i]][[models[2]]]
    m3 = stock.model.list[[v.i]][[models[3]]]
    comp.ij = data.frame(loo_compare(m1, m2, m3, criterion = criteria[j],
                                     model_names=models))
    colnames(comp.ij) = c("elpd_diff","se_diff","elpd","se_elpd","p","se_p","ic","se_ic")
    comp.ij$variable = v.i
    comp.ij$variable.label = var.labels[i]
    comp.ij$criterion = criteria.labels[j]
    comp.ij$model = row.names(comp.ij)
    comp.ij$best.model = c(1,0,0)
    comp.df = rbind(comp.df, comp.ij)
  }
}

# write information criterion comparisons to file
comp.df$model.label = rep(0, nrow(comp.df))
for (i in 1:n.m) { comp.df$model.label[which(comp.df$model == models[i])] = model.labels[i] }
comp.df$best.model = c("No","Yes")[comp.df$best.model + 1]
write.csv(comp.df, "Tree_Analysis/Posteriors/Stock_Model_Information_Criteria.csv", row.names=F)

# get posterior intervals and write to file
df.int = data.frame(matrix(nrow=0, ncol=13))
int.cols = c("model","model.label",
             "variable","variable.label","variable.mean",
             "treatment","full.treatment.name",
             "posterior.mean","posterior.se",
             "5","95","25","75")
colnames(df.int) = int.cols
for (i in 1:n.v) {
  v.i = vars[i]
  for (j in 1:n.m) {
    m.j = models[j]
    fit.ij = stock.model.list[[v.i]][[m.j]]
    df.ij = data.frame(matrix(nrow=n.t, ncol=13))
    colnames(df.ij) = int.cols
    v = vars[i]
    hdi.90 = hdi(fit.ij, ci=0.90)
    hdi.50 = hdi(fit.ij, ci=0.50)
    df.ij$model = models[j]
    df.ij$model.label = model.labels[j]
    df.ij$variable = v.i
    df.ij$variable.label = var.labels[i]
    df.ij$variable.mean = all.means[[v]]
    df.ij$treatment = trt.letters
    df.ij$full.treatment.name = trt.names
    df.ij$posterior.mean = fixef(fit.ij)[,"Estimate"]
    df.ij$posterior.se = fixef(fit.ij)[,"Est.Error"]
    df.ij[1:n.t,"5"] = hdi.90[,"CI_low"]
    df.ij[1:n.t,"95"] = hdi.90[,"CI_high"]
    df.ij[1:n.t,"25"] = hdi.50[,"CI_low"]
    df.ij[1:n.t,"75"] = hdi.50[,"CI_high"]
    df.int = rbind(df.int, df.ij)
  }
}
df.int$model.label = rep(0, nrow(df.int))
for (i in 1:n.m) { df.int$model.label[which(df.int$model == models[i])] = model.labels[i] }
write.csv(df.int, "Tree_Analysis/Posteriors/Carbon_Stocks_Richness_Means_Intervals_5Chains_BRMS.csv", row.names=F)

################################################################################
### statistical analyses on carbon stock data by species

## coarse woody debris C stocks

# read in CWD data
cwd.melt = read.csv("Tree_Analysis/Clean_Data_By_Species/CWD_Carbon_Stocks_By_Species.csv", header=T)

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
colnames(post.df.cwd) = c("tid","spid","treatment","species",
                          "posterior.mean","5","95","25","75")
k = 1
for (i in 1:n.t) {
  for (j in 1:n.spp.cwd) {
    a.mean = mean(a.cwd[,i,j], 2)
    a.hpdi.90 = HPDI(a.cwd[,i,j], prob=0.90)
    a.hpdi.50 = HPDI(a.cwd[,i,j], prob=0.50)
    post.df.cwd[k,"tid"] = i
    post.df.cwd[k,"spid"] = j
    post.df.cwd[k,"treatment"] = trt.names[i]
    post.df.cwd[k,"species"] = spp.cwd[j]
    post.df.cwd[k,"posterior.mean"] = a.mean
    post.df.cwd[k,"5"] = a.hpdi.90[1]
    post.df.cwd[k,"95"] = a.hpdi.90[2]
    post.df.cwd[k,"25"] = a.hpdi.50[1]
    post.df.cwd[k,"75"] = a.hpdi.50[2]
    k = k + 1
  }
}

## snag C stocks

# read in snag data
snag.melt = read.csv("Tree_Analysis/Clean_Data_By_Species/Snag_Carbon_Stocks_By_Species.csv", header=T)

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
colnames(post.df.snag) = c("tid","spid","treatment","species",
                           "posterior.mean","5","95","25","75")
k = 1
for (i in 1:n.t) {
  for (j in 1:n.spp.snag) {
    a.mean = mean(a.snag[,i,j], 2)
    a.hpdi.90 = HPDI(a.snag[,i,j], prob=0.90)
    a.hpdi.50 = HPDI(a.snag[,i,j], prob=0.50)
    post.df.snag[k,"tid"] = i
    post.df.snag[k,"spid"] = j
    post.df.snag[k,"treatment"] = trt.names[i]
    post.df.snag[k,"species"] = spp.snag[j]
    post.df.snag[k,"posterior.mean"] = a.mean
    post.df.snag[k,"5"] = a.hpdi.90[1]
    post.df.snag[k,"95"] = a.hpdi.90[2]
    post.df.snag[k,"25"] = a.hpdi.50[1]
    post.df.snag[k,"75"] = a.hpdi.50[2]
    k = k + 1
  }
}

## large woody biomass C stocks

# prep dataframe
treeC.dat = read.csv("Tree_Analysis/Clean_Data_By_Species/WoodyBiomass_C_Stocks_By_Species.csv")

# make full species name column
treeC.dat = unite(treeC.dat, col='spp', c('genus','species'), sep=' ')

# make wide df
treeC.wide = pivot_wider(treeC.dat, 
                         id_cols = c(trt, num),
                         names_from = spp, 
                         values_from = bmC.ha3)
treeC.wide[is.na(treeC.wide)] = 0

# reshape dataframe
treeC.melt = melt(treeC.wide, id.vars=c("treatment","plot"), variable.name="spp")

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
#traceplot(m.treeC)
#trankplot(m.treeC)

# get posterior distributions
post.treeC = extract.samples(m.treeC)
a.scale = exp(post.treeC$a) * mean.treeC
spp.treeC = levels(factor(sp.C.melt$spp))
n.spp.treeC = length(spp.treeC)

# create df with means and HPDIs by treatment and species
post.df.treeC = data.frame(matrix(nrow=n.t*n.spp.treeC, ncol=9))
colnames(post.df.treeC) = c("tid","spid","treatment","species",
                            "posterior.mean","5","95","25","75")
k = 1
for (i in 1:n.t) {
  for (j in 1:n.spp.bmC) {
    a.mean = mean(a.scale[,i,j], 2)
    a.sd = sd(a.scale[,i,j])
    a.hpdi.90 = HPDI(a.scale[,i,j], prob=0.90)
    a.hpdi.50 = HPDI(a.scale[,i,j], prob=0.50)
    post.df.treeC[k,"tid"] = i
    post.df.treeC[k,"spid"] = j
    post.df.treeC[k,"treatment"] = trt.names[i]
    post.df.treeC[k,"species"] = spp.treeC[j]
    post.df.treeC[k,"posterior.mean"] = a.mean
    post.df.treeC[k,"5"] = a.hpdi.90[1]
    post.df.treeC[k,"95"] = a.hpdi.90[2]
    post.df.treeC[k,"25"] = a.hpdi.50[1]
    post.df.treeC[k,"75"] = a.hpdi.50[2]
    k = k + 1
  }
}

## herbaceous biomass C stocks

# read in biomass C stock estimates
herbC.dat = read.csv("Tree_Analysis/Clean_Data_By_Plot/Clean_Veg_Soil_C_Stocks_Richness_by_Plot.csv")

# isolate biomass data for each species
herbC.vars = c("mixed.biomass.c.stock","p.arundinacea.c.stock","h.japonicus.c.stock")
herbC.spp = c("Mixed community","Phalaris arundinacea","Humulus japonicus")
herbC.spp.df = herbC.dat[,c("treatment","full.treatment.name","plot",herbC.vars)]
colnames(herbC.spp.df)[4:6] = herbC.spp

# melt dataframe
herbC.spp.melt = melt(herbC.spp.df, 
                      id.vars=c("treatment","full.treatment.name","plot"),
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
colnames(post.df.herbC) = c("tid","spid","treatment","species",
                            "posterior.mean","5","95","25","75")
k = 1
for (i in 1:n.t) {
  for (j in 1:n.spp.herbC) {
    a.mean = mean(a.herbC.sp[,i,j], 2)
    a.hpdi.90 = HPDI(a.herbC.sp[,i,j], prob=0.90)
    a.hpdi.50 = HPDI(a.herbC.sp[,i,j], prob=0.50)
    post.df.herbC[k,"tid"] = i
    post.df.herbC[k,"spid"] = j
    post.df.herbC[k,"treatment"] = trt.names[i]
    post.df.herbC[k,"species"] = spp.herbC[j]
    post.df.herbC[k,"posterior.mean"] = a.mean
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
post.df.all = rbind(post.df.snag, post.df.treeC, post.df.cwd, post.df.herbC)
write.csv(post.df.all,
          "Tree_Analysis/Posteriors/Vegetation_Carbon_Stocks_Species_Means_Intervals_5Chains.csv",
          row.names=F)


