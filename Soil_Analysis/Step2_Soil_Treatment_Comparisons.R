path_to_soil_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch2_Floodplain_Experiment/Floodplain-Experiment-Repo"
setwd(path_to_soil_folder)

library(rethinking)
library(dplyr)
library(reshape2)
library(brms)
library(bayestestR)

### script that reads quadrat-level soil data, runs statistical analyses on the
### variables of interest, and then writes posterior intervals to a file

# set seed 
set.seed(2718)

# treatment names
trt.letters = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")
n.t = length(trt.letters)
n.p = 3

# read in soil data
soil.df = read.csv("Soil_Analysis/Clean_Data/Soil_Data_by_Quadrat_June2023.csv", header=T)

# variables to omit
leave.out = c("texture.class","volumetric.moisture","coarse.material","h.meq",
              "elevation.mean","vegetation.height.mean","fpom")
soil.df = soil.df[,-which(colnames(soil.df) %in% leave.out)]

# scale moisture to percentage
soil.df$gravimetric.moisture = soil.df$gravimetric.moisture * 100

# make simple list for Bayes analysis
vars = colnames(soil.df)[6:ncol(soil.df)]
n.v = length(vars)
all.lists = list()
all.means = list()
for (i in 1:n.v) {
  v = vars[i]
  all.means[[v]] = mean(soil.df[,v])
  d.v =  list(treatment=as.factor(soil.df$treatment), 
              strip=as.factor(soil.df$strip), 
              plot=as.factor(soil.df$plot), 
              y=soil.df[,v])
  d.v$y = d.v$y/mean(d.v$y)
  all.lists[[v]] = d.v    
}

# run three models for each variable and store results in list
seed = 3141; n.iter=5000; n.chain=5
soil.model.list = list()
models = c("simple","plot.random","strip.plot.random")
model.labels = c("Fixed effects only","Plot random effects","Strip and plot random effects")
n.m = length(models)
for (i in 1:n.v) {
  v.i = vars[i]
  soil.model.list[[v.i]] = list()
  for (j in 1:n.m) {
    m.j = models[j]
    if (m.j == "simple") {
      model.fit.ij = brm(data = all.lists[[v.i]], 
                         y ~ treatment,
                         family=gaussian(link="log"),
                         chains=n.chain, seed=seed, iter=n.iter)
    } else if (m.j == "plot.random") {
      model.fit.ij = brm(data = all.lists[[v.i]], 
                         y ~ treatment + (1|treatment:plot),
                         family=gaussian(link="log"),
                         chains=n.chain, seed=seed, iter=n.iter)
    } else if (m.j == "strip.plot.random") {
      model.fit.ij = brm(data = all.lists[[v.i]], 
                         y ~ treatment + (1|treatment:strip) + (1|treatment:strip:plot),
                         family=gaussian(link="log"),
                         chains=n.chain, seed=seed, iter=n.iter)
    }
    model.fit.ij = add_criterion(model.fit.ij,"loo")
    model.fit.ij = add_criterion(model.fit.ij,"waic")
    soil.model.list[[v.i]][[m.j]] = model.fit.ij
  }
}

# save model fits
#for (i in 1:n.v) {
#  v.i = vars[i]
#  for (j in 1:n.m) {
#    m.j = models[j]
#    file.name = paste(v.i, m.j, "model.rds", sep="_")
#    saveRDS(soil.model.list[[v.i]][[m.j]], file = paste("Soil_Analysis/Posteriors/BRMS_Objects/",file.name,sep=""))
#  }
#}

# variable labels
soil.var.labs = read.csv("Soil_Analysis/Clean_Data/quadrat.variable.metadata.csv")
soil.var.list = list()
for (i in 1:nrow(soil.var.labs)) {
  v1 = soil.var.labs[i,"variable"]
  v2 = soil.var.labs[i,"label"]
  soil.var.list[[v1]] = v2
}

# save Rhat statistic for each variable and model
rhat.df = data.frame(matrix(nrow=0,ncol=7))
colnames(rhat.df) = c("variable","variable.label",
                      "model","model.label",
                      "treatment","full.treatment.name","Rhat")
for (i in 1:n.v) {
  v.i = vars[i]
  for (j in 1:n.m) {
    m.j = models[j]
    model.sum.ij = summary(soil.model.list[[v.i]][[m.j]])
    model.df.ij = data.frame(matrix(nrow=6,ncol=7))
    colnames(model.df.ij) = c("variable","variable.label",
                              "model","model.label",
                              "treatment","full.treatment.name","Rhat")
    model.df.ij$variable = v.i
    model.df.ij$variable.label = soil.var.list[[v.i]]
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
write.csv(rhat.df, "Soil_Analysis/Posteriors/Soil_Rhat_Statistic.csv", row.names=F)

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
    m1 = soil.model.list[[v.i]][[models[1]]]
    m2 = soil.model.list[[v.i]][[models[2]]]
    m3 = soil.model.list[[v.i]][[models[3]]]
    comp.ij = data.frame(loo_compare(m1, m2, m3, criterion = criteria[j],
                          model_names=models))
    colnames(comp.ij) = c("elpd_diff","se_diff","elpd","se_elpd","p","se_p","ic","se_ic")
    comp.ij$variable = v.i
    comp.ij$variable.label = soil.var.list[[v.i]]
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
write.csv(comp.df, "Soil_Analysis/Posteriors/Soil_Model_Information_Criteria.csv", row.names=F)

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
    fit.ij = soil.model.list[[v.i]][[m.j]]
    df.ij = data.frame(matrix(nrow=n.t, ncol=13))
    colnames(df.ij) = int.cols
    v = vars[i]
    hdi.90 = hdi(fit.ij, ci=0.90)
    hdi.50 = hdi(fit.ij, ci=0.50)
    df.ij$model = models[j]
    df.ij$model.label = model.labels[j]
    df.ij$variable = v.i
    df.ij$variable.label = soil.var.list[[v.i]]
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
write.csv(df.int, "Soil_Analysis/Posteriors/Soil_Posterior_Intervals_ChainCount5_BRMS.csv", row.names=F)
