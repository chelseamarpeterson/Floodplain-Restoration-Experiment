path_to_soil_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo/Soil_Analysis"
setwd(path_to_soil_folder)

library(rethinking)
library(dplyr)
library(reshape2)

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
soil.df = read.csv("Clean_Data/Soil_Data_by_Quadrat_June2023.csv", header=T)

# variables to omit
leave.out = c("texture.class","volumetric.moisture","coarse.material","h.meq",
              "elevation.mean","vegetation.height.mean")
soil.df = soil.df[,-which(colnames(soil.df) %in% leave.out)]

# scale moisture to percentage
soil.df$gravimetric.moisture = soil.df$gravimetric.moisture * 100

# make index for treatment, quadrat, and treatment-plot combination
soil.df$tid = as.numeric(factor(soil.df$treatment))
soil.df$qid = as.numeric(factor(soil.df$quadrat))
soil.df$tpid = rep(0, nrow(soil.df))
k = 1
for (t in 1:n.t) {
  for (p in 1:n.p) {
    tp.ind = which(soil.df$tid == t & soil.df$plot == p)
    soil.df$tpid[tp.ind] = k
    k = k + 1
  }
}

# isolate fpom data
fpom.dat = soil.df[-which(is.na(soil.df$fpom)),
                   c("plot","treatment","full.treatment.name","plot","quadrat","fpom")]

# make index for treatment, quadrat, and plot number
fpom.dat$tid = as.numeric(factor(fpom.dat$treatment))
fpom.dat$qid = as.numeric(factor(fpom.dat$quadrat))
fpom.dat$tpid = rep(0, nrow(fpom.dat))
k = 1
for (t in 1:n.t) {
  for (p in 1:n.p) {
    tp.ind = which(fpom.dat$tid == t & fpom.dat$plot == p)
    fpom.dat$tpid[tp.ind] = k
    k = k + 1
  }
}

# make simple list for Bayes analysis
vars = colnames(soil.df)[5:(ncol(soil.df)-3)]
n.v = length(vars)
all.lists = list()
all.means = list()
for (i in 1:n.v) {
  v = vars[i]
  if (v == "fpom") {
    all.means[[v]] = mean(fpom.dat[,v])
    d.v =  list(tid=as.integer(fpom.dat$tid), 
                tpid=as.integer(fpom.dat$tpid), 
                y=fpom.dat[,v])
  } else {
    all.means[[v]] = mean(soil.df[,v])
    d.v =  list(tid=as.integer(soil.df$tid), 
                tpid=as.integer(soil.df$tpid), 
                y=soil.df[,v])
  }
  d.v$y = d.v$y/mean(d.v$y)
  all.lists[[v]] = d.v    
}

# model code w/ random effect
soil.models = list()
for (i in 1:n.v) {
  v = vars[i]
  m.v = ulam(alist(y ~ normal(mu, sigma),
                    log(mu) <- a[tid] + g[tid, tpid],
                    a[tid] ~ dnorm(0, 1),
                    matrix[tid, tpid]: g ~ dnorm(0, sigma_g),
                    sigma_g ~ dexp(1),
                    sigma ~ dexp(1)),
              data=all.lists[[v]], chains=5, log_lik=T)
  soil.models[[v]] = m.v
}

# MCMC checks
for (i in 1:n.v) { traceplot(soil.models[[vars[i]]]) }
for (i in 1:n.v) { trankplot(soil.models[[vars[i]]]) }

# extract posterior samples
all.samples = list()
all.a = list()
for (i in 1:n.v) {
  v = vars[i]
  all.samples[[v]] = extract.samples(soil.models[[v]])
  all.a[[v]] = exp(all.samples[[v]]$a) * all.means[[v]]
}

# variable labels
soil.var.labs = read.csv("Clean_Data/quadrat.variable.metadata.csv")
soil.var.list = list()
for (i in 1:n.v) {
  v1 = soil.var.labs[i,"variable"]
  v2 = soil.var.labs[i,"label"]
  soil.var.list[[v1]] = v2
}

# get posterior intervals and write to file
df.int = data.frame(matrix(nrow=0, ncol=10))
int.cols = c("variable","variable.label","variable.mean",
             "treatment","posterior.mean","posterior.sd",
             "5","95","25","75")
colnames(df.int) = 
for (i in 1:n.v) {
  df.v = data.frame(matrix(nrow=n.t, ncol=12))
  colnames(df.v) = int.cols
  v = vars[i]
  a.v = all.a[[v]]
  a.means = apply(a.v, 2, mean)
  a.sds = apply(a.v, 2, sd)
  a.hpdi.90 = apply(a.v, 2, HPDI, prob=0.90)
  a.hpdi.70 = apply(a.v, 2, HPDI, prob=0.70)
  a.hpdi.50 = apply(a.v, 2, HPDI, prob=0.50)
  df.v[1:n.t,"variable"] = v
  df.v[1:n.t,"variable.label"] = soil.var.list[[v]]
  df.v[1:n.t,"variable.mean"] = all.means[[v]]
  df.v[1:n.t,"treatment"] = trt.names
  df.v[1:n.t,"posterior.mean"] = a.means
  df.v[1:n.t,"posterior.sd"] = a.sds
  df.v[1:n.t,"5"] = a.hpdi.90[1,]
  df.v[1:n.t,"95"] = a.hpdi.90[2,]
  df.v[1:n.t,"25"] = a.hpdi.50[1,]
  df.v[1:n.t,"75"] = a.hpdi.50[2,]
  df.int = rbind(df.int, df.v)
}

write.csv(df.int, "Posteriors/All_Soil_Posterior_Intervals_ChainCount5.csv", row.names=F)
