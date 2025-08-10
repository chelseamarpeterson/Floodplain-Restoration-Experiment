path_to_soil_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Public-Repo/Soil_Analysis"
setwd(path_to_soil_folder)

library(rethinking)
library(dplyr)
library(reshape2)

# set seed 
set.seed(2718)

# treatment names
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

# read in soil data
soil.dat = read.csv("Clean_Data/All_Soil_Data_2023.csv", header=T)

# column names for groups of related variables
ag.cols = c("fpom","g475mm","g2mm","g250um","g53um","l53um","mwd")
meq.cols = c("k.meq","ca.meq","mg.meq","h.meq")
text.cols  = c("sand","silt","clay")

# variables to omit
leave.out = c("ton","pon")
soil.dat = soil.dat[,-which(colnames(soil.dat) %in% leave.out)]

# update column names
colnames(soil.dat)[6:9] = c("temp","gmoist","vmoist","bd") 
soil.dat$gmoist = soil.dat$gmoist * 100

# make index for treatment, quadrat, and plot number
soil.dat$tid = as.numeric(factor(soil.dat$trt))
soil.dat$qid = as.numeric(factor(soil.dat$quad))
soil.dat$pid = rep(0, nrow(soil.dat))
k = 1
for (t in 1:6) {
  for (p in 1:3) {
    tp.ind = which(soil.dat$tid == t & soil.dat$num == p)
    soil.dat$pid[tp.ind] = k
    k = k + 1
  }
}

# isolate fpom data
fpom.dat = soil.dat[-which(is.na(soil.dat$fpom)),c("plot","trt","trt.full","num","quad","fpom")]

# make index for treatment, quadrat, and plot number
fpom.dat$tid = as.numeric(factor(fpom.dat$trt))
fpom.dat$qid = as.numeric(factor(fpom.dat$quad))
fpom.dat$pid = rep(0, nrow(fpom.dat))
k = 1
for (t in 1:6) {
  for (p in 1:3) {
    tp.ind = which(fpom.dat$tid == t & fpom.dat$num == p)
    fpom.dat$pid[tp.ind] = k
    k = k + 1
  }
}

# add maoc:poc ratio
soil.dat$maoc_poc = soil.dat$maoc/soil.dat$poc

# make simple list for Bayes analysis
vars = c("som","bulk.n","bulk.c","bulk.cn","toc","poc","tic","maoc","temp",
         "gmoist","vmoist","bd","ph","p","k","ca","mg","no3","nh4","cec",
         text.cols,meq.cols[1:3],ag.cols,"maoc_poc") 
n.v = length(vars)
all.lists = list()
all.means = list()
for (i in 1:n.v) {
  v = vars[i]
  if (v == "fpom") {
    all.means[[v]] = mean(fpom.dat[,v])
    d.v =  list(tid=as.integer(fpom.dat$tid), 
                pid=as.integer(fpom.dat$num), 
                y=fpom.dat[,v])
  } else {
    all.means[[v]] = mean(soil.dat[,v])
    d.v =  list(tid=as.integer(soil.dat$tid), 
                pid=as.integer(soil.dat$num), 
                y=soil.dat[,v])
  }
  d.v$y = d.v$y/mean(d.v$y)
  all.lists[[v]] = d.v    
}

# model code w/ random effect
soil.models = list()
for (i in 1:n.v) {
  v = vars[i]
  m.v = ulam(alist(y ~ normal(mu, sigma),
                    log(mu) <- a[tid] + g[tid, pid],
                    a[tid] ~ dnorm(0, 1),
                    matrix[tid, pid]: g ~ dnorm(0, sigma_g),
                    sigma_g ~ dexp(1),
                    sigma ~ dexp(1)),
              data=all.lists[[v]], chains=5, log_lik=T)
  soil.models[[v]] = m.v
}

# extract posterior samples
all.samples = list()
all.a = list()
for (i in 1:n.v) {
  v = vars[i]
  all.samples[[v]] = extract.samples(soil.models[[v]])
  all.a[[v]] = exp(all.samples[[v]]$a) * all.means[[v]]
}

# variable labels
var.order = c("som","toc","tic","poc","maoc","bulk.c","bulk.n","bulk.cn",
              "temp","gmoist","bd",text.cols,"no3","nh4","p","k","ca","mg",
               meq.cols[1:3],"cec","ph",ag.cols,"maoc_poc")
var.labels = c("SOM (%)","TOC (%)","TIC (%)","POC (%)","MAOC (%)","TC (%)","TN (%)","C:N Ratio",
               "Temperature (C)","Gravitational moisture (%)","Bulk density (g/cm3)",
               "Sand (%)","Silt (%)","Clay (%)","NO3 (ppm)","NH4 (ppm)","P (ppm)",
               "K (ppm)","Ca (ppm)","Mg (ppm)","K (meq)","Ca (meq)","Mg (meq)",
               "CEC (meq/100 g)","pH","Floating POM",
               ">4.75 mm","2 - 4.75 mm","250 \u03bcm - 2 mm","53 - 250 \u03bcm",
               "<53 \u03bcm","Mean-weight diameter (mm)","MAOC:POC ratio")
n.v.plot = length(var.order)

# get posterior intervals and write to file
df.int = data.frame(matrix(nrow=0, ncol=11))
colnames(df.int) = c("v","v.mean","t","mean","sd","5","95","15","85","25","75")
for (i in 1:n.v.plot) {
  df.v = data.frame(matrix(nrow=6, ncol=11))
  colnames(df.v) = c("v","v.mean","t","mean","sd","5","95","15","85","25","75")
  v = var.order[i]
  a.v = all.a[[v]]
  a.means = apply(a.v, 2, mean)
  a.sds = apply(a.v, 2, sd)
  a.hpdi.90 = apply(a.v, 2, HPDI, prob=0.90)
  a.hpdi.70 = apply(a.v, 2, HPDI, prob=0.70)
  a.hpdi.50 = apply(a.v, 2, HPDI, prob=0.50)
  df.v[1:6,"v"] = var.labels[i]
  df.v[1:6,"v.mean"] = all.means[[v]]
  df.v[1:6,"t"] = trt.names
  df.v[1:6,"mean"] = a.means
  df.v[1:6,"sd"] = a.sds
  df.v[1:6,"5"] = a.hpdi.90[1,]
  df.v[1:6,"95"] = a.hpdi.90[2,]
  df.v[1:6,"15"] = a.hpdi.70[1,]
  df.v[1:6,"85"] = a.hpdi.70[2,]
  df.v[1:6,"25"] = a.hpdi.50[1,]
  df.v[1:6,"75"] = a.hpdi.50[2,]
  df.int = rbind(df.int, df.v)
}
write.csv(df.int, "Posteriors/All_Soil_Posterior_Intervals_ChainCount5.csv", row.names=F)
