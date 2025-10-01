path_to_soil_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo/Soil_Analysis"
setwd(path_to_soil_folder)

library(rethinking)
library(dplyr)
library(reshape2)

# set seed 
set.seed(2718)

# treatment names
trt.letters = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

# read in soil data
soil.df = read.csv("Clean_Data/Soil_Data_by_Quadrat_June2023.csv", header=T)

# variables to omit
leave.out = c("ton.percent","pon.percent","texture.class","n.release",
              "volumetric.moisture","coarse.material",
              "k.sat","ca.sat","mg.sat","h.sat","h.meq",
              "elevation.mean","vegetation.height.mean")
soil.df = soil.df[,-which(colnames(soil.df) %in% leave.out)]

# scale moisture to percentage
soil.df$gravimetric.moisture = soil.df$gravimetric.moisture * 100

# make index for treatment, quadrat, and treatment-plot combination
soil.df$tid = as.numeric(factor(soil.df$treatment))
soil.df$qid = as.numeric(factor(soil.df$quadrat))
soil.df$tpid = rep(0, nrow(soil.df))
k = 1
for (t in 1:6) {
  for (p in 1:3) {
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
for (t in 1:6) {
  for (p in 1:3) {
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
for (i in n.v:n.v) {
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

# extract posterior samples
all.samples = list()
all.a = list()
for (i in n.v:n.v) {
  v = vars[i]
  all.samples[[v]] = extract.samples(soil.models[[v]])
  all.a[[v]] = exp(all.samples[[v]]$a) * all.means[[v]]
}

# variable labels
var.labels = list("temperature"="Temperature (C)",
                  "gravimetric.moisture"="Gravitational moisture (%)",
                  "bulk.density"="Bulk density (g/cm3)",
                  "bulk.n.percent"="TN (%)","bulk.c.percent"="TC (%)",
                  "bulk.cn.ratio"="C:N Ratio","toc.percent"="TOC (%)",
                  "sand"="Sand (%)","silt"="Silt (%)","clay"="Clay (%)",
                  "ph"="pH","p"="P (ppm)","k"="K (ppm)","ca"="Ca (ppm)",
                  "mg"="Mg (ppm)","som.percent"="SOM (%)",
                  "no3.n"="NO3 (ppm)","nh4.n"="NH4 (ppm)","cec"="CEC (meq/100 g)",
                  "k.meq"="K (meq)","ca.meq"="Ca (meq)","mg.meq"="Mg (meq)",
                  "poc.percent"="POC (%)","fpom"="Floating POM",
                  "g475mm"=">4.75 mm","g2mm"="2-4.75 mm","g250um"="0.250-2 mm",
                  "g53um"="0.053-0.250 mm","l53um"="<0.053 mm","mwd"="Mean-weight diameter (mm)",
                  "tic.percent"="TIC (%)","maoc.percent"="MAOC (%)",
                  "poc.maoc.ratio"="POC:MAOC ratio")

# get posterior intervals and write to file
df.int = data.frame(matrix(nrow=0, ncol=11))
colnames(df.int) = c("v","v.mean","t","mean","sd","5","95","15","85","25","75")
for (i in n.v:n.v) {
  df.v = data.frame(matrix(nrow=6, ncol=11))
  colnames(df.v) = c("v","v.mean","t","mean","sd","5","95","15","85","25","75")
  v = vars[i]
  a.v = all.a[[v]]
  a.means = apply(a.v, 2, mean)
  a.sds = apply(a.v, 2, sd)
  a.hpdi.90 = apply(a.v, 2, HPDI, prob=0.90)
  a.hpdi.70 = apply(a.v, 2, HPDI, prob=0.70)
  a.hpdi.50 = apply(a.v, 2, HPDI, prob=0.50)
  df.v[1:6,"v"] = var.labels[[v]]
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
