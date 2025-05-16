path_to_soil_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Public-Repo/Soil_Data"
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

# labels
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
fpom.dat = soil.dat[-which(is.na(soil.dat$fpom)),
                    c("plot","trt","trt.full","num","quad","fpom")]

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

# make simple list for Bayes analysis
vars = c("som","bulk.n","bulk.c","bulk.cn","toc","poc","tic","maoc",
         "temp","gmoist","vmoist","bd",
         "ph","p","k","ca","mg","no3","nh4","cec",
         text.cols,meq.cols[1:3],ag.cols) 
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

# plot histograms of transformed variables
#par(mfrow=c(5,3), mar=c(2.5,3,1.5,1))
#for (i in 1:n.v) {
#  v = vars[i]
  #hist(log(all.lists1[[v]]$y/mean(all.lists1[[v]]$y)), main=v, xlab="")
#  hist(all.lists1[[v]]$y, main=v, xlab="")
#}

# model code w/o random effect
#set.seed(14)
#soil.models1 = list()
#soil.models1[[vars[1]]] = model1
#for (i in 2:n.v) {
#  v = vars[i]
#  m.v1 = ulam(alist(y ~ normal(mu, sigma),
#                    mu <- a[tid],
#                    a[tid] ~ dnorm(0, 1),
#                    sigma ~ dexp(1)), 
#              data=all.lists1[[v]], chains=5, log_lik=T)
#  soil.models1[[v]] = m.v1
#}

# model code w/ random effect
soil.models2 = list()
for (i in 1:n.v) {
  v = vars[i]
  m.v = ulam(alist(y ~ normal(mu, sigma),
                    log(mu) <- a[tid] + b[tid, pid],
                    a[tid] ~ dnorm(0, 1),
                    matrix[tid, pid]: b ~ dnorm(0, sigma_b),
                    sigma_b ~ dexp(1),
                    sigma ~ dexp(1)),
              data=all.lists[[v]], chains=5, log_lik=T)
  soil.models2[[v]] = m.v
}

# extract posterior samples
all.samples = list()
all.a = list()
all.sig = list()
for (i in 1:n.v) {
  v = vars[i]
  all.samples[[v]] = extract.samples(soil.models2[[v]])
  all.a[[v]] = exp(all.samples[[v]]$a) * all.means[[v]]
  all.sig[[v]] = exp(all.samples[[v]]$sigma) * all.means[[v]]
}

# variable labels
var.order = c("som","toc","tic","poc","maoc","bulk.c","bulk.n","bulk.cn",
              "temp","gmoist","bd",text.cols,"no3","nh4","p","k","ca","mg",
               meq.cols[1:3],"cec","ph",ag.cols)
var.labels = c("SOM (%)","TOC (%)","TIC","POM-C (%)","MAOM-C (%)","TC (%)","TN (%)","C:N Ratio",
               "Temperature (C)","Gravitational moisture (%)","Bulk density (g/cm3)",
               "Sand (%)","Silt (%)","Clay (%)","NO3 (ppm)","NH4 (ppm)","P (ppm)",
               "K (ppm)","Ca (ppm)","Mg (ppm)","K (meq)","Ca (meq)","Mg (meq)",
               "CEC (meq/100 g)","pH","Floating POM",
               ">4.75 mm","2 - 4.75 mm","250 \u03bcm - 2 mm","53 - 250 \u03bcm",
               "<53 \u03bcm","Mean-weight diameter (mm)")
n.v.plot = length(var.order)

# get posterior intervals and write to file
df.int = data.frame(matrix(nrow=0, ncol=10))
colnames(df.int) = c("v","t","mean","sd","5","95","15","85","25","75")
for (i in 1:n.v.plot) {
  df.v = data.frame(matrix(nrow=6, ncol=10))
  colnames(df.v) = c("v","t","mean","sd","5","95","15","85","25","75")
  v = var.order[[i]]
  a.means = apply(all.a[[v]], 2, mean)
  a.sds = apply(all.a[[v]], 2, sd)
  a.hpdi.90 = apply(all.a[[v]], 2, HPDI, prob=0.90)
  a.hpdi.70 = apply(all.a[[v]], 2, HPDI, prob=0.70)
  a.hpdi.50 = apply(all.a[[v]], 2, HPDI, prob=0.50)
  df.v[1:6,"v"] = var.labels[i]
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

#### compare WAIC and PSIS for all soil variables
#df.waic = data.frame(matrix(nrow=0, ncol=8))
#df.psis = data.frame(matrix(nrow=0, ncol=8))
#colnames(df.waic)[1:7] = c(colnames(compare(soil.models2[[vars[1]]], soil.models1[[vars[1]]], func="WAIC")),"var")
#colnames(df.psis)[1:7] = c(colnames(compare(soil.models2[[vars[1]]], soil.models1[[vars[1]]], func="PSIS")),"var")
#colnames(df.waic)[8] = "model"
#colnames(df.psis)[8] = "model"
#for (i in 1:n.v) {
#  v = var.new[i]
#  m1 = soil.models1[[v]]
#  m2 = soil.models2[[v]]
#  v.waic = compare(m1, m2, func="WAIC")
#  v.psis = compare(m1, m2, func="PSIS")
#  v.waic$var = var.labels[i]
#  v.psis$var = var.labels[i]
#  v.waic$model[which(row.names(v.waic) == "m1")] = "No random effects"
#  v.waic$model[which(row.names(v.waic) == "m2")] = "Random effects"
#  v.psis$model[which(row.names(v.psis) == "m1")] = "No random effects"
#  v.psis$model[which(row.names(v.psis) == "m2")] = "Random effects"
#  df.waic = rbind(df.waic, v.waic)
#  df.psis = rbind(df.psis, v.psis)
#}
#df.waic$metric = "Widely Applicable Information Criterion (WAIC)"
#df.psis$metric = "Pareto-smoothed Importance Sampling Cross-validation (PSIS))"
#df.all.met = rbind(df.waic, df.psis)

#mod.labs = c("No random effects","Random effects")
#df.waic$model = factor(df.waic$model, levels=mod.labs)
#p1 = ggplot(df.waic[which(df.waic$model=="Random effects"),], aes(y=factor(var, levels=var.labels), x=1-weight)) + 
#      geom_point() + xlim(c(0,1)) +
#      theme(axis.title.y=element_blank(), legend.position="none") + 
#      labs(x = "1 - (Random effects weight)") + scale_y_discrete(limits = rev)
#p2 = ggplot(df.waic, aes(y=factor(var, levels=var.labels), x=WAIC, color=model)) + 
#            geom_point(position=position_dodge(width=0.5)) + 
#            geom_errorbar(aes(y=var, xmin=WAIC-SE, xmax=WAIC+SE), width=.1, position=position_dodge(width=0.5)) +
#            theme(axis.title.y=element_blank(), legend.title=element_blank(), axis.text.y=element_blank()) + 
#            labs(x = "WAIC") + scale_y_discrete(limits = rev)
#p1 + p2
#p1 = ggplot(df.psis, aes(y=factor(var, levels=var.labels), x=PSIS, color=factor(model, levels=mod.labs))) + 
#            geom_point(position=position_dodge(width=0.5)) + 
#            geom_errorbar(aes(y=var, xmin=PSIS-SE, xmax=PSIS+SE), width=.1, position=position_dodge(width=0.5)) +
#            theme(axis.title.y=element_blank(), legend.title=element_blank(), axis.text.y=element_blank()) + 
#            labs(x = "PSIS") + scale_y_discrete(limits = rev)
#p2 = ggplot(df.waic[which(df.waic$model=="Random effects"),], aes(y=factor(var, levels=var.labels), x=1-weight)) + 
#            geom_point() + xlim(c(0,1)) +
#            theme(axis.title.y=element_blank(), legend.position="none") + 
#            labs(x = "1- (Random effects weight)") + scale_y_discrete(limits = rev)
#p2 + p1

## test simple hierarchical models
#model3.code = "
#data{
#     vector[90] y;
#     array[90] int pid;
#     array[90] int tid;
#}
#parameters{
#     matrix[6,3] b;
#     vector[6] a;
#     real<lower=0> sigma_t;
#     real<lower=0> sigma;
#}
#model{
#     vector[90] mu;
#     matrix[6,3] mu_t;
#     sigma_t ~ exponential( 1 );
#     sigma ~ exponential( 1 );
#     a ~ normal( 0 , 1 );
#     for ( i in 1:6 ) {
#        for (j in 1:3) {
#           mu_t[i, j] = a[i];
#        }
#     }
#     to_vector( b ) ~ normal( to_vector(mu_t) , sigma_t);
#     for ( i in 1:90 ) {
#         mu[i] = b[tid[i], pid[i]];
#     }
#     y ~ normal( mu , sigma );
#}"
#
#stan_test <- stan( model_code=model3.code , data=all.lists[["som"]] , chains=1)
#precis(stan_test, depth=3)
#compare( stan_test , model1 )


