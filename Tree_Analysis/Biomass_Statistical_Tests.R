setwd("C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/TreeData")

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

### statistical analyses on live and dead C stock components by plot ----

# read in biomass C stock estimates
bm.dat = read.csv("All_Biomass_Variables.csv")

# convert treatments into indices
bm.dat$tid = as.numeric(factor(bm.dat$trt))

# make simplified data lists
vars.short = c("fwd","cwd","snag","liveC","lstems","dstems","fwdC","fwdCN","litterC","litterCN","herbC","herbCN")
vars.long = c("fwd.count","cwd.area.sum","snag.area","live.woody.C.stock","live.stem.count","dead.stem.count",
              "FWD.C.g.m2","FWD.CN.ratio","HL.C.g.m2","HL.CN.ratio","total.C.g.m2","bm.CN.ratio")
n.v = length(vars.short)
all.df = list()
all.mean = list()
all.max = list()
all.sd = list()
for (i in 1:n.v) {
  v.short = vars.short[i]
  v.long = vars.long[i]
  all.df[[v.short]] = list(tid = bm.dat$tid, y = bm.dat[,v.long]/mean(bm.dat[,v.long]))
  all.mean[[v.short]] = mean(bm.dat[,v.long])
  all.sd[[v.short]] = sd(bm.dat[,v.long])
}

# look at histograms
par(mfrow=c(3,4))
for (i in 1:n.v) {
  v = vars.short[i]
  #hist(all.df[[v]]$y)
  hist(log(all.df[[v]]$y), main="")
}

par(mfrow=c(1,2))
v = vars.short[1]
hist(all.df[[v]]$y)
hist(log(all.df[[v]]$y))

# fit simple models
set.seed(12)
m.all = list()
for (i in 1:n.v) {
  v = vars.short[i]
  m.all[[v]] = ulam(alist(y ~ dnorm(mu, sigma),
                          log(mu) <- a[tid],
                          a[tid] ~ dnorm(0, 1),
                          sigma ~ dexp(1)),
                    data = all.df[[v]], chains=4)
}

# MCMC checks
m.all = m.all3
for (i in 1:n.v) { traceplot(m.all[[vars.short[i]]]) }
for (i in 1:n.v) { trankplot(m.all[[vars.short[i]]]) }

# get posterior distributions
post.all = list()
for (i in 1:n.v) {
  v = vars.short[[i]]
  post.all[[v]] = extract.samples(m.all[[v]])
}

# put posterior samples into list
vars.labels = c("Live woody C stock (Mg/ha)","Snag area (m2/ha)","CWD [>= 7.6 cm] area (m2/km)",
                "Live stems (#/ha)","Dead stems (#/ha)","FWD [2.5-7.5 cm] counts (#/m)",
                "Herbaceous C stock (Mg/ha)","Litter C stock (Mg/ha)","FWD [<2.5 cm] C stock (Mg/ha)",
                "Herbaceous C:N ratio","Litter C:N ratio","FWD C:N ratio")
vars = c("liveC","snag","cwd","fwd","lstems","dstems","herbC","litterC","fwdC","herbCN","litterCN","fwdCN")
n.v = length(vars)
all.a = list()
for (i in 1:n.v) {
  v = vars[i]
  all.a[[vars.labels[i]]] = exp(post.all[[v]]$a)*all.mean[[v]]
}

# melt posterior samples into dataframe
new.df = all.a[["Live woody C stock (Mg/ha)"]]
colnames(new.df) = trt.names
all.melt = melt(new.df, id.vars=trt.names)
all.melt$var = "Live woody C stock (Mg/ha)"
all.melt = all.melt[,2:4]
colnames(all.melt) = c("trt","value","var")
for (v in vars.labels[2:n.v]) {
  new.df = all.a[[v]]
  colnames(new.df) = trt.names
  melt.df = melt(new.df, id.vars=trt.names) 
  melt.df$var = v
  melt.df = melt.df[,2:4]
  colnames(melt.df) = c("trt","value","var")
  all.melt = rbind(all.melt, melt.df)
}
all.melt$var = factor(all.melt$var, levels=vars.labels)

# plot posterior distributions
ggplot(all.melt, aes(x=value, fill=trt)) + 
       geom_density(alpha=0.5, linewidth=0.2) + 
       facet_wrap(.~var, scales = "free") +
       labs(y="Density",x="") + 
       theme(legend.title=element_blank())

# create df with means and HPDIs by treatment and species
post.df.bm = data.frame(matrix(nrow=n.trt*n.v, ncol=8))
colnames(post.df.bm) = c("tid","trt","var","mean","5","95","25","75")
k = 1
for (j in 1:n.v) {
  for (i in 1:n.trt) {
    v = vars.labels[j]
    a.v = all.a[[v]]
    a.mean = mean(a.v[,i], 2)
    a.hpdi.90 = HPDI(a.v[,i], prob=0.90)
    a.hpdi.70 = HPDI(a.v[,i], prob=0.70)
    a.hpdi.50 = HPDI(a.v[,i], prob=0.50)
    post.df.bm[k,"tid"] = i
    post.df.bm[k,"trt"] = trt.names[i]
    post.df.bm[k,"var"] = vars.labels[j]
    post.df.bm[k,"mean"] = a.mean
    post.df.bm[k,"5"] = a.hpdi.90[1]
    post.df.bm[k,"95"] = a.hpdi.90[2]
    post.df.bm[k,"25"] = a.hpdi.50[1]
    post.df.bm[k,"75"] = a.hpdi.50[2]
    #post.df.bm[k,"15"] = a.hpdi.70[1]
    #post.df.bm[k,"85"] = a.hpdi.70[2]
    k = k + 1
  }
}

# plot posteriors
ggplot(post.df.bm, aes(y=factor(trt, levels=trt.names), x=mean)) + 
       geom_errorbar(aes(xmin=`5`, xmax=`95`), width=0.25, color="black") +
       geom_errorbar(aes(xmin=`25`, xmax=`75`), width=0.25, color="blue") +
       geom_point() + labs(x="", y="") +
       facet_wrap(.~factor(var, levels=vars.labels), scales="free_x", ncol=3) 

## convert posterior estimates back to natural scale
all.post.df = data.frame(matrix(nrow=6, ncol=6))
colnames(all.post.df) = trts

all.post.df[1,] = round(exp(apply(post.liveC$a, 2, mean))*sd(bm.dat$Live.woody.C.stock), 1)
all.post.df[2,] = round(exp(apply(post.snag$a, 2, mean))*sd(bm.dat$Snag.area), 1)
all.post.df[3,] = round(exp(apply(post.cwd$a, 2, mean))*sd(bm.dat$CWD.area), 1)
all.post.df[4,] = round(exp(apply(post.fwd$a, 2, mean))*sd(bm.dat$FWD.count), 1)
all.post.df[5,] = round(exp(apply(post.lstems$a, 2, mean))*sd(bm.dat$Live.stems), 1)
all.post.df[6,] = round(exp(apply(post.dstems$a, 2, mean))*sd(bm.dat$Dead.stems), 1)

post.fwd.HPDI = round(exp(apply(post.fwd$a, 2, HPDI, prob=0.90))*sd(bm.dat$FWD.count),1)
post.cwd.HPDI = round(exp(apply(post.cwd$a, 2, HPDI, prob=0.90))*sd(bm.dat$CWD.area),1)
post.snag.HPDI = round(exp(apply(post.snag$a, 2, HPDI, prob=0.90))*sd(bm.dat$Snag.area),1)
post.liveC.HPDI = round(exp(apply(post.liveC$a, 2, HPDI, prob=0.90))*sd(bm.dat$Live.woody.C.stock),1)
post.lstems.HPDI = round(exp(apply(post.lstems$a, 2, HPDI, prob=0.90))*sd(bm.dat$Live.stems),1)
post.dstems.HPDI = round(exp(apply(post.dstems$a, 2, HPDI, prob=0.90))*sd(bm.dat$Dead.stems),1)

## posterior contrasts for biomass variables v. reference
diff.all = data.frame(matrix(nrow=0, ncol=3))
colnames(diff.all) = c("trt","value","var")
diff.mean = data.frame(matrix(nrow=5*n.v, ncol=7))
colnames(diff.mean) = c("trt","var","mean","5","95","25","75")
k = 1
for (v in vars) {
  v.post = all_a[[v]]
  diff.post = data.frame(matrix(nrow=nrow(v.post), ncol=5))
  colnames(diff.post) = trt.names[1:5]
  for (i in 1:5) {
    diff.post[,trt.names[i]] = v.post[,i]-v.post[,6]
    diff.mean[k,"trt"] = trt.names[i]
    diff.mean[k,"var"] = v
    diff.mean[k,"mean"] = mean(v.post[,i]-v.post[,6])
    hpdi.90 = HPDI(v.post[,i]-v.post[,6], prob=0.90)
    hpdi.50 = HPDI(v.post[,i]-v.post[,6], prob=0.50)
    diff.mean[k,"5"] = hpdi.90[1]
    diff.mean[k,"95"] = hpdi.90[2]
    diff.mean[k,"25"] = hpdi.50[1]
    diff.mean[k,"75"] = hpdi.50[2]
    k = k + 1
  }
  diff.melt = melt(diff.post, variable.name="trt")
  diff.melt$var = v
  diff.all = rbind(diff.all, diff.melt)
}
  
ggplot(diff.all, aes(x=value, fill=trt)) + 
       geom_density(alpha=0.5, linewidth=0.2) + 
       labs(y="Density",x="") + 
       facet_wrap(.~factor(var, levels=vars), scales="free") + 
       theme(legend.title=element_blank())

ggplot(diff.mean, aes(y=factor(trt, levels=trt.names), x=mean)) + 
       geom_point() + labs(x="", y="") +
       geom_errorbar(aes(xmin=`5`, xmax=`95`), width=0.25, color="black") +
       geom_errorbar(aes(xmin=`25`, xmax=`75`), width=0.25, color="blue") +
       facet_wrap(.~factor(var, levels=vars), scales="free_x") +
       scale_fill_manual(values = c("5"="blue","25"="black"),
                         labels = c("50","90"))
