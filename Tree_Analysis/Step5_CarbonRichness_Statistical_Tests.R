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
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")
n.trt = length(trts)

################################################################################
## statistical analyses on live and dead C stock components by plot

# read in biomass C stock estimates
veg.soil.df = read.csv("Tree_Analysis/Clean_Data/Clean_Veg_Soil_Cstocks_Richness.csv")

# variable names
vars = colnames(veg.soil.df)[4:38]
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
set.seed(12)
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
var.labels = c("Fine woody debris (2.5-7.5 cm)","Coarse woody debris (\u2265 7.6 cm)",
               "Standing dead trees (\u2265 2.5 cm)","Live F. pennsylvanica (\u2265 2.5 cm)",
               "Standing dead F. pennsylvanica (\u2265 2.5 cm)","Diff. F. pennsylvanica (Live-Dead)",
               "Dead stems (< 2.5 cm)","Live stems (< 2.5 cm)","Live trees (\u2265 2.5 cm)",
               "Herbaceous biomass","Mixed biomass","P. arundinacea biomass",
               "H. japonicus biomass","Herbaceous litter","Fine woody debris",
               "Fine woody debris C:N ratio","Herbaceous litter C:N ratio","Herbaceous biomass C:N ratio",
               "MAOC","POC","SOC","TIC","TC",
               "Herbaceous species","Tree species","Total richness",
               "Total live vegetation","Total dead vegetation","Total vegetation","Total ecosystem",
               "Live trees w/ F. pennsylvanica","Standing dead trees w/o F. pennsylvanica",
               "Total live w/ F. pennsylvanica","Total dead w/o F. pennsylvanica",
               "Total vegetation w/o F. pennsylvanica death","Total ecosystem w/o F. pennsylvanica death")
length(vars)
length(var.labels)

# extract posterior samples and create dataframes
all.samples = list()
all.a = list()
all.b = 
a.df = data.frame(matrix(nrow=0, ncol=7))
colnames(a.df) = c("var",trts)
for (i in 1:n.v) {
  v1 = vars[i]
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
post.df = data.frame(matrix(nrow=n.v*n.trt, ncol=8))
colnames(post.df) = c("tid","trt","var","mean","5","95","25","75")
k = 1
for (j in 1:length(vars)) {
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
write.csv(post.df,"Tree_Analysis/Posteriors/CarbonRichness_Means_Intervals_5Chains.csv",
          row.names=F)

################################################################################
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

################################################################################
### statistical analysis on carbon stock data by species

## CWD C stocks

# read in CWD data
cwd.melt = read.csv("Clean_Data/CWD_Carbon_Stocks_By_Species.csv", header=T)

# convert treatment & species into indices
cwd.melt$tid = as.numeric(factor(cwd.melt$trt))
cwd.melt$spid = as.numeric(factor(cwd.melt$species))

# plot histogram
#par(mfrow=c(1,2))
#hist(log(cwd.melt$value))
#hist(log(cwd.melt$value/mean(cwd.melt$value)))

# make simplified data list
d.cwd.sp = list(tid = cwd.melt$tid, 
                spid = cwd.melt$spid,
                y = cwd.melt$cwd.carbon/mean(cwd.melt$cwd.carbon)) 
mean.cwd = mean(cwd.melt$cwd.carbon)

# fit simple model
set.seed(12)
m.cwd.sp = ulam(alist(y ~ dnorm(mu, sigma),
                      log(mu) <- a[tid, spid],
                      matrix[tid, spid]:a ~ dnorm(0, 1),
                      sigma ~ dexp(1)),
                data = d.cwd.sp, chains=5)

# check model
#traceplot(m.cwd.sp)
#trankplot(m.cwd.sp)

# get posterior distributions
post.cwd.sp = extract.samples(m.cwd.sp)
a.cwd = exp(post.cwd.sp$a)*mean.cwd
n.trt = 6
cwd.spp = levels(factor(cwd.melt$species))
n.spp.cwd = length(cwd.spp)

# create df with means and HPDIs by treatment and species
int.df.cwd = data.frame(matrix(nrow=n.trt*n.spp.cwd, ncol=9))
colnames(int.df.cwd) = c("tid","spid","trt","spp","mean","5","95","25","75")
k = 1
for (i in 1:n.trt) {
  for (j in 1:n.spp.cwd) {
    a.mean = mean(a.cwd[,i,j], 2)
    a.hpdi.90 = HPDI(a.cwd[,i,j], prob=0.90)
    a.hpdi.50 = HPDI(a.cwd[,i,j], prob=0.50)
    int.df.cwd[k,"tid"] = i
    int.df.cwd[k,"spid"] = j
    int.df.cwd[k,"trt"] = trt.names[i]
    int.df.cwd[k,"spp"] = cwd.spp[j]
    int.df.cwd[k,"mean"] = a.mean
    int.df.cwd[k,"5"] = a.hpdi.90[1]
    int.df.cwd[k,"95"] = a.hpdi.90[2]
    int.df.cwd[k,"25"] = a.hpdi.50[1]
    int.df.cwd[k,"75"] = a.hpdi.50[2]
    k = k + 1
  }
}

## snag C stocks

# read in snag data
snag.melt = read.csv("Clean_Data/Snag_Carbon_Stocks_By_Species.csv", header=T)

# convert treatment and species into id
snag.melt$tid = as.numeric(factor(snag.melt$trt))
snag.melt$spid = as.numeric(factor(snag.melt$species))

# make simplified data list
d.snag.sp = list(tid = snag.melt$tid, 
                 spid = snag.melt$spid,
                 y = snag.melt$snag.carbonmin/mean(snag.melt$snag.carbonmin)) 
snag.mean = mean(snag.melt$snag.carbonmin)

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
n.trt = 6
snag.spp = levels(factor(snag.melt$species))
n.spp = length(snag.spp)

# create df with means and HPDIs by treatment and species
post.df.snag = data.frame(matrix(nrow=n.trt*n.spp, ncol=9))
colnames(post.df.snag) = c("tid","spid","trt","spp","mean","5","95","25","75")
k = 1
for (i in 1:n.trt) {
  for (j in 1:n.spp) {
    a.mean = mean(a.snag[,i,j], 2)
    a.hpdi.90 = HPDI(a.snag[,i,j], prob=0.90)
    a.hpdi.50 = HPDI(a.snag[,i,j], prob=0.50)
    post.df.snag[k,"tid"] = i
    post.df.snag[k,"spid"] = j
    post.df.snag[k,"trt"] = trt.names[i]
    post.df.snag[k,"spp"] = snag.spp[j]
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
sp.C.dat = read.csv("Clean_Data/WoodyBiomass_C_Stocks_By_Species.csv")
sp.C.dat.2022 = sp.C.dat[which(sp.C.dat$year == 2022), 
                         c("trt","num","genus","species","abC_ha3")]

# make full species name column
sp.C.dat.2022 = unite(sp.C.dat.2022, col='spp', c('genus','species'), sep=' ')

# make wide df
sp.C.wide = pivot_wider(sp.C.dat.2022, 
                        id_cols = c(trt, num),
                        names_from = spp, 
                        values_from = abC_ha3)
sp.C.wide[is.na(sp.C.wide)] = 0

# reshape dataframe
sp.C.melt = melt(sp.C.wide, id.vars=c("trt","num"), variable.name="spp")

# convert treatments and species into indices
sp.C.melt$tid = as.numeric(factor(sp.C.melt$trt))
sp.C.melt$spid = as.numeric(factor(sp.C.melt$spp))

# make simplified data list
d.bmC = list(tid = sp.C.melt$tid, 
             spid = sp.C.melt$spid,
             y = sp.C.melt$value/mean(sp.C.melt$value))
mean.bmC = mean(sp.C.melt$value)

# fit simple model
set.seed(12)
m.bmC = ulam(alist(y ~ dnorm(mu, sigma),
                   log(mu) <- a[tid, spid],
                   matrix[tid, spid]:a ~ dnorm(0, 1),
                   sigma ~ dexp(1)),
             data = d.bmC, chains=5)

# model checks
#traceplot(m.bmC)
#trankplot(m.bmC)

# get posterior distributions
post.bmC = extract.samples(m.bmC)
a.scale = exp(post.bmC$a) * mean.bmC
all.spp = levels(factor(sp.C.melt$spp))
n.spp = length(all.spp)

# create df with means and HPDIs by treatment and species
post.df.woodC = data.frame(matrix(nrow=n.trt*n.spp, ncol=9))
colnames(post.df.woodC) = c("tid","spid","trt","spp","mean","5","95","25","75")
k = 1
for (i in 1:n.trt) {
  for (j in 1:n.spp) {
    a.mean = mean(a.scale[,i,j], 2)
    a.sd = sd(a.scale[,i,j])
    a.hpdi.90 = HPDI(a.scale[,i,j], prob=0.90)
    a.hpdi.50 = HPDI(a.scale[,i,j], prob=0.50)
    post.df.woodC[k,"tid"] = i
    post.df.woodC[k,"spid"] = j
    post.df.woodC[k,"trt"] = trt.names[i]
    post.df.woodC[k,"spp"] = all.spp[j]
    post.df.woodC[k,"mean"] = a.mean
    post.df.woodC[k,"5"] = a.hpdi.90[1]
    post.df.woodC[k,"95"] = a.hpdi.90[2]
    post.df.woodC[k,"25"] = a.hpdi.50[1]
    post.df.woodC[k,"75"] = a.hpdi.50[2]
    k = k + 1
  }
}

### statistical analyses on biomass C stocks by species ----

# read in biomass C stock estimates
bm.dat = read.csv("Clean_Data/Clean_Veg_Soil_Cstocks_Richness.csv")

# isolate biomass data for each species
bm.vars = c("PA.C.Mg.ha","HJ.C.Mg.ha","mix.C.Mg.ha")
bm.sp = c("Phalaris arundinacea","Humulus japonicus","Mixed community")
bm.sp.C.dat = bm.dat[,c("trt","trt.full","num",bm.vars)]
colnames(bm.sp.C.dat)[4:6] = bm.sp

# melt dataframe
bm.sp.C.melt = melt(bm.sp.C.dat, 
                    id.vars=c("trt","trt.full","num"),
                    variable.name="species")

# convert treatmentinto indices
bm.sp.C.melt$tid = as.numeric(factor(bm.sp.C.melt$trt))
bm.sp.C.melt$spid = as.numeric(factor(bm.sp.C.melt$species))

# make simple dataframe
d.bm = list(tid = as.integer(bm.sp.C.melt$tid),
            spid = as.integer(bm.sp.C.melt$spid),
            y = bm.sp.C.melt$value/mean(bm.sp.C.melt$value))

# fit model
set.seed(12)
m.bm.sp = ulam(alist(y ~ dnorm(mu, sigma),
                     log(mu) <- a[tid, spid],
                     matrix[tid, spid]:a ~ dnorm(0, 1),
                     sigma ~ dexp(1)),
               data = d.bm, chains=5)

# check model
#traceplot(m.bm.sp)
#trankplot(m.bm.sp)

# get posterior distributions
post.bm.sp = extract.samples(m.bm.sp)
a.bm.sp = exp(post.bm.sp$a) * mean(bm.sp.C.melt$value)
n.trt = 6
bm.sp = levels(factor(bm.sp.C.melt$species))
n.spp.bm = 3

# create df with means and HPDIs by treatment and species
post.df.bm.sp = data.frame(matrix(nrow=n.trt*n.spp.bm, ncol=9))
colnames(post.df.bm.sp) = c("tid","spid","trt","spp","mean","5","95","25","75")
k = 1
for (i in 1:n.trt) {
  for (j in 1:n.spp.bm) {
    a.mean = mean(a.bm.sp[,i,j], 2)
    a.hpdi.90 = HPDI(a.bm.sp[,i,j], prob=0.90)
    a.hpdi.50 = HPDI(a.bm.sp[,i,j], prob=0.50)
    post.df.bm.sp[k,"tid"] = i
    post.df.bm.sp[k,"spid"] = j
    post.df.bm.sp[k,"trt"] = trt.names[i]
    post.df.bm.sp[k,"spp"] = bm.sp[j]
    post.df.bm.sp[k,"mean"] = a.mean
    post.df.bm.sp[k,"5"] = a.hpdi.90[1]
    post.df.bm.sp[k,"95"] = a.hpdi.90[2]
    post.df.bm.sp[k,"25"] = a.hpdi.50[1]
    post.df.bm.sp[k,"75"] = a.hpdi.50[2]
    k = k + 1
  }
}

# plot posteriors
p.sp1 = ggplot(post.df.bm.sp, aes(y=factor(trt, levels=trt.names), 
                          fill=factor(spp, levels=bm.sp), x=mean)) + 
               geom_bar(stat = "identity") + 
               labs(y="", x="Posterior mean C stock (Mg/ha)", 
                    title = "",
                    fill="Species group") + 
               theme(text = element_text(size=14))
p.sp1
setwd(path_to_repo)
ggsave("Figures/Herbaceous_C_Stocks_By_Species.jpeg", 
       plot=p.sp1, width=18, height=10, units="cm")

## make combined color plot with live woody, snag/CWD area, and biomass stocks by species
post.df.woodC$type = "Live trees (\u2265 2.5 cm)"
post.df.snag$type = "Standing dead trees (\u2265 2.5 cm)"
int.df.cwd$type = "Coarse woody debris (\u2265 7.6 cm)"
types = c("Live trees (\u2265 2.5 cm)",
          "Standing dead trees (\u2265 2.5 cm)",
          "Coarse woody debris (\u2265 7.6 cm)")
post.df.sp.all = rbind(post.df.woodC, post.df.snag, int.df.cwd)
post.df.sp.all$spp = trimws(post.df.sp.all$spp)
all.sp = sort(unique(post.df.sp.all$spp))
p.sp2 = ggplot(post.df.sp.all, aes(x=mean, y=factor(trt, levels=trt.names),
                           fill=factor(spp, levels=all.sp))) + 
               geom_bar(stat="identity") +
               facet_wrap(.~factor(type, levels=types), scales="free_x", ncol=1) +
               labs(x="Posterior mean C stock (Mg/ha)", y="", fill="Species") +
               guides(fill=guide_legend(ncol=1)) + 
               theme(text = element_text(size=14))
ggsave("Figures/Woody_C_Stocks_By_Species.jpeg", 
       plot=p.sp2, width=17, height=19, units="cm")
