setwd('C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/VegetationData')

library(ggplot2)
library(readxl)
library(ggfortify)
library(patchwork)
library(reshape2)
library(lme4)
library(dplyr)
library(tidyr)

## read in biomass data for 2022 and 2023
bm.2023 = read_excel("Joslin_Vegetation_Cover_Biomass_Aug2023.xlsx", sheet="Biomass")
bm.2022 = read_excel('HerbaceousBiomass_Sep2022.xlsx', sheet="ActualBiomass")

# remove boat mass columns
bm.2023 = bm.2023[, -which(colnames(bm.2023) %in% c("Mass of boat (g)","Mass of boat + vegetation (g)","Note"))]
bm.2022 = bm.2022[, -which(colnames(bm.2022) %in% c("Container mass (g)","Actual dry mass + container w/o bag (g)","Max diameter (cm)"))]

# update column names
colnames(bm.2023) = c("plot","quad","spp","mass")
colnames(bm.2022) = c("plot","quad","spp","mass")

# combine duplicate herbaceous litter and fine woody debris data within same quadrat
bm.2023 = bm.2023 %>% group_by(plot, quad, spp) %>% summarise(sum.mass = sum(mass))
bm.2022 = bm.2022 %>% group_by(plot, quad, spp) %>% summarise(sum.mass = sum(mass))

# add year column 
bm.2023$year = 2023
bm.2022$year = 2022

# combine dataframes
bm.all = rbind(bm.2022, bm.2023)

# normalize biomass data by quadrat size
bm.all$mass.g.m2 = bm.all$sum.mass/(30^2)*(100^2)

# separate plot column into treatment and plot number 
bm.all = bm.all %>% separate(plot, into = c("trt", "num"), sep = 1, remove=F)

# treatment vectors
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

# isolate biomass data
bm.dat = bm.all[!(bm.all$spp == "Herbaceous litter" | bm.all$spp == "Fine woody debris"),]
bm.dat$trt.full = rep(0, nrow(bm.dat))
for (i in 1:6) {bm.dat$trt.full[which(bm.dat$trt == trts[i])] = trt.names[i]}
bm.sp.sum = bm.dat %>% group_by(year, plot, trt, trt.full, spp) %>% summarise(mean.mass = sum(mass.g.m2)/3)

ggplot(bm.sp.sum, aes(y=factor(trt.full, levels=trt.names), x=mean.mass, fill=spp)) + 
       geom_bar(position="stack", stat="identity") +
       labs(y="Treatment", x="Herbaceous biomass (g/m2)", fill="Species") +
       facet_wrap(.~year)

# make data frame for total herbaceous biomass, herbaceous litter, and fine woody debris in each quadrat 
n.yr = 2; n.trt = 6; n.plt = 3; n.quad = 5
df.quad.sum = data.frame(matrix(nrow = n.yr*n.trt*n.plt*n.quad, ncol=11))
colnames(df.quad.sum) = c("year","trt","trt.full","num","quad","Herbaceous biomass","Herbaceous litter","Fine woody debris","Mixed species","Phalaris arundinacea","Humulus japonicus")
k = 1
for (y in c(2022, 2023)) {
  for (l in 1:n.trt) {
    t = trts[l]
    for (i in 1:n.plt) {
      for (j in 1:n.quad) {
        # get data for given year, treatment, plot, and quadrat
        quad.id = which((bm.all$year == y &bm.all$trt == t) & (bm.all$num == as.character(i) & bm.all$quad == paste("Q", j, sep="")))
        quad.df = bm.all[quad.id,]
        
        # get HL and FWD id
        hl.id = which(quad.df$spp == "Herbaceous litter")
        fwd.id = which(quad.df$spp == "Fine woody debris")
        pa.id = which(quad.df$spp == "Phalaris arundinacea")
        hj.id = which(quad.df$spp == "Humulus japonicus")
        bm.id = which(!(quad.df$spp %in% c("Herbaceous litter","Fine woody debris")))
        m.id = which(!(quad.df$spp %in% c("Herbaceous litter","Fine woody debris","Phalaris arundinacea","Humulus japonicus")))
                        
        # sum mass for each category
        df.quad.sum[k, "Herbaceous litter"] =  sum(quad.df[hl.id, "mass.g.m2"])
        df.quad.sum[k, "Fine woody debris"] = sum(quad.df[fwd.id, "mass.g.m2"])
        df.quad.sum[k, "Herbaceous biomass"] = sum(quad.df[bm.id, "mass.g.m2"])
        
        if (length(pa.id) > 0) { 
          df.quad.sum[k, "Phalaris arundinacea"] = sum(quad.df[pa.id, "mass.g.m2"]) 
        } else {
          df.quad.sum[k, "Phalaris arundinacea"] = 0
        }
        
        if (length(hj.id) > 0) { 
          df.quad.sum[k, "Humulus japonicus"] = sum(quad.df[hj.id, "mass.g.m2"]) 
        } else {
          df.quad.sum[k, "Humulus japonicus"] = 0
        }
        
        if (length(m.id) > 0) { 
          df.quad.sum[k, "Mixed species"] = sum(quad.df[m.id, "mass.g.m2"]) 
        } else {
          df.quad.sum[k, "Mixed species"] = 0
        }
        
        # fill in identification columns
        df.quad.sum[k, c("year","trt","trt.full","num","quad")] = c(y, t, trt.names[l], i, j)
        
        # increment counter
        k = k + 1
      }
    }
  }
}

df.quad.sum[-which(df.quad.sum[,"Herbaceous biomass"] == df.quad.sum[,"Mixed species"] + df.quad.sum[,"Phalaris arundinacea"] + df.quad.sum[,"Humulus japonicus"]),]

# melt biomass dataframe
df.melt = melt(df.quad.sum, id.vars = c("year","trt","trt.full","num","quad"), variable.name = "type", value.name="mass.g.m2")

# plot biomass components for each treatment by quadrat
ggplot(df.melt[which(df.melt$type %in% c("Herbaceous biomass","Herbaceous litter","Fine woody debris")),], aes(y=factor(trt.full, levels=trt.names), x=mass.g.m2, fill=factor(year))) + 
       geom_boxplot() + labs(y="Treatment", x="Total mass (g/m2)", fill="") + 
       theme(legend.title = element_blank(),legend.text = element_text(size=14), 
             axis.title = element_text(size=16), axis.text = element_text(size=12)) +
       facet_wrap(.~type, scales="free_x")

ggplot(df.melt[which(df.melt$type %in% c("Mixed species","Phalaris arundinacea","Humulus japonicus")),], aes(y=factor(trt.full, levels=trt.names), x=mass.g.m2, fill=factor(type))) + 
       geom_boxplot() + labs(y="Treatment", x="Total mass (g/m2)", fill="") + 
       theme(legend.title = element_blank(), legend.text = element_text(size=14), 
             axis.title = element_text(size=16), axis.text = element_text(size=12)) +
       facet_wrap(.~year, scales="free_x")

# nested ANOVAs comparing biomass components across treatments
aov.bm.nest <- aov(log(bm.mass.g.m2+1) ~ factor(trt) / factor(plot), data=df.quad.sum)
summary(aov.bm.nest)
plot(aov.bm.nest)

aov.hl.nest <- aov(log(hl.mass.g.m2+1) ~ factor(trt) / factor(plot), data=df.quad.sum)
summary(aov.hl.nest)
plot(aov.hl.nest)

aov.fwd.nest <- aov(log(fwd.mass.g.m2+1) ~ factor(trt) / factor(plot), data=df.quad.sum)
summary(aov.fwd.nest)

p1 = ggplot(df.quad.sum, aes(x=factor(trt_full,levels=trt_names), y=bm.mass.g.m2, fill=factor(plot))) + geom_boxplot(show.legend = FALSE) +
            labs(x = "", y = "Herbaceous\nbiomass (g/m2)") + theme(axis.text.x = element_blank())
p2 = ggplot(df.quad.sum, aes(x=factor(trt_full,levels=trt_names), y=hl.mass.g.m2, fill=factor(plot))) + geom_boxplot() + 
            labs(x = "", y = "Herbaceous\nlitter (g/m2)", fill='Plot') + theme(axis.text.x = element_blank()) 
p3 = ggplot(df.quad.sum, aes(x=factor(trt_full,levels=trt_names), y=fwd.mass.g.m2, fill=factor(plot))) + geom_boxplot(show.legend = FALSE) +
            labs(x = "", y = "Fine woody\ndebris (g/m2)")
p1 / p2 / p3

## linear mixed models comparing biomass components across treatments

# herbaceous biomass
df.quad.sum$log.bm = log(df.quad.sum$bm.mass.g.m2 + 1)
lme.bm = lmer(log.bm ~ trt + (1|trt:plot), data=df.quad.sum, REML=FALSE)
lme.bm.null = lmer(log.bm ~ (1|trt:plot), data=df.quad.sum, REML=FALSE)

summary(lme.bm)
plot(lme.bm)
hist(residuals(lme.bm))

library(pbkrtest)
drop1(lme.bm, test="Chisq")
KRmodcomp(lme.bm, lme.bm.null)

coef(lme.bm)
rand.eff = ranef(lme.bm)
fixed.eff = fixef(lme.bm)


# herbaceous litter
df.quad.sum$log.hl = log(df.quad.sum$hl.mass.g.m2 + 1)
lme.hl = lmer(log.hl ~ 1 + trt + (1|trt:plot), data=df.quad.sum, REML=FALSE)
lme.hl.null = lmer(log.hl ~ (1|trt:plot), data=df.quad.sum, REML=FALSE)
summary(lme.hl)
drop1(lme.hl, test="Chisq")
KRmodcomp(lme.hl, lme.hl.null)

# fine woody debris
df.quad.sum$log.fwd = log(df.quad.sum$fwd.mass.g.m2 + 1)
lme.fwd = lmer(log.fwd ~ 1 + trt + (1|trt:plot), data=df.quad.sum)
lme.fwd.null = lmer(log.fwd ~ (1|trt:plot), data=df.quad.sum, REML=FALSE)
summary(lme.fwd)
drop1(lme.fwd, test="Chisq")
KRmodcomp(lme.fwd, lme.fwd.null)


# dataframes for herb biomass, herb lit, and fwd in each plot summed over quadrats
bm.plt.total = data.frame(matrix(nrow=n.trt*n.plt, ncol=3))
hl.plt.total = data.frame(matrix(nrow=n.trt*n.plt, ncol=3))
fwd.plt.total = data.frame(matrix(nrow=n.trt*n.plt, ncol=3))
colnames(bm.plt.total) = c("trt","plt","mass.g.m2")
colnames(hl.plt.total) = c("trt","plt","mass.g.m2")
colnames(fwd.plt.total) = c("trt","plt","mass.g.m2")
k = 1
for (t in trts) {
  for (i in 1:n.plt) {
    bm.plt.id = which(bm.dat$plot == paste(t, i , sep=""))
    hl.plt.id = which(hl.dat$plot == paste(t, i , sep=""))
    fwd.plt.id = which(fwd.dat$plot == paste(t, i , sep=""))
    bm.plt.total[k, c("trt","plt")] = c(t, i)
    hl.plt.total[k, c("trt","plt")] = c(t, i)
    fwd.plt.total[k, c("trt","plt")] = c(t, i)
    bm.plt.total[k, "mass.g.m2"] = sum(bm.dat$dry.mass[bm.plt.id]/(30^2)*(100^2))/5
    hl.plt.total[k, "mass.g.m2"] = sum(hl.dat$dry.mass[hl.plt.id]/(30^2)*(100^2))/5
    fwd.plt.total[k, "mass.g.m2"] = sum(fwd.dat$dry.mass[fwd.plt.id]/(30^2)*(100^2))/5
    k = k + 1
  }
}

# bind biomass, litter, and fwd dfs
colnames(bm.plt.total)[3] = "Herbaceous biomass"
colnames(hl.plt.total)[3] = "Herbaceous litter"
colnames(fwd.plt.total)[3] = "Fine woody debris"
bm.melt = melt(bm.plt.total, id.vars = c("trt","plt"), variable.name = "Type", value.name="mass.g.m2")
hl.melt = melt(hl.plt.total, id.vars = c("trt","plt"), variable.name = "Type", value.name="mass.g.m2")
fwd.melt = melt(fwd.plt.total, id.vars = c("trt","plt"), variable.name = "Type", value.name="mass.g.m2")
fwd.melt[is.na(fwd.melt$mass.g.m2), "mass.g.m2"] = 0
df.total = rbind(bm.melt,hl.melt,fwd.melt)

ggplot(df.total, aes(x=trt, y=mass.g.m2, fill=factor(Type))) + geom_boxplot() +
       labs(x="Treatment", y="Total herbaceous biomass (g/m2)", fill="")

# sort biomass dataframe
bm.sort = bm.plt.sum[order(bm.plt.sum$plot,-bm.plt.sum$mass.g),]

# data frame of most abundant species in each plot
df.sp.sort = data.frame(matrix(nrow=n.trt*n.plt, ncol=max.n.sp+2))
colnames(df.sp.sort) = c("trt","plot",seq(1,max.n.sp))
k = 1
for (t in trts) {
  for (i in 1:n.plt) {
      plt.ind = which(bm.sort$trt == t & bm.sort$plot == i)
      plt.sp = bm.sort$spp[plt.ind]
      df.sp.sort[k,"trt"] = t
      df.sp.sort[k,"plot"] = i
      df.sp.sort[k,3:(length(plt.sp)+2)] = plt.sp
      k = k + 1
  }
}

# dataframe for number of unique species in each plot
df.uni.sp = data.frame(matrix(nrow=n.plots, ncol=5))
colnames(df.uni.sp) = c("trt","plot","n.sp","n.nat.sp","n.inv.sp")
native.sp = sp.2022$Spp[which(sp.2022$Native == 1)] 
k = 1
for (t in trts) {
  for (i in 1:n.plt) {
    plt.ind = which(bm.sort$trt == t & bm.sort$plot == i)
    plt.sp = bm.sort$spp[plt.ind]
    df.uni.sp[k,"trt"] = t
    df.uni.sp[k,"plot"] = i
    df.uni.sp[k,"n.sp"] = length(plt.sp)
    df.uni.sp[k,"n.nat.sp"] = sum(plt.sp %in% native.sp)
    df.uni.sp[k,"n.inv.sp"] = length(plt.sp) - sum(plt.sp %in% native.sp)
    k = k + 1
  }
}

# plot number of species for biomass samples
p1 = ggplot(df.uni.sp, aes(y=n.sp, x=trt)) + 
     geom_boxplot() +
     labs(x="Treatment", y="Total species", title="2022")
p2 = ggplot(df.uni.sp, aes(y=n.nat.sp, x=trt)) + 
     geom_boxplot() +
     labs(x="", y="# native species", title="2022")
p3 = ggplot(df.uni.sp, aes(y=n.inv.sp, x=trt)) + 
     geom_boxplot() +
     labs(x="Treatment", y="# non-native species", title="2022")
p1 + p2/p3

# sum up litter across quadrats by species
lit.plt.sum = data.frame(matrix(nrow=n.plots, ncol=2))
colnames(lit.plt.sum) = c("trt","mass.g")
k = 1
for (t in trts) {
  for (i in 1:n.plt) {
    plt.ind = which(lit.dat$plot == paste(t, i, sep=""))
    lit.plt.sum[k,"trt"] = t
    lit.plt.sum[k,"mass.g"] = sum(lit.dat$dry.mass[plt.ind])
    k = k + 1
  }
}
lit.plt.sum$mass.g.m2 = lit.plt.sum$mass.g/5/(30^2)*(100^2)

# litter biomass boxplots
ggplot(lit.plt.sum, aes(y=mass.g.m2, x=trt)) + 
  geom_boxplot() +
  labs(x="Treatment", y="Total mass (g)", title="2022")
  