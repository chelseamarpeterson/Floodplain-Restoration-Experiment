library(ggplot2)
library(readxl)
setwd('C:/Users/Chels/OneDrive - University of Illinois - Urbana/VegetationData')

library(ggfortify)
library(lme4)
library(comprehenr)
library(patchwork)
library(tidyverse)
library(reshape2)
library(lemon)
library(WorldFlora)

## script that analyzes vegetation cover data

# load data

# read in cover data
all.cov = read.csv("PlantCoverClean_2022.csv")

# important site variables
yrs = c(2013,2022); n.y = length(yrs)
trts = c("A","B","C","D","E","R"); n.t = length(trts); n.p = 3; n.q = 5
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")
forms = c("forb","vine","tree","graminoid"); n.f = length(forms)
form.vars = to_vec(for(i in 1:n.f) paste("n",forms[i],sep="."))

################################################################################
# plot vegetation summary data

# make dataframe for values by treatment, plot, quad
my.cols = c("year","trt","plot","quad",
            "n.sp","n.fam","n.ord",
            "n.nat.sp","n.inv.sp","n.per",
            "rcg.cov","hum.cov","inv.cov","nat.cov",
            "meanC","FQI","PI","mean.wis",
             form.vars)
df.veg = data.frame(matrix(nrow = n.y*n.t*n.p*n.q, ncol = length(my.cols)))
colnames(df.veg) = my.cols
form.ind = which(colnames(df.veg) %in% form.vars)
k = 1
for (y in yrs) {
  for (t in trts) {
    for (p in 1:n.p) {
      for (q in 1:n.q) {
        id = which((all.cov$year == y & all.cov$trt == t) & (all.cov$plot == p & all.cov$quad == q))
        df.veg[k, c("year","trt","plot","quad")] = c(y, t, p, q)
        
        # species richness
        n.sp = length(id)
        df.veg[k,"n.sp"] = n.sp

        # family richness
        n.fam = length(unique(all.cov[id,"family"]))
        df.veg[k,"n.fam"] = n.fam

        # order richness
        n.ord = length(unique(all.cov[id,"order"]))
        df.veg[k,"n.ord"] = n.ord
        
        # number of perennials
        df.veg[k,"n.per"] = sum(all.cov[id,"life.dur"])

        # native & invasive richness
        df.veg[k,"n.nat.sp"] = sum(all.cov[id,"native"])
        df.veg[k,"n.inv.sp"] = n.sp - sum(all.cov[id,"native"])
                
        # mean C, prevalence index, mean wis, physiognomy
        if (length(id) > 0) {
          df.veg[k,"meanC"] = sum(all.cov[id,"c.value"])/n.sp
          df.veg[k,"FQI"] = df.veg[k,"meanC"]*sqrt(n.sp)
          df.veg[k,"PI"] = sum(all.cov[id,"wis.num"]*all.cov[id,"cover"])/sum(all.cov[id,"cover"]) 
          df.veg[k,"mean.wis"] = sum(all.cov[id,"wis.num"])/n.sp
          
          phys = all.cov[id,"growth.form"]
          phys.vec = rep(0,n.sp) 
          for (i in 1:n.sp) {
            phys.split = strsplit(phys[i],"[ -]+")[[1]]
            phys.vec[i] = tolower(phys.split[length(phys.split)])
          }
          df.veg[k,"n.forb"] = sum(phys.vec=="forb")
          df.veg[k,"n.tree"] = sum(phys.vec %in% c("tree","shrub"))
          df.veg[k,"n.vine"] = sum(phys.vec=="vine")
          df.veg[k,"n.graminoid"] = sum(phys.vec %in% c("grass","sedge"))
          
        } else {
          df.veg[k,"meanC"] = 0
          df.veg[k,"FQI"] = 0
          df.veg[k,"PI"] = 0
          df.veg[k,"mean.wis"] = 0
          for (i in 1:n.f) {df.veg[k,form.ind[i]] = 0}
        }
        
        # rcg cover
        rcg.id = which(all.cov[id,"spp"] == "Phalaris arundinacea")
        if (length(rcg.id) > 0) {
          df.veg[k,"rcg.cov"] = all.cov[id[rcg.id],"cover"]          
        } else {
          df.veg[k,"rcg.cov"] = 0
        }

        # humulus cover
        hum.id = which(all.cov[id,"spp"] == "Humulus japonicus")
        if (length(hum.id) > 0) {
          df.veg[k,"hum.cov"] = all.cov[id[hum.id],"cover"]          
        } else {
          df.veg[k,"hum.cov"] = 0
        }
        
        # invasive cover
        inv.id = which(all.cov[id,"native"] == 0)
        if (length(inv.id) > 0) {
          df.veg[k,"inv.cov"] = sum(all.cov[id[inv.id],"cover"])         
        } else {
          df.veg[k,"inv.cov"] = 0
        }
        
        # native cover
        nat.id = which(all.cov[id,"native"] == 1)
        if (length(nat.id) > 0) {
          df.veg[k,"nat.cov"] = sum(all.cov[id[nat.id],"cover"])          
        } else {
          df.veg[k,"nat.cov"] = 0
        }
        k = k + 1
      }
    }
  }
}
df.melt = melt(df.veg, id.vars = c("year","trt","plot","quad"), variable.name = "var", value.name="value")
hum.label = expression(paste(italic("H. japonicus "), "percent cover", sep=" "))
rcg.label = expression(paste(italic("P. arundinacea "), "percent cover", sep=" "))

# plot vegetation survey results for 2022
df.melt.2022 = df.melt[which(df.melt$year==2022),]
p1 = ggplot(df.melt.2022[which(df.melt.2022$var %in% c("rcg.cov")),], aes(x=trt, y=value)) + 
            geom_boxplot() + labs(x="Treatment",y=rcg.label, fill="") + 
            scale_x_discrete(labels=trt.names) + coord_flip() +
            theme(legend.title = element_blank(),legend.text = element_text(size=14), 
                  axis.title = element_text(size=16), axis.text = element_text(size=12))
p2 = ggplot(df.melt.2022[which(df.melt.2022$var %in% c("hum.cov")),], aes(x=trt, y=value)) + 
            geom_boxplot() + labs(x = "", y = hum.label, fill="") + 
            scale_x_discrete(labels=trt.names) + coord_flip() +
            theme(legend.title = element_blank(),legend.text = element_text(size=14), 
                  axis.title = element_text(size=16), axis.text.y = element_blank(), axis.text = element_text(size=12))
p3 = ggplot(df.melt.2022[which(df.melt.2022$var %in% c("inv.cov")),], aes(x=trt, y=value)) + 
            geom_boxplot() + labs(x="",y="Invasive cover", fill="") + 
            scale_x_discrete(labels=trt.names) + coord_flip() +
            theme(legend.title = element_blank(),legend.text = element_text(size=14), 
                  axis.title = element_text(size=16), axis.text.y = element_blank(), axis.text = element_text(size=12))
p4 = ggplot(df.melt.2022[which(df.melt.2022$var %in% c("nat.cov")),], aes(x=trt, y=value)) + 
            geom_boxplot() + labs(x="Treatment",y="Native cover", fill="") + 
            scale_x_discrete(labels=trt.names) + coord_flip() +
            theme(legend.title = element_blank(),legend.text = element_text(size=14), 
                  axis.title = element_text(size=16), axis.text = element_text(size=12))
p5 = ggplot(df.melt.2022[which(df.melt.2022$var %in% c("meanC")),], aes(x=trt, y=value)) + 
            geom_boxplot() + labs(x="",y="Mean C", fill="") + 
            scale_x_discrete(labels=trt.names) + coord_flip() +
            theme(legend.title = element_blank(), legend.text = element_text(size=14), 
                  axis.title = element_text(size=16), axis.text.y = element_blank(), axis.text = element_text(size=12))
p6 = ggplot(df.melt.2022[which(df.melt.2022$var %in% c("n.sp")),], aes(x=trt, y=value)) + 
            geom_boxplot() + labs(x="",y="Richness", fill="") + 
            scale_x_discrete(labels=trt.names) + coord_flip() +
            theme(legend.title = element_blank(),legend.text = element_text(size=14), 
                  axis.title = element_text(size=16), axis.text.y = element_blank(), axis.text = element_text(size=12))
(p1+p2+p3)/(p4+p5+p6)

# plot 2013 & 2022 data
plot.vars = c("rcg.cov","hum.cov","inv.cov","nat.cov",
              "n.nat.sp","n.sp","n.per","meanC",
              "FQI","PI")
plot.labels=c("Phalaris cover","Humulus cover","Invasive cover","Native cover",
              "Native species richness","Total species richness","Perennial count","Mean C",
              "FQI","Prevalence Index")
df.plot = df.melt[which(df.melt$var %in% plot.vars),]
ggplot(df.plot, aes(y=value, x=trt, fill=factor(year))) + geom_boxplot() + 
       facet_wrap(~factor(var, levels=plot.vars, labels=plot.labels), scales="free_x", ncol=4) + coord_flip() + 
       labs(x="",y="", fill="") + scale_x_discrete(labels=trt.names)

# plot different taxonomic levels by year
tax.levs = c("Species richness","Family richness","Order richness")
tax.vars = c("n.sp","n.fam","n.ord")
df.tax = df.melt[which(df.melt$var %in% tax.vars),]
ggplot(df.tax, aes(x=trt, y=value, fill=factor(var))) + 
       geom_boxplot() + facet_grid(rows=vars(year)) +
       labs(y="Count",x="",fill="") + scale_x_discrete(labels=trt.names) +
       scale_color_discrete(breaks=tax.vars,labels=tax.levs)

# plot different functional types by year
df.form = df.melt[which(df.melt$var %in% form.vars),]
ggplot(df.form, aes(x=trt, y=value, fill=factor(var))) + 
       geom_boxplot() + facet_grid(rows=vars(year)) +
       labs(y="Count",x="",fill="") + scale_x_discrete(labels=trt.names)

# plot differences between years
df.form.diff = df.form[,-which(colnames(df.form)=="year")]
df.form.diff[,"value"] = df.form[which(df.form$year==2022),"value"] - df.form[which(df.form$year==2013),"value"]
ggplot(df.form.diff, aes(x=trt, y=value, fill=factor(var))) + geom_boxplot()

################################################################################
# statistical analyses on vegetation summary data

## linear mixed models on 2022 cover data
df.veg$log.rcg = log(1+df.veg$rcg.cov)
lme.rcg = lmer(log.rcg ~ trt + (1|trt:plot), data=df.veg[which(df.veg$year==2022),])
#summary(lme.rcg)
plot(lme.rcg)
hist(residuals(lme.rcg))
drop1(lme.rcg, test="Chisq")

df.veg$log.hum = log(1+df.veg$hum.cov)
lme.hum1 = lmer(hum.cov ~ trt + (1|trt:plot), data=df.veg[which(df.veg$year==2022),])
lme.hum2 = lmer(log.hum ~ trt + (1|trt:plot), data=df.veg[which(df.veg$year==2022),])
lme.hum2.null = lmer(log.hum ~ (1|trt:plot), data=df.veg[which(df.veg$year==2022),])
#summary(lme.hum)
plot(lme.hum1)
plot(lme.hum2)
hist(residuals(lme.hum1))
hist(residuals(lme.hum2))
hist(df.veg$hum.cov)
drop1(lme.hum2, test="Chisq")

df.veg$log.meanC = log(1+df.veg$meanC)
lme.meanC1 = lmer(meanC ~ trt + (1|trt:plot), data=df.veg[which(df.veg$year==2022),])
lme.meanC2 = lmer(log.meanC ~ trt + (1|trt:plot), data=df.veg[which(df.veg$year==2022),])
plot(lme.meanC1)
plot(lme.meanC2)
hist(residuals(lme.meanC1))
hist(residuals(lme.meanC2))
drop1(lme.meanC2, test="Chisq")

df.veg$log.richness = log(1+df.veg$richness)
lme.richness1 = lmer(richness ~ trt + (1|trt:plot), data=df.veg[which(df.veg$year==2022),])
lme.richness2 = lmer(log.richness ~ trt + (1|trt:plot), data=df.veg[which(df.veg$year==2022),])
plot(lme.richness1)
plot(lme.richness2)
hist(residuals(lme.richness1))
hist(residuals(lme.richness2))
drop1(lme.richness2, test="Chisq")

df.veg$log.nat = log(1+df.veg$nat.cov)
lme.nat = lmer(log.nat ~ trt + (1|trt:plot), data=df.veg[which(df.veg$year==2022),])
plot(lme.nat)
hist(residuals(lme.nat))
drop1(lme.nat, test="Chisq")

df.veg$log.inv = log(1+df.veg$inv.cov)
lme.inv = lmer(log.inv ~ trt + (1|trt:plot), data=df.veg[which(df.veg$year==2022),])
plot(lme.inv)
hist(residuals(lme.inv))
drop1(lme.inv, test="Chisq")

library(pbkrtest)
KRmodcomp(lme.hum2, lme.hum2.null)

## calculate change in variables from 2013 to 2022
df.veg.del = data.frame(matrix(nrow = n.t*n.p*n.q, ncol = 9))
colnames(df.veg.del) = colnames(df.veg)[2:10]
df.veg.del[,1:3] = df.veg[which(df.veg$year == 2013),2:4]
df.veg.del[,4:9] = df.veg[which(df.veg$year == 2022),5:10] - df.veg[which(df.veg$year == 2013),5:10]
df.del.melt = melt(df.veg.del, id.vars = c("trt","plot","quad"), variable.name = "var", value.name="cover")

# plot change in phalaris and humulus cover
ggplot(df.del.melt[which(df.del.melt$var %in% c("rcg.cov","hum.cov")),], aes(x=trt, y=cover, fill=var)) + 
       geom_boxplot() + labs(x="Treatment",y="Change in cover from 2013-2022", fill="") + 
       scale_fill_discrete(labels=c("Phalaris arundinacea","Humulus japonicus")) + 
       scale_x_discrete(labels=trt.names) + 
       theme(legend.title = element_blank(),legend.text = element_text(size=14), 
             axis.title = element_text(size=16), axis.text = element_text(size=12))

# plot change in native and invasive cover
ggplot(df.del.melt[which(df.del.melt$var %in% c("inv.cov","nat.cov")),], aes(x=trt, y=cover, fill=var))
       geom_boxplot()

# change in richness
p1 = ggplot(df.veg, aes(x=trt, y=richness)) + geom_boxplot(show.legend = FALSE) +
            labs(x = "", y = "Change in richness") + theme(axis.text.x = element_blank())
p2 = ggplot(df.veg, aes(x=trt, y=meanC)) + geom_boxplot(show.legend = FALSE) +
            labs(x = "Treatment", y = "Change in Mean C") + scale_x_discrete(labels=trt.names)
p1 / p2

# linear mixed models on change in cover data
df.veg.del$log.rcg = sign(df.veg.del$rcg.cov)*log(1+abs(df.veg.del$rcg.cov))
lme.rcg = lmer(log.rcg ~ trt + (1|trt:plot), data=df.veg.del)
summary(lme.rcg)
plot(lme.rcg)
hist(residuals(lme.rcg))
drop1(lme.rcg, test="Chisq")

df.veg.del$log.hum = sign(df.veg.del$hum.cov)*log(1+abs(df.veg.del$hum.cov))
lme.hum = lmer(log.hum ~ trt + (1|trt:plot), data=df.veg.del)
summary(lme.hum)
plot(lme.hum)
hist(residuals(lme.hum))
drop1(lme.hum, test="Chisq")

#df.veg.del$log.meanC = sign(df.veg.del$hum.cov)*log(1+abs(df.veg.del$hum.cov))
lme.meanC = lmer(meanC ~ trt + (1|trt:plot), data=df.veg.del)
summary(lme.meanC)
plot(lme.meanC)
hist(residuals(lme.meanC))
drop1(lme.meanC, test="Chisq")

lme.richness = lmer(richness ~ trt + (1|trt:plot), data=df.veg.del)
summary(lme.richness)
plot(lme.richness)
hist(residuals(lme.richness))
drop1(lme.richness, test="Chisq")

################################################################################ 
## calculate average changes across plots

# data frame for average changes
df.ave.del = data.frame(matrix(nrow=n.t*n.p, ncol=length(colnames(df.veg.del))-1)) # all columns except quad
colnames(df.ave.del) = colnames(df.veg.del)[-which(colnames(df.veg.del) == "quad")]
num.cols = colnames(df.veg.del)[-which(colnames(df.veg.del) %in% c("trt","plot","quad"))]
k = 1
for (t in trts) {
  for (i in 1:n.p) {
    p.ind = which(df.veg.del$trt == t & df.veg.del$plot == i)
    df.ave.del[k, c("trt","plot")] = c(t, i)
    df.ave.del[k, num.cols] = colMeans(df.veg.del[p.ind, num.cols])
    k = k + 1
  }
}

ggplot(df.ave.del, aes(x=hum.cov, y=rcg.cov, color=factor(trt))) + geom_point()
ggplot(df.ave.del, aes(x=hum.cov, y=richness, color=trt)) + geom_point()
ggplot(df.ave.del, aes(x=hum.cov, y=meanC, color=trt)) + geom_point()
## compare species richness & C-values across plots

# data frame for richness, mean C, and native/nonnative cover by plot
df.plot = data.frame(matrix(nrow=n.y*n.t*n.p, ncol=10))
colnames(df.plot) = c("year","trt","trt.full","plot","richness","meanC","ave.nat.cov","ave.inv.cov","ave.hum.cov","ave.rcg.cov")
k = 1
for (y in yrs) {
  for (i in 1:n.t) {
    t = trts[i]
    for (j in 1:n.p) {
      plt.ind1 = which(df.ave.cov$year == y & df.ave.cov$trt == t & df.ave.cov$plot == j)
      plt.ind2 = which(all.cov$year == y & all.cov$Plot == paste(t, j, sep=""))
      df.plot$year[k] = y
      df.plot$trt[k] = t
      df.plot$trt.full[k] = trt.names[i]
      df.plot$plot[k] = j
      df.plot$richness[k] = length(df.ave.cov$spp[plt.ind1])
      df.plot$meanC[k] = mean(all.cov$C.value[plt.ind2])
      df.plot$ave.nat.cov[k] = sum(df.ave.cov$ave.cov[plt.ind1]*df.ave.cov$native[plt.ind1])
      df.plot$ave.inv.cov[k] = sum(df.ave.cov$ave.cov[plt.ind1]*(1-df.ave.cov$native[plt.ind1]))
      rcg.ind = which(df.ave.cov$spp[plt.ind1] == "Phalaris arundinacea")
      hum.ind = which(df.ave.cov$spp[plt.ind1] == "Humulus japonicus")
      df.plot$ave.rcg.cov[k] = sum(df.ave.cov$ave.cov[plt.ind1][rcg.ind])
      df.plot$ave.hum.cov[k] = sum(df.ave.cov$ave.cov[plt.ind1][hum.ind])
      k = k + 1
    }
  }
}


# richness boxplot by year
ggplot(df.plot, aes(x=factor(trt.full,level=trt.names), y=richness, fill=factor(year))) + 
  geom_boxplot() + geom_point(position=position_dodge(width=0.75),aes(group=year)) +
  labs(x="Treatment", y="Species Richness") +
  theme(legend.title = element_blank())

# run anova on richness
df.plot$year = as.factor(df.plot$year)
df.plot$trt = as.factor(df.plot$trt)
aov.richness = aov(richness~trt*year, data=df.plot)
summary(aov.richness)

# run TSD on richness
t.hsd.rich = TukeyHSD(aov.richness)
par(mar=c(5,14,4,1))
plot(TukeyHSD(aov.richness, conf.level=.95), las=1)

# meanC boxplot by year
ggplot(df.plot, aes(x=factor(trt.full,level=trt.names), y=meanC, fill=factor(year))) + 
  geom_boxplot() + geom_point(position=position_dodge(width=0.75),aes(group=year)) +
  labs(x="Treatment", y="Mean C") +
  theme(legend.title = element_blank())

# anova on mean C
aov.meanC = aov(meanC~trt*year, data=df.plot)
summary(aov.meanC)

# compare invasive v. native cover by year
ggplot(df.plot, aes(x=factor(trt.full,level=trt.names), y=ave.nat.cov, fill=factor(year))) + 
  geom_boxplot() + geom_point(position=position_dodge(width=0.75),aes(group=year)) +
  labs(x="Treatment", y="Average Native Cover") +
  theme(legend.title = element_blank())

ggplot(df.plot, aes(x=factor(trt.full,level=trt.names), y=ave.inv.cov, fill=factor(year))) + 
  geom_boxplot() + geom_point(position=position_dodge(width=0.75),aes(group=year)) +
  labs(x="Treatment", y="Average Invasive Cover") +
  theme(legend.title = element_blank())

aov.inv.cover = aov(ave.inv.cov~trt*year, data=df.plot)
summary(aov.inv.cover)

aov.nat.cover = aov(ave.nat.cov~trt*year, data=df.plot)
summary(aov.nat.cover)

### calculate average cover by year

# calculate average cover for each species per plot

df.ave.cov = data.frame(matrix(nrow=0, ncol=7))
colnames(df.ave.cov) = c("year","trt","trt.full","plot","spp","ave.cov","native")
for (i in 1:n.y) {
  y = yrs[i]
  for (j in 1:n.t) {
    t = trts[j]
    for (k in 1:n.p) {
      plt.ind = which(all.cov$year == y & all.cov$Plot == paste(t, k, sep=""))
      uni.sp = unique(all.cov$Spp[plt.ind])
      n.sp = length(uni.sp)
      new.df = data.frame(matrix(nrow=n.sp, ncol=7))
      colnames(new.df) = c("year","trt","trt.full","plot","spp","ave.cov","native")
      new.df$year = y
      new.df$trt = t
      new.df$trt.full = trt.names[j]
      new.df$plot = k
      new.df$spp = uni.sp
      l = 1
      for (s in uni.sp) {
        s.ind = which(all.cov$Spp[plt.ind] == s)
        new.df$ave.cov[l] = sum(all.cov$Cover[plt.ind][s.ind])/5
        new.df$native[l] = all.cov$Native[plt.ind][s.ind][1]
        l = l + 1
      }
      df.ave.cov = rbind(df.ave.cov,new.df)
    }
  }
}

# compare species average cover across plots
ggplot(df.ave.cov, aes(fill=spp, y=ave.cover, x=trt)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x="Treatment", y="Average cover", title="2022")


################################################################################
# PCA analysis on cover data

# df for cover by species
uni.sp = unique(df.ave.cov$spp)
n.sp = length(uni.sp)
df.cover.sp = data.frame(matrix(nrow=n.y*n.t*n.p, ncol=n.sp+4))
colnames(df.cover.sp) = c("year","trt","trt.full","plot",uni.sp)
df.cover.sp[,uni.sp] = 0
k = 1
for (n in 1:n.y) {
  for (i in 1:n.t) {
    for (j in 1:n.p) {
      plt.ind = which(df.ave.cov$year == yrs[n] & df.ave.cov$trt == trts[i] & df.ave.cov$plot == j)
      sp.ind = which(colnames(df.cover.sp) %in% df.ave.cov$spp[plt.ind])
      df.cover.sp$year[k] = yrs[n]
      df.cover.sp$trt[k] = trts[i]
      df.cover.sp$trt.full[k] = trt.names[i]
      df.cover.sp$plot[k] = j
      df.cover.sp[k,sp.ind] = df.ave.cov$ave.cov[plt.ind]
      k = k + 1
    }
  }
}

# separate 2013 & 2022 data
df.cover.2013 = df.cover.sp[which(df.cover.sp$year==2013),]
df.cover.2022 = df.cover.sp[which(df.cover.sp$year==2022),]

# remove zero columns
zero.2013 = which(apply(df.cover.2013[,-seq(1,4)], 2, FUN = function(x) sum(x)) == 0) + 4
zero.2022 = which(apply(df.cover.2022[,-seq(1,4)], 2, FUN = function(x) sum(x)) == 0) + 4
df.new.2013 = df.cover.2013[,-zero.2013]
df.new.2022 = df.cover.2022[,-zero.2022]

# run pcas
pca.cov.2013 = prcomp(df.new.2013[,-seq(1,4)], center=TRUE, scale=TRUE)
pca.cov.2022 = prcomp(df.new.2022[,-seq(1,4)], center=TRUE, scale=TRUE)
#pca.C$rotation[,1:3]
#pcs = data.frame(as.matrix(df.c[,-c(1,2,3)]) %*% as.matrix(pca.C$rotation))
#pcs$trt = df.c$trt
#df.c$year = factor(df.c$year)

# 2013
thres=0
sp.thres = which(abs(pca.cov.2013$rotation[,"PC1"])>thres | abs(pca.cov.2013$rotation[,"PC2"]) >thres)
n.sp.2013 = length(colnames(df.new.2013))-4
lab.sizes = rep(0,n.sp.2013)
lab.sizes[sp.thres] = 2
lab.cols = rep("white",n.sp.2013)
lab.cols[sp.thres] = "black"
autoplot(pca.cov.2013, data=df.cover.2013, colour="trt.full", size=3,
         loadings=TRUE, loadings.label=TRUE, loadings.label.size=lab.sizes,
         frame = TRUE, repel=TRUE, loadings.label.colour=lab.cols) + 
         labs(colour="Treatment", fill='Treatment', title="2013")

# 2022
thres=0.1
sp.thres = which(abs(pca.cov.2022$rotation[,"PC1"])>thres | abs(pca.cov.2022$rotation[,"PC2"]) >thres)
n.sp.2022 = length(colnames(df.new.2022))-4
lab.sizes = rep(0,n.sp.2022)
lab.sizes[sp.thres] = 2
lab.cols = rep("white",n.sp.2022)
lab.cols[sp.thres] = "black"
autoplot(pca.cov.2022, data=df.cover.2022, colour="trt.full", size=3,
         loadings=TRUE, loadings.label=TRUE, loadings.label.size=lab.sizes,
         frame = TRUE, repel=TRUE, loadings.label.colour=lab.cols) + 
         labs(colour="Treatment", fill='Treatment', title="2022")

# both years together
pca.cov = prcomp(df.cover.sp[,-seq(1,4)], center=TRUE, scale=TRUE)

thres=0.1
sp.thres = which(abs(pca.cov$rotation[,"PC1"])>thres | abs(pca.cov$rotation[,"PC2"]) >thres)
n.sp = length(colnames(df.cover.sp))-4
lab.sizes = rep(0,n.sp)
lab.sizes[sp.thres] = 2
lab.cols = rep("white",n.sp)
lab.cols[sp.thres] = "black"
df.cover.sp$year = as.factor(df.cover.sp$year)
autoplot(pca.cov, data=df.cover.sp, colour="trt.full", shape="year", size=3,
         loadings=TRUE, loadings.label=TRUE, loadings.label.size=lab.sizes,
         repel=TRUE, loadings.label.colour=lab.cols) + 
         labs(colour="Treatment", fill='Treatment', title="2013 & 2022")

################################################################################
## calculate species/family-level replacement by treatment

# make dataframe for number of species, families, and orders replaced
df.cols = c("year","trt",
            "n.sp","n.sp.rep","n.sp.lost","n.sp.gain",
            "n.fam","n.fam.rep","n.fam.lost","n.fam.gain",
            "n.ord","n.ord.rep","n.ord.lost","n.ord.gain")
df.rep = data.frame(matrix(nrow=n.y*n.t, ncol=length(df.cols)))
colnames(df.rep) = df.cols
tax.levels = c("sp","fam","ord")
rep.list = list()
for (i in 1:n.y) {
  rep.list[[i]] = list()
  for (j in 1:n.t) {
    rep.list[[i]][[trts[j]]] = list()
    for (k in 1:3) {rep.list[[i]][[trts[j]]][[tax.levels[k]]] = list("rep"=list(),"lost"=list(),"gain"=list())}
  }
}

k = 1
for (i in 1:n.y) {
  y = yrs[i]
  ref.id = which(all.cov$year == y & all.cov$trt == "R")
  ref.sp = unique(all.cov[ref.id,"spp"])
  ref.fam = unique(all.cov[ref.id,"family"])
  ref.ord = unique(all.cov[ref.id,"order"])
  for (j in 1:n.t) {
      t = trts[j]
      trt.id = which(all.cov$year == y & all.cov$trt == t)
      plt.sp = unique(all.cov[trt.id,"spp"])
      plt.fam = unique(all.cov[trt.id,"family"])
      plt.ord = unique(all.cov[trt.id,"order"])
      df.rep[k,c("year","trt")] = c(y, t)
      
      # species
      df.rep[k,"n.sp"] = length(plt.sp)
      df.rep[k,"n.sp.rep"] = sum(plt.sp %in% ref.sp)
      df.rep[k,"n.sp.lost"] = length(ref.sp) - sum(ref.sp %in% plt.sp)
      df.rep[k,"n.sp.gain"] = length(plt.sp) - sum(plt.sp %in% ref.sp)
      rep.list[[i]][[trts[j]]][["sp"]][["rep"]] = plt.sp[which(plt.sp %in% ref.sp)]
      rep.list[[i]][[trts[j]]][["sp"]][["lost"]] = ref.sp[-which(ref.sp %in% plt.sp)]
      rep.list[[i]][[trts[j]]][["sp"]][["gain"]] = plt.sp[-which(plt.sp %in% ref.sp)]
      
      # families
      df.rep[k,"n.fam"] = length(plt.fam)
      df.rep[k,"n.fam.rep"] = sum(plt.fam %in% ref.fam)
      df.rep[k,"n.fam.lost"] = length(ref.fam) - sum(ref.fam %in% plt.fam)
      df.rep[k,"n.fam.gain"] = length(plt.fam) - sum(plt.fam %in% ref.fam)
      rep.list[[i]][[trts[j]]][["fam"]][["rep"]] = plt.fam[which(plt.fam %in% ref.fam)]
      rep.list[[i]][[trts[j]]][["fam"]][["lost"]] = ref.fam[-which(ref.fam %in% plt.fam)]
      rep.list[[i]][[trts[j]]][["fam"]][["gain"]] = plt.fam[-which(plt.fam %in% ref.fam)]
      
      # orders
      df.rep[k,"n.ord"] = length(plt.ord)
      df.rep[k,"n.ord.rep"] = sum(plt.ord %in% ref.ord)
      df.rep[k,"n.ord.lost"] = length(ref.ord) - sum(ref.ord %in% plt.ord)
      df.rep[k,"n.ord.gain"] = length(plt.ord) - sum(plt.ord %in% ref.ord)
      rep.list[[i]][[trts[j]]][["ord"]][["rep"]] = plt.ord[which(plt.ord %in% ref.ord)]
      rep.list[[i]][[trts[j]]][["ord"]][["lost"]] = ref.ord[-which(ref.ord %in% plt.ord)]
      rep.list[[i]][[trts[j]]][["ord"]][["gain"]] = plt.ord[-which(plt.ord %in% ref.ord)]
      k = k + 1
  }
}

## plot number of species/families/orders replaced, gained, and lost

# species, family, and order in separate plots
my.colors <- c("cornflowerblue","indianred2","olivedrab3")
df.rep.melt = melt(df.rep, id.vars = c("year","trt"), variable.name = "var", value.name="value")
sp.vars = c("n.sp.rep","n.sp.gain","n.sp.lost")
sp.ind = which(df.rep.melt$var %in% sp.vars & df.rep.melt$trt %in% trts[1:5])
sp.rep = ggplot(df.rep.melt[sp.ind,], aes(x=trt, y=value, fill=var)) + 
                geom_bar(position="dodge", stat="identity") + facet_grid(rows=vars(year)) +
                scale_x_discrete(labels=trt.names) + labs(x="",y="Count",fill="",title="Species")+
                scale_fill_manual(labels=c("Replaced","Lost","Gained"),values=my.colors)

fam.vars = c("n.fam.rep","n.fam.gain","n.fam.lost")
fam.ind = which(df.rep.melt$var %in% fam.vars & df.rep.melt$trt %in% trts[1:5])
fam.rep = ggplot(df.rep.melt[fam.ind,], aes(x=trt, y=value, fill=var)) + 
                 geom_bar(position="dodge", stat="identity") + facet_grid(rows=vars(year)) +
                 scale_x_discrete(labels=trt.names) + labs(x="",y="",fill="",title="Families")+
                 scale_fill_manual(labels=c("Replaced","Lost","Gained"),values=my.colors)

ord.vars = c("n.ord.rep","n.ord.gain","n.ord.lost")
ord.ind = which(df.rep.melt$var %in% ord.vars & df.rep.melt$trt %in% trts[1:5])
ord.rep = ggplot(df.rep.melt[ord.ind,], aes(x=trt, y=value, fill=var)) + 
                 geom_bar(position="dodge", stat="identity") + facet_grid(rows=vars(year)) +
                 scale_x_discrete(labels=trt.names) + labs(x="",y="",fill="",title="Orders") + 
                 scale_fill_manual(labels=c("Replaced","Lost","Gained"),values=my.colors)
grid_arrange_shared_legend(sp.rep, fam.rep, ord.rep, ncol=3, nrow=1)

# replaced, gained, and lost in separate plots
rep.vars = c("n.sp.rep","n.fam.rep","n.ord.rep")
rep.ind = which(df.rep.melt$var %in% rep.vars & df.rep.melt$trt %in% trts[1:5])
p.rep = ggplot(df.rep.melt[rep.ind,], aes(x=trt, y=value, fill=var)) + 
               geom_bar(position="dodge", stat="identity") + facet_grid(rows=vars(year)) +
                        scale_x_discrete(labels=trt.names) + labs(x="",y="Count",fill="",title="Replaced")+
                        scale_fill_manual(labels=c("Species","Family","Order"),values=my.colors)

gain.vars = c("n.sp.gain","n.fam.gain","n.ord.gain")
gain.ind = which(df.rep.melt$var %in% gain.vars & df.rep.melt$trt %in% trts[1:5])
p.gain = ggplot(df.rep.melt[gain.ind,], aes(x=trt, y=value, fill=var)) + 
                geom_bar(position="dodge", stat="identity") + facet_grid(rows=vars(year)) +
                scale_x_discrete(labels=trt.names) + labs(x="",y="",fill="",title="Gained")+
                scale_fill_manual(labels=c("Species","Family","Order"),values=my.colors)

lost.vars = c("n.sp.lost","n.fam.lost","n.ord.lost")
lost.ind = which(df.rep.melt$var %in% lost.vars & df.rep.melt$trt %in% trts[1:5])
p.lost = ggplot(df.rep.melt[lost.ind,], aes(x=trt, y=value, fill=var)) + 
                geom_bar(position="dodge", stat="identity") + facet_grid(rows=vars(year)) +
                scale_x_discrete(labels=trt.names) + labs(x="",y="",fill="",title="Lost") + 
                scale_fill_manual(labels=c("Species","Family","Order"),values=my.colors)
grid_arrange_shared_legend(p.rep, p.gain, p.lost, ncol=3, nrow=1)

# plot family replacement by treatment
fam.2013 = unique(all.cov[all.cov$year==2013,"family"]); n.fam.2013 = length(fam.2013)
fam.2022 = unique(all.cov[all.cov$year==2022,"family"]); n.fam.2022 = length(fam.2022)
n.fam = n.fam.2013 + n.fam.2022
df.fam = data.frame(matrix(nrow=n.fam,ncol=n.t+1))
colnames(df.fam) = c("year","family",trts[1:5])
k = 1
for (i in 1:n.y) {
  y = yrs[i]
  y.id = which(all.cov$year == y)
  fam.yr = unique(all.cov[y.id,"family"])
  for (f in fam.yr) {
    df.fam[k, c("year","family")] = c(y, f)
    for (j in 1:5) {
      t = trts[j]
      if  (f %in% rep.list[[i]][[trts[j]]][["fam"]][["gain"]]) {
        df.fam[k, t] = 1
      } else if (f %in% rep.list[[i]][[trts[j]]][["fam"]][["rep"]]) {
        df.fam[k, t] = 2
      } else if (f %in% rep.list[[i]][[trts[j]]][["fam"]][["lost"]]) {
        df.fam[k, t] = 3        
      } else {
        df.fam[k, t] = 0
      }
    }
    k = k + 1
  }
}

df.fam.2022 = df.fam[which(df.fam$year==2022),c("family",trts[1:5])]
row.totals = rowSums(df.fam.2022[,trts[1:5]])
row.sort = sort(row.totals, index.return=TRUE, decreasing=F)
fam.melt = melt(df.fam.2022, id.vars = c("family"), variable.name = "trt", value.name="value")
fam.melt$family <- factor(fam.melt$family, levels=(fam.melt$family)[row.sort$ix])
group.colors <- c("white","olivedrab3", "cornflowerblue", "indianred2")
ggplot(fam.melt, aes(y=trt, x=family, fill=factor(value))) + geom_tile() + labs(fill="") +
       scale_fill_manual(values=group.colors,labels=c("None","Gained","Replaced","Lost")) +
       scale_y_discrete(labels=trt.names) + labs(y="",x="Family",title="2022")+
       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# plot species replacement by treatment
sp.2013 = unique(all.cov[all.cov$year==2013,"spp"]); n.sp.2013 = length(sp.2013)
sp.2022 = unique(all.cov[all.cov$year==2022,"spp"]); n.sp.2022 = length(sp.2022)
n.sp = n.sp.2013 + n.sp.2022
df.sp = data.frame(matrix(nrow=n.sp, ncol=n.t+1))
colnames(df.sp) = c("year","spp",trts[1:5])
k = 1
for (i in 1:n.y) {
  y = yrs[i]
  y.id = which(all.cov$year == y)
  sp.yr = unique(all.cov[y.id,"spp"])
  for (s in sp.yr) {
    df.sp[k, c("year","spp")] = c(y, s)
    for (j in 1:5) {
      t = trts[j]
      if  (s %in% rep.list[[i]][[trts[j]]][["sp"]][["gain"]]) {
        df.sp[k, t] = 1
      } else if (s %in% rep.list[[i]][[trts[j]]][["sp"]][["rep"]]) {
        df.sp[k, t] = 2
      } else if (s %in% rep.list[[i]][[trts[j]]][["sp"]][["lost"]]) {
        df.sp[k, t] = 3        
      } else {
        df.sp[k, t] = 0
      }
    }
    k = k + 1
  }
}

df.sp.2022 = df.sp[which(df.sp$year==2022), c("spp",trts[1:5])]
row.totals = rowSums(df.sp.2022[,trts[1:5]])
row.sort = sort(row.totals, index.return=TRUE, decreasing=F)
sp.melt = melt(df.sp.2022, id.vars = c("spp"), variable.name = "trt", value.name="value")
sp.melt$spp <- factor(sp.melt$spp, levels=(sp.melt$spp)[row.sort$ix])
ggplot(sp.melt, aes(y=trt, x=spp, fill=factor(value))) + geom_tile() + labs(fill="") +
       scale_fill_manual(values=group.colors,labels=c("None","Gained","Replaced","Lost")) +
       scale_y_discrete(labels=trt.names) + labs(y="",x="Species",title="2022") +
       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# plot order replacement by treatment
ord.2013 = unique(all.cov[all.cov$year==2013,"order"]); n.ord.2013 = length(ord.2013)
ord.2022 = unique(all.cov[all.cov$year==2022,"order"]); n.ord.2022 = length(ord.2022)
n.ord = n.ord.2013 + n.ord.2022
df.ord = data.frame(matrix(nrow=n.ord, ncol=n.t+1))
colnames(df.ord) = c("year","order",trts[1:5])
k = 1
for (i in 1:n.y) {
  y = yrs[i]
  y.id = which(all.cov$year == y)
  ord.yr = unique(all.cov[y.id,"order"])
  for (o in ord.yr) {
    df.ord[k, c("year","order")] = c(y, o)
    for (j in 1:5) {
      t = trts[j]
      if  (o %in% rep.list[[i]][[trts[j]]][["ord"]][["gain"]]) {
        df.ord[k, t] = 1
      } else if (o %in% rep.list[[i]][[trts[j]]][["ord"]][["rep"]]) {
        df.ord[k, t] = 2
      } else if (o %in% rep.list[[i]][[trts[j]]][["ord"]][["lost"]]) {
        df.ord[k, t] = 3        
      } else {
        df.ord[k, t] = 0
      }
    }
    k = k + 1
  }
}

df.ord.2022 = df.ord[which(df.ord$year==2022), c("order",trts[1:5])]
row.totals = rowSums(df.ord.2022[,trts[1:5]])
row.sort = sort(row.totals, index.return=TRUE, decreasing=F)
ord.melt = melt(df.ord.2022, id.vars = c("order"), variable.name = "trt", value.name="value")
ord.melt$order <- factor(ord.melt$order, levels=(ord.melt$order)[row.sort$ix])
ggplot(ord.melt, aes(y=trt, x=order, fill=factor(value))) + geom_tile() + labs(fill="") +
       scale_fill_manual(values=group.colors,labels=c("None","Gained","Replaced","Lost")) +
       scale_y_discrete(labels=trt.names) + labs(y="",x="Order",title="2022") +
       theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
