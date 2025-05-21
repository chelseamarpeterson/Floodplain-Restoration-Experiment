path_to_tree_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Public-Repo/Tree_Analysis"
path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Public-Repo"
setwd(path_to_tree_folder)

library(rethinking)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(RColorBrewer)

# read in posterior intervals
int.df = read.csv("Posteriors/CarbonRichness_Means_Intervals_5Chains.csv")

# variables
var.order = c("n.total","N.herb","N.tree",
              "MAOM.C","POM.C","SOC","TIC","TC",
              "int.fwd.carbon","cwd.carbon",
              "snag.carbonmin","dead.stem.carbonmin",
              "live.woody.c.stock","live.stem.carbon",
              "FWD.C.Mg.ha","HL.C.Mg.ha","total.C.Mg.ha",
              "total.live.carbon","total.dead.carbon",
              "total.veg.carbon","total.ecosystem.carbon")
var.labels = c("Total richness","Herbaceous species","Tree species",
               "MAOM-C (< 53 \U00B5m)","POM-C (\u2265 53 \U00B5m)",
               "TOC","Carbonates","TC",
               "Fine woody debris (2.5-7.5 cm)","Coarse woody debris (\u2265 7.6 cm)",
               "Standing dead trees (\u2265 2.5 cm)","Dead stems (< 2.5 cm)",
               "Live trees (\u2265 2.5 cm)","Live stems (< 2.5 cm)",
               "Fine woody debris (< 2.5 cm)","Herbaceous litter","Herbaceous biomass",
               "Total live vegetation","Total dead vegetation",
               "Total vegetation","Total ecosystem")

# plot intervals for all variables
ggplot(int.df, aes(y=factor(trt, levels=trt.names), x=mean)) +
       geom_errorbar(aes(xmin=X5, xmax=X95), width=0.25, color="black") +
       geom_errorbar(aes(xmin=X25, xmax=X75), width=0.25, color="blue") +
       geom_point() + labs(y="",x="Posterior mean") +
       facet_wrap(.~factor(var, levels=var.labels), 
                  ncol=3, scales="free_x")

# print intervals
i = 21 # left off on 13
print(cbind(int.df[which(int.df$var == var.labels[i]),c("var","trt")],
            round(int.df[which(int.df$var == var.labels[i]),c("mean","X5","X95")],1)))
print(cbind(int.df[which(int.df$var == var.labels[i]),c("var","trt")],
            round(int.df[which(int.df$var == var.labels[i]),c("mean","X5","X95")],2)))
print(cbind(int.df[which(int.df$var == var.labels[i]),c("var","trt")],
            round(int.df[which(int.df$var == var.labels[i]),c("mean","X5","X95")],3)))

## plot stacked means of biomass variables
SOC.est.df = data.frame(matrix(nrow=3,ncol=2))
SOC.est.df$`IPCC SOC` = c("Natural wetland","Revegetated cropland","Annual crops")
SOC.est.df$value = c(88.0, 77.2, 55.9)

ABC.est.df = data.frame(matrix(nrow=2,ncol=2))
ABC.est.df$`IPCC woody biomass` = c("Restored temperate forest","Mature temperate forest")
ABC.est.df$value = c(46.1, 62.4)

dead.est.df = data.frame(matrix(nrow=2,ncol=2))
dead.est.df$`IPCC dead wood and litter` = c("Restored temperate forest","Mature temperate forest")
dead.est.df$value = c(14.2, 29.7)

TC.est.df = data.frame(matrix(nrow=2,ncol=2))
TC.est.df$`IPCC total organic C` = c("Restored forested wetland","Mature forested wetland")
TC.est.df$value = c(132.5, 180.1)

# total ecosystem
stack.vars1 = c("TIC","Total dead vegetation","TOC","Total live vegetation")
v.palette <- brewer.pal(11,"RdYlGn")[c(7,9,11)]
s.palette <- c(brewer.pal(9,"Greys")[5],
               brewer.pal(11,"BrBG")[c(2,1)])
d.palette <- brewer.pal(9,"YlOrBr")[seq(3,8)]
p1 = ggplot(int.df[which(int.df$var %in% stack.vars1),], 
            aes(y=factor(trt, levels=trt.names), x=mean, 
                fill=factor(var, levels=stack.vars1))) +
            geom_bar(stat="identity",position="stack") + 
            labs(fill="",y="",
                 x="Posterior mean (Mg/ha)",title="") +
            geom_vline(data=TC.est.df, aes(xintercept=value, 
                                           color=`IPCC total organic C`,
                                           linetype=`IPCC total organic C`), 
                       linewidth=1) +
            scale_color_manual(values=c("blue","red")) +
            scale_linetype_manual(values=c("dashed","solid")) +
            scale_fill_manual(values=c(brewer.pal(9,"Greys")[5],d.palette[4],s.palette[3],v.palette[3]),
                              labels=c("Carbonates","Litter and woody debris","Soil organic matter","Living biomass")) +
            theme(text=element_text(size=14), legend.key.size=unit(0.7,'cm'))
p1

#  live vegetation plot
stack.vars2 = c("Herbaceous biomass","Live stems (< 2.5 cm)","Live trees (\u2265 2.5 cm)")
ipcc.vars2 = c("Restored temperate forest","Mature temperate forest")
p2 = ggplot(int.df[which(int.df$var %in% stack.vars2),], 
            aes(y=factor(trt, levels=trt.names), x=mean, 
                fill=factor(var, levels=stack.vars2))) +
            geom_bar(stat="identity",position="stack") + 
            labs(fill="",y="",
                 x="Posterior mean (Mg/ha)",title="") + 
            geom_vline(data=ABC.est.df,
                       aes(xintercept=value, 
                           color=factor(`IPCC woody biomass`,levels=ipcc.vars2),
                           linetype=factor(`IPCC woody biomass`,levels=ipcc.vars2)), 
                       linewidth=1) +
            scale_color_manual(values=c("red","blue")) +
            scale_linetype_manual(values=c("solid","dashed")) + 
            scale_fill_manual(values=c("olivedrab3","olivedrab4",v.palette[3])) +
            guides(color=guide_legend(title="IPCC woody biomass"),
                   linetype=guide_legend(title="IPCC woody biomass")) +
            theme(text=element_text(size=14), legend.key.size=unit(0.7,'cm'))
p2

# debris plot
stack.vars3 = c("Herbaceous litter","Fine woody debris (< 2.5 cm)",
                "Fine woody debris (2.5-7.5 cm)","Coarse woody debris (\u2265 7.6 cm)",
                "Dead stems (< 2.5 cm)","Standing dead trees (\u2265 2.5 cm)")
ipcc.vars3 = c("Restored temperate forest","Mature temperate forest")
p3 = ggplot(int.df[which(int.df$var %in% stack.vars3),], 
            aes(y=factor(trt, levels=trt.names), x=mean, 
                fill=factor(var, levels=stack.vars3))) +
            geom_bar(stat="identity",position="stack") + 
            labs(fill="",y="",
                 x="Posterior mean (Mg/ha)",title="") + 
            geom_vline(data=dead.est.df, 
                       aes(xintercept=value,
                           color=factor(`IPCC dead wood and litter`,levels=ipcc.vars3),
                           linetype=factor(`IPCC dead wood and litter`,levels=ipcc.vars3)), 
                       linewidth=1) +
            scale_color_manual(values=c("red","blue")) +
            scale_linetype_manual(values=c("solid","dashed")) + 
            scale_fill_manual(values=d.palette) + 
            guides(color=guide_legend(title="IPCC dead wood and litter"),
                   linetype=guide_legend(title="IPCC dead wood and litter")) +
            theme(text=element_text(size=14), legend.key.size=unit(0.7,'cm'))
p3

# soil plot
stack.vars4 = c("Carbonates","POM-C (\u2265 53 \U00B5m)","MAOM-C (< 53 \U00B5m)")
ipcc.vars4 = c("Natural wetland","Revegetated cropland","Annual crops")
p4 = ggplot(int.df[which(int.df$var %in% stack.vars4),], 
            aes(y=factor(trt, levels=trt.names), 
                x=mean, fill=factor(var, levels=stack.vars4))) +
            geom_bar(stat="identity",position="stack") + 
            scale_fill_manual(values=s.palette) +
            labs(fill="",y="",
                 x="Posterior mean (Mg/ha)",title="") +
            geom_vline(data=SOC.est.df, 
                       aes(xintercept=value, 
                           color=factor(`IPCC SOC`,levels=ipcc.vars4),
                           linetype=factor(`IPCC SOC`,levels=ipcc.vars4)), 
                       linewidth=1) +
            scale_color_manual(values=c("blue","purple","red")) + 
            scale_linetype_manual(values=c("dashed","dotted","solid")) +
            guides(color=guide_legend(title="IPCC soil organic carbon"),
                   linetype=guide_legend(title="IPCC soil organic carbon")) +
            theme(text=element_text(size=14), legend.key.size=unit(0.7,'cm'))
p4
p.c.stocks = (p4+theme(plot.margin=unit(c(0,32,0,0),"pt"))+p2)/(p3+theme(plot.margin=unit(c(0,2,0,0),"pt"))+p1)
p.c.stocks
setwd(path_to_repo)
ggsave("Figures/Soil_Veg_Ecosystem_Cstocks.jpeg", 
       plot = p.c.stocks, width = 38, height = 20, units="cm")

# richness plot with error bars
stack.vars5 = c("Total richness","Herbaceous species","Tree species")
ggplot(int.df[which(int.df$var %in% stack.vars5),], 
       aes(y=factor(trt, levels=trt.names), 
           x=mean, 
           fill=factor(var, levels=stack.vars5))) +
       geom_bar(stat="identity", position=position_dodge()) +
       #geom_errorbarh(aes(xmin=`5`, xmax=`95`,y=factor(t,levels=trt.names)), 
       #                height=0.2, position=position_dodge(0.9)) + 
       labs(fill="",y="",x="Posterior mean",title="Species richness") + 
       scale_fill_manual(values=c("sienna","darkolivegreen4","purple")) +
       theme(text = element_text(size=14))

################################################################################
### statistical analysis on carbon stock data by species

## CWD C stocks

# read in CWD data
cwd.melt = read.csv("CWD_Carbon_Stocks_By_Species.csv", header=T)

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
                area = cwd.melt$cwd.carbon.sum/mean(cwd.melt$cwd.carbon.sum)) 
mean.cwd = mean(cwd.melt$cwd.carbon.sum)

# fit simple model
set.seed(12)
m.cwd.sp = ulam(alist(area ~ dnorm(mu, sigma),
                      log(mu) <- a[tid, spid],
                      matrix[tid, spid]:a ~ dnorm(0, 1),
                      sigma ~ dexp(1)),
                data = d.cwd.sp, chains=10)

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
int.dfwd = data.frame(matrix(nrow=n.trt*n.spp.cwd, ncol=9))
colnames(int.dfwd) = c("tid","spid","trt","spp","mean","5","95","25","75")
k = 1
for (i in 1:n.trt) {
  for (j in 1:n.spp.cwd) {
    a.mean = mean(a.cwd[,i,j], 2)
    a.hpdi.90 = HPDI(a.cwd[,i,j], prob=0.90)
    a.hpdi.50 = HPDI(a.cwd[,i,j], prob=0.50)
    int.dfwd[k,"tid"] = i
    int.dfwd[k,"spid"] = j
    int.dfwd[k,"trt"] = trt.names[i]
    int.dfwd[k,"spp"] = cwd.spp[j]
    int.dfwd[k,"mean"] = a.mean
    int.dfwd[k,"5"] = a.hpdi.90[1]
    int.dfwd[k,"95"] = a.hpdi.90[2]
    int.dfwd[k,"25"] = a.hpdi.50[1]
    int.dfwd[k,"75"] = a.hpdi.50[2]
    k = k + 1
  }
}

# plot posteriors
ggplot(int.dfwd, aes(x=factor(trt, levels=trt.names), 
                        fill=factor(spp, levels=cwd.spp), y=mean)) + 
       geom_bar(stat = "identity") + 
       labs(x="", y="Posterior mean of CWD [>= 7.6 cm] area (m2/km)", 
            fill="Species")

## snag C stocks

# read in snag data
snag.melt = read.csv("Snag_Carbon_Stocks_By_Species.csv", header=T)

# convert treatment and species into id
snag.melt$tid = as.numeric(factor(snag.melt$trt))
snag.melt$spid = as.numeric(factor(snag.melt$species))

# make simplified data list
d.snag.sp = list(tid = snag.melt$tid, 
                 spid = snag.melt$spid,
                 area = snag.melt$snag.carbonmin/mean(snag.melt$snag.carbonmin)) 
snag.mean = mean(snag.melt$snag.carbonmin)
#hist(log(d.snag.sp$area))

# fit simple model
set.seed(12)
m.snag.sp = ulam(alist(area ~ dnorm(mu, sigma),
                       log(mu) <- a[tid, spid],
                       matrix[tid, spid]:a ~ dnorm(0, 1),
                       sigma ~ dexp(1)),
                 data = d.snag.sp, chains=10)

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

# plot posteriors
ggplot(post.df.snag, aes(x=factor(trt, levels=trt.names), 
                         fill=factor(spp, levels=snag.spp), y=mean)) + 
       geom_bar(stat = "identity") + 
       labs(x="", y="Posterior mean of snag basal area (m2/ha)", fill="Species")

## large woody biomass C stocks

# prep dataframe
sp.C.dat = read.csv("Biomass_C_stocks_by_species.csv")
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
             bmC = sp.C.melt$value/mean(sp.C.melt$value)) 
mean.bmC = mean(sp.C.melt$value)

# fit simple model
set.seed(12)
m.bmC = ulam(alist(bmC ~ dnorm(mu, sigma),
                   log(mu) <- a[tid, spid],
                   matrix[tid, spid]:a ~ dnorm(0, 1),
                   sigma ~ dexp(1)),
             data = d.bmC, chains=10)

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

# plot posteriors
all.spp.sort = sort(levels(factor(sp.C.melt$spp)))
ggplot(post.df.woodC, aes(y=factor(trt, levels=trt.names), 
                           fill=factor(spp, levels=all.spp.sort), x=mean)) + 
       geom_bar(stat = "identity") + 
       labs(y="", x="Posterior mean of live woody C stocks (Mg/ha)", 
            fill="Species")

### statistical analyses on biomass C stocks by species ----

# read in biomass C stock estimates
bm.dat = read.csv("Clean_Cstocks_Richness.csv")

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

# plot histogram
par(mfrow=c(1,1))
hist(log(bm.sp.C.melt$value/sd(bm.sp.C.melt$value)))

# make simple dataframe
d.bm = list(tid = as.integer(bm.sp.C.melt$tid),
            spid = as.integer(bm.sp.C.melt$spid),
            bm = bm.sp.C.melt$value/mean(bm.sp.C.melt$value))

# fit model
set.seed(12)
m.bm.sp = ulam(alist(bm ~ dnorm(mu, sigma),
                     log(mu) <- a[tid, spid],
                     matrix[tid, spid]:a ~ dnorm(0, 1),
                     sigma ~ dexp(1)),
               data = d.bm, chains=10)

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
ggplot(post.df.bm.sp, aes(y=factor(trt, levels=trt.names), 
                          fill=factor(spp, levels=bm.sp), x=mean)) + 
       geom_bar(stat = "identity") + 
       labs(y="", x="Posterior mean C stock (Mg/ha)", 
            title = "Herbaceous biomass",
            fill="Species group") + 
       theme(text = element_text(linewidth=14))

## make combined color plot with live woody, snag/CWD area, and biomass stocks by species
post.df.woodC$type = "Live trees (\u2265 2.5 cm)"
post.df.snag$type = "Standing dead trees (\u2265 2.5 cm)"
int.dfwd$type = "Coarse woody debris (\u2265 7.6 cm)"
types = c("Live trees (\u2265 2.5 cm)",
          "Standing dead trees (\u2265 2.5 cm)",
          "Coarse woody debris (\u2265 7.6 cm)")
post.df.sp.all = rbind(post.df.woodC, post.df.snag, int.dfwd)
post.df.sp.all$spp = trimws(post.df.sp.all$spp)
all.sp = sort(unique(post.df.sp.all$spp))
ggplot(post.df.sp.all, aes(x=mean, y=factor(trt, levels=trt.names),
                           fill=factor(spp, levels=all.sp))) + 
       geom_bar(stat="identity") +
       facet_wrap(.~factor(type, levels=types), scales="free_x", ncol=3) +
       labs(x="Posterior mean C stock (Mg/ha)", y="", fill="Species") +
       guides(fill=guide_legend(ncol=1)) + 
       theme(text = element_text(linewidth=14))

################################################################################
### run PCA on biomass data
library(ggfortify)
library(ggrepel)  
library(ggalt)
library(plyr)

# call prcomp()
bm.pc = bm.df
omit.vars = c("num",
              "H. jap C stock","P. arun C stock","Mixed biomass C stock",
              "FWD C:N ratio","Litter C:N ratio","Herbaceous C:N ratio")
bm.pc = bm.pc[-which(bm.pc$trt=="R"), -which(colnames(bm.pc) %in% omit.vars)]
pca.biomass = prcomp(bm.pc[,-which(colnames(bm.pc) %in% c("trt","trt.full"))], 
                     center=T, scale=T)
#calculate total variance explained by each principal component
var_explained = pca.biomass$sdev^2 / sum(pca.biomass$sdev^2)

#create scree plot
plot(seq(1:length(var_explained)), var_explained, type="b")

# get loadings matrix
CAloadings <- data.frame(Variables = rownames(pca.biomass$rotation), pca.biomass$rotation)

# add PCA scores to the dataset
bm.pc[, c('PC1', 'PC2')] = pca.biomass$x[, 1:2]

# save variable loadings in a separate dataset
rot = as.data.frame(pca.biomass$rotation[, 1:2])
rot$var = rownames(pca.biomass$rotation)

# rescale the loadings to fit nicely within the scatterplot of our data
mult = max(abs(bm.pc[, c('PC1', 'PC2')])) / max(abs(rot[, 1:2])) / 2
rot[, 1:2] = rot[, 1:2] * mult

# if there are many variables to plot, you can play with ggrepel 
var1 = paste("PC1", paste("(", 100*round(var_explained[1], 4), "%)", sep=""), sep=" ")
var2 = paste("PC2", paste("(", 100*round(var_explained[2], 4), "%)", sep=""), sep=" ")
ggplot(data = rot, aes(x=0, y=0, xend=PC1, yend=PC2, label=var)) +
  geom_point(data = bm.pc, aes(x=PC1, y=PC2, color=factor(trt.full, levels=trt.names[1:5])), inherit.aes=FALSE, linewidth=4) +
  geom_segment(color = 'red', arrow = arrow(length = unit(0.03, "npc"))) +
  geom_text_repel(aes(PC1 * 1, PC2 * 1),color = 'red') +
  labs(x = var1, y = var2, color="Treatment") + 
  theme(legend.title = element_text(linewidth=14),
        legend.text = element_text(linewidth=12), 
        axis.title = element_text(linewidth=14), 
        axis.text = element_text(linewidth=10))
barplot((pca.biomass$sdev^2) / sum(pca.biomass$sdev^2))

################################################################################
# linear analysis code

# fit model with reference for each y var
#lm.mods = list()
#for (i in 1:(n.v-2)) {
#  v = c.vars[i]
#  m.v = ulam(alist(y ~ normal(mu, sigma),
#                   log(mu) <- a + b*N,
#                   a ~ dnorm(0, 1),
#                   b ~ dnorm(0, 1),
#                   sigma ~ dexp(1)),
#             data=c.lists[[v]], chains=10, log_lik=T)
#  lm.mods[[v]] = m.v
#}

# make plot of effects
c.var.labs2 = c("Total ecosystem","Bulk soil","Soil organic carbon","Soil inorganic carbon",
                "Total biomass","Woody biomass","Herbaceous biomass",
                "Herbaceous litter","Fine woody debris")
n.v2 = n.v-2
all.b = data.frame(matrix(nrow=0, ncol=6))
all.b.nr = data.frame(matrix(nrow=0, ncol=6))
colnames(all.b) = c("var","mean","25%","75%","5%","95%")
colnames(all.b.nr) = c("var","mean","25%","75%","5%","95%")
for (i in 1:n.v2) {
  v = c.vars[i]
  b.df = data.frame(matrix(nrow=0, ncol=6))
  b.df.nr = data.frame(matrix(nrow=0, ncol=6))
  colnames(b.df) = c("var","mean","25%","75%","5%","95%")
  colnames(b.df.nr) = c("var","mean","25%","75%","5%","95%")
  b.df[1,c("mean","25%","75%")] = precis(lm.mods[[v]], depth=2, prob=0.50)[2, c("mean","25%","75%")]
  b.df.nr[1,c("mean","25%","75%")] = precis(lm.mods.nr[[v]], depth=2, prob=0.50)[2, c("mean","25%","75%")]
  b.df[1,c("5%","95%")] = precis(lm.mods[[v]], depth=2, prob=0.90)[2, c("5%","95%")]
  b.df.nr[1,c("5%","95%")] = precis(lm.mods.nr[[v]], depth=2, prob=0.90)[2, c("5%","95%")]
  b.df[1,"var"] = c.var.labs2[i]
  b.df.nr[1,"var"] = c.var.labs2[i]
  all.b = rbind(all.b, b.df)
  all.b.nr = rbind(all.b.nr, b.df.nr)
}
all.b$Reference = "Included"
all.b.nr$Reference = "Excluded"
b.r.nr = rbind(all.b, all.b.nr)
b.r.nr$Reference = factor(b.r.nr$Reference)
ggplot(b.r.nr, aes(y=factor(var, levels=c.var.labs2), x=mean, shape=factor(Reference))) + 
  geom_vline(xintercept = 0, linetype=1, color = "darkgray", linewidth=0.5) + 
  geom_errorbar(aes(xmin=`5%`,xmax=`95%`), width=.1, color="black", position=position_dodge(width=0.5)) +
  geom_errorbar(aes(xmin=`25%`,xmax=`75%`), width=.1, color="blue", position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width = 0.5), linewidth=1.2) +
  labs(x="Effect linewidth", y="") +
  theme(legend.title = element_blank()) + 
  xlim(-2.7,1.6) 

# try plotting distributions
post.b = data.frame(matrix(nrow=0, ncol=2))
post.b.nr = data.frame(matrix(nrow=0, ncol=2))
colnames(post.b) = c("var","value")
colnames(post.b.nr) = c("var","value")
for (i in 1:n.v2) {
  v = c.vars[i]
  b.samples = extract.samples(lm.mods[[v]])$b
  b.samples.nr = extract.samples(lm.mods.nr[[v]])$b
  b.df = data.frame(matrix(nrow=length(b.samples), ncol=2))
  b.df.nr = data.frame(matrix(nrow=length(b.samples.nr), ncol=2))
  colnames(b.df) = c("var","value")
  colnames(b.df.nr) = c("var","value")
  b.df$var = v
  b.df.nr$var = v
  b.df$value = as.numeric(b.samples)
  b.df.nr$value = as.numeric(b.samples.nr)
  post.b = rbind(post.b, b.df)
  post.b.nr = rbind(post.b.nr, b.df.nr)
}
post.b$Reference = "Included"
post.b.nr$Reference = "Excluded"
post.b.all = rbind(post.b, post.b.nr)

ggplot(post.b.all, aes(x=value, fill=factor(Reference))) + 
  geom_density() + facet_wrap(.~var, ncol=3, scales="free_y") +
  geom_vline(xintercept = 0, linetype=1, color = "darkgray", linewidth=0.5)

# plot posterior intervals for C stocks v. species richness
n = 500
df.all = data.frame(matrix(nrow=0, ncol=7))
df.all.nr = data.frame(matrix(nrow=0, ncol=7))
colnames(df.all) = c("var","N","mu","lb.5","ub.95","lb.25","ub.75")
colnames(df.all.nr) = c("var","N","mu","lb.5","ub.95","lb.25","ub.75")
for (i in 1:n.v2) {
  v = c.vars[i]
  N.seq = seq(0.34, 1.65, length.out=n)
  
  df.v = data.frame(matrix(nrow=n, ncol=7))
  df.v.nr = data.frame(matrix(nrow=n, ncol=7))
  colnames(df.v) = c("var","N","mu","lb.5","ub.95","lb.25","ub.75")
  colnames(df.v.nr) = c("var","N","mu","lb.5","ub.95","lb.25","ub.75")
  
  df.v$N = N.seq*var.means[["N.herb"]]
  df.v.nr$N = N.seq*var.means.nr[["N.herb"]]
  df.v$var = c.var.labs2[i]
  df.v.nr$var = c.var.labs2[i]
  
  l.v <- link( lm.mods[[v]] , data=list( N=N.seq ) )
  l.v.nr <- link( lm.mods.nr[[v]] , data=list( N=N.seq ) )
  
  mu.v = exp(apply( l.v , 2 , mean ))*var.means[[v]]
  mu.v.nr = exp(apply( l.v.nr , 2 , mean ))*var.means.nr[[v]]
  ci.50 <- apply( l.v , 2 , HPDI, prob=0.50)
  ci.50.nr <- apply( l.v.nr , 2 , HPDI, prob=0.50)
  ci.90 <- apply( l.v , 2 , HPDI, prob=0.90)
  ci.90.nr <- apply( l.v.nr , 2 , HPDI, prob=0.90)
  
  df.v$mu = mu.v
  df.v$lb.5 = exp(ci.90[1,])*var.means[[v]]
  df.v$ub.95 = exp(ci.90[2,])*var.means[[v]]
  df.v$lb.25 = exp(ci.50[1,])*var.means[[v]]
  df.v$ub.75 = exp(ci.50[2,])*var.means[[v]]
  
  df.v.nr$mu = mu.v.nr
  df.v.nr$lb.5 = exp(ci.90.nr[1,])*var.means.nr[[v]]
  df.v.nr$ub.95 = exp(ci.90.nr[2,])*var.means.nr[[v]]
  df.v.nr$lb.25 = exp(ci.50.nr[1,])*var.means.nr[[v]]
  df.v.nr$ub.75 = exp(ci.50.nr[2,])*var.means.nr[[v]]
  
  df.all = rbind(df.all, df.v)
  df.all.nr = rbind(df.all.nr, df.v.nr)
}
df.all$Reference = "Included"
df.all.nr$Reference = "Excluded"
df.r.nr = rbind(df.all, df.all.nr)
df.r.nr$Reference = factor(df.r.nr$Reference)

ggplot(df.r.nr, aes(x=N, y=mu, color=Reference)) + 
       geom_line() +
       geom_ribbon(aes(ymin=lb.5,ymax=ub.95),alpha=0.1,linewidth=0) +
       #geom_ribbon(aes(ymin=lb.25,ymax=ub.75),alpha=0.1,linewidth=0,fill="blue") +
       facet_wrap(.~factor(var, levels=c.var.labs2), scales="free_y") + xlim(1.8,8.2) + 
       labs(y="C stock (Mg/ha)",x="Understory Species Richness")

par(mfrow=c(1,1))
df.all$var = factor(df.all$var, levels=c.var.labs)
bm.plot = bm.df[,c("trt","N",c.vars)]
colnames(bm.plot)[3:10] = c.var.labs
bm.plot$N = factor(bm.plot$N)
bm.melt = melt(bm.plot, id.variables=c("trt","N"))
ggplot(df.all, aes(x=N, y=mu)) + geom_line() +
       geom_ribbon(aes(ymin=lb.5,ymax=ub.95),alpha=0.1,linewidth=0) +
       geom_ribbon(aes(ymin=lb.25,ymax=ub.75),alpha=0.1,linewidth=0,fill="blue") +
       facet_wrap(.~factor(var, levels=c.var.labs), scales="free_y") + 
       labs(y="C stock (Mg/ha)",x="Understory Species Richness")
