setwd('C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/VegetationData')

library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stats)
library(tidyverse)
library(ggthemes)
library(multcompView)

################################################################################
## clean and analyze C & N concentration data

# read in data
cn.data = read_excel("Vegetation_CHN_Results_060523.xlsx")
colnames(cn.data) = c("Well","ID","N","C")

# separate Id column
cn.data = cn.data %>% separate(ID, c("plot", "type"), remove=T)
cn.data = cn.data %>% separate(plot, into = c('trt', 'num'), sep = 1, remove=T)

# calculate CN ratio
cn.data$CN = as.numeric(cn.data$C)/as.numeric(cn.data$N)

# add column for full treatment name
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")
cn.data$trt.full = rep(0, nrow(cn.data))
for (i in 1:6) { cn.data$trt.full[which(cn.data$trt == trts[i])] = trt.names[i] }

# separate litter and biomass data
lit.cn.data = cn.data[which(cn.data$type %in% c("FWD","HL")),]
fwd.cn.data = cn.data[which(cn.data$type %in% c("FWD")),]
hl.cn.data = cn.data[which(cn.data$type %in% c("HL")),]
bm.cn.data = cn.data[which(cn.data$type %in% c("HJ","PA","mix")),]

### all plots 

# plot litter data by treatment
p1 = ggplot(lit.cn.data, aes(x=factor(trt.full,levels=trt.names), y=C, color=type)) + 
            geom_boxplot() + 
            labs(y="%C") +
            theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p2 = ggplot(lit.cn.data, aes(x=factor(trt.full,levels=trt.names), y=N, color=type)) + 
            geom_boxplot() + 
            labs(y="%N") + 
            theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p3 = ggplot(lit.cn.data, aes(x=factor(trt.full,levels=trt.names), y=CN, color=type)) + 
            geom_boxplot() + 
            labs(y="C/N") +
            theme(axis.title.x=element_blank()) +  
            ylim(0,70)
p1/p2/p3

# plot litter data by type
p1 = ggplot(lit.cn.data, aes(x=type, y=C)) + 
            geom_boxplot() + 
            labs(y="%C") +
            theme(axis.title.x=element_blank())
p2 = ggplot(lit.cn.data, aes(x=type, y=N)) + 
            geom_boxplot() + 
            labs(y="%N") +
            theme(axis.title.x=element_blank())
p3 = ggplot(lit.cn.data, aes(x=type, y=CN)) + 
            geom_boxplot() + 
            labs(y="C/N") +
            theme(axis.title.x=element_blank())
p1+p2+p3

################################################################################
## estimate carbon stocks 

# read in biomass data 
bm.all = read_excel('HerbaceousBiomass_Sep2022.xlsx', sheet="ActualBiomass")
colnames(bm.all) = c("plot","quad","spp","dry.mass.w.container","container.mass","dry.mass","diameter")
bm.all = bm.all[,-which(colnames(bm.all) %in% c("dry.mass.w.container","container.mass","diameter"))]

# split plot column
bm.all = bm.all %>% separate(plot, into = c('trt', 'num'), sep = 1)

# add column for full treatment name
bm.all$trt.full = rep(0, nrow(bm.all))
for (i in 1:6) { bm.all$trt.full[which(bm.all$trt == trts[i])] = trt.names[i] }

# separate biomass & litter data
bm.data = bm.all[!(bm.all$spp == "Herbaceous litter" | bm.all$spp == "Fine woody debris"),]
hl.data = bm.all[bm.all$spp == "Herbaceous litter",]
fwd.data = bm.all[bm.all$spp == "Fine woody debris",]

# sum hl & fwd data over quadrats
hl.data$num = as.factor(hl.data$num)
hl.sum = hl.data %>% 
         group_by(trt, trt.full, num) %>% 
         summarise(sum.mass = sum(dry.mass))
fwd.sum = fwd.data %>% 
          group_by(trt, trt.full, num) %>% 
          summarise(sum.mass = sum(dry.mass))
fwd.sum = rbind(fwd.sum, data.frame(trt="E", trt.full="Seedbank", num="1", sum.mass=0))

# sum biomass data over quadrats by type
bm.cn.data$sum.mass = rep(0, nrow(bm.cn.data))
for (i in 1:6) {
  for (j in 1:3) {
    # find row for specific treatment and plot in bm.data
    t = trts[i]
    bm.trt.plot = bm.data[which(bm.data$trt == t & bm.data$num ==j),]
    
    # find rows for specific treatment and plot in 
    bm.pa.id = which((bm.cn.data$trt == t & bm.cn.data$num == j) & bm.cn.data$type == "PA")
    bm.hj.id = which((bm.cn.data$trt == t & bm.cn.data$num == j) & bm.cn.data$type == "HJ")
    bm.mix.id = which((bm.cn.data$trt == t & bm.cn.data$num == j) & bm.cn.data$type == "mix")
    
    # get masses of phalaris and humulus
    bm.cn.data[bm.pa.id,"sum.mass"] = sum(bm.trt.plot[which(bm.trt.plot$spp == "Phalaris arundinacea"),"dry.mass"])
    bm.cn.data[bm.hj.id,"sum.mass"] = sum(bm.trt.plot[which(bm.trt.plot$spp == "Humulus japonicus"),"dry.mass"])
    
    # sum mass of all species except phalaris and humulus
    bm.pa.hj.id = which(bm.trt.plot$spp %in% c("Phalaris arundinacea","Humulus japonicus"))
    if (length(bm.pa.hj.id) > 0) {
      bm.cn.data[bm.mix.id,"sum.mass"] = sum(bm.trt.plot[-bm.pa.hj.id,"dry.mass"])
    } else {
      bm.cn.data[bm.mix.id,"sum.mass"] = sum(bm.trt.plot[,"dry.mass"])
    }
  }
}

# divide sums by quadrat area
hl.sum$mass.g.m2 = hl.sum$sum.mass/5/(30^2)*(100^2)
fwd.sum$mass.g.m2 = fwd.sum$sum.mass/5/(30^2)*(100^2)
bm.cn.data$mass.g.m2 = bm.cn.data$sum.mass/5/(30^2)*(100^2)

# estimate C and N stocks for fwd
fwd.sum$C.g.m2 = rep(0, nrow(fwd.sum))
fwd.sum$N.g.m2 = rep(0, nrow(fwd.sum))
for (i in 1:nrow(fwd.sum)) {
  t = fwd.sum$trt[i]
  p = fwd.sum$num[i]
  percent.C = as.numeric(fwd.cn.data[which(fwd.cn.data$trt == t & fwd.cn.data$num == j), "C"])
  percent.N = as.numeric(fwd.cn.data[which(fwd.cn.data$trt == t & fwd.cn.data$num == j), "N"])
  fwd.sum$C.g.m2[i] = fwd.sum$mass.g.m2[i] * percent.C/100
  fwd.sum$N.g.m2[i] = fwd.sum$mass.g.m2[i] * percent.N/100
}

# estimate C and N stocks for hl
hl.sum$C.g.m2 = rep(0, nrow(hl.sum))
hl.sum$N.g.m2 = rep(0, nrow(hl.sum))
for (i in 1:nrow(hl.sum)) {
  t = hl.sum$trt[i]
  p = hl.sum$num[i]
  percent.C = as.numeric(hl.cn.data[which(hl.cn.data$trt == t & hl.cn.data$num == j), "C"])
  percent.N = as.numeric(hl.cn.data[which(hl.cn.data$trt == t & hl.cn.data$num == j), "N"])
  hl.sum$C.g.m2[i] = hl.sum$mass.g.m2[i] * percent.C/100
  hl.sum$N.g.m2[i] = hl.sum$mass.g.m2[i] * percent.N/100
}

hl.sum$type = "HL"
fwd.sum$type = "FWD"
lit.sum = rbind(hl.sum, fwd.sum)

# estimate C and N stocks for mixed biomass
bm.cn.data$C.g.m2 = bm.cn.data$mass.g.m2 * bm.cn.data$C / 100
bm.cn.data$N.g.m2 = bm.cn.data$mass.g.m2 * bm.cn.data$N / 100

# fill in zeros for biomass data
bm.cn.fill = data.frame(matrix(nrow=0, ncol=ncol(bm.cn.data)))
colnames(bm.cn.fill) = colnames(bm.cn.data)
types = c("mix","HJ","PA")
for (type in types) {
  for (i in 1:6) {
    trt = trts[i]
    for (j in 1:3) {
      trt.plot.id = which(bm.cn.data$type == type & (bm.cn.data$trt == trt & bm.cn.data$num == j))
      if (length(trt.plot.id) > 0) {
        bm.cn.fill = rbind(bm.cn.fill, bm.cn.data[trt.plot.id,])
      } else {
        bm.cn.new = data.frame(matrix(nrow=1, ncol=ncol(bm.cn.data)))
        colnames(bm.cn.new) = colnames(bm.cn.data)
        bm.cn.new[1, c("trt","num","type", "trt.full")] = c(trt, j, type, trt.names[i])
        bm.cn.new[1, c("sum.mass","mass.g.m2","C.g.m2","N.g.m2")] = rep(0, 4)
        bm.cn.fill = rbind(bm.cn.fill, bm.cn.new)
      }
    }
  }
}


# write all C and N stock data to csv files
all.lit.data = right_join(lit.cn.data, lit.sum, by=c("trt","trt.full","num","type"))
all.biomass.data = rbind(all.lit.data, bm.cn.fill)
all.biomass.data = all.biomass.data[,-which(colnames(all.biomass.data) == "Well")]
write.csv(all.biomass.data, "CN_Summary_2022.csv", row.names=F)

## C & N stock plots

# plot litter C and N stocks by treatment
p1 = ggplot(lit.sum, aes(x=factor(trt.full, levels=trt.names), y=mass.g.m2, color=type)) + 
  geom_boxplot() + 
  labs(y="Biomass (g/m2)") + 
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p2 = ggplot(lit.sum, aes(x=factor(trt.full, levels=trt.names), y=C.g.m2, color=type)) + 
  geom_boxplot() + 
  labs(y="C stock (g/m2)") + 
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p3 = ggplot(lit.sum, aes(x=factor(trt.full, levels=trt.names), y=N.g.m2, color=type)) + 
  geom_boxplot() + 
  labs(y="N stock (g/m2)") +
  theme(axis.title.x=element_blank())
p1/p2/p3

# plot biomass C and N stocks by trt and type
p1 = ggplot(bm.cn.fill, aes(x=factor(trt.full, levels=trt.names), y=mass.g.m2, color=type)) + 
  geom_boxplot() + 
  labs(y="Biomass (g/m2)") + 
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p2 = ggplot(bm.cn.fill, aes(x=factor(trt.full, levels=trt.names), y=C.g.m2, color=type)) + 
  geom_boxplot() + 
  labs(y="C stock (g/m2)") +
  theme(axis.text.x=element_blank(), axis.title.x=element_blank())
p3 = ggplot(bm.cn.fill, aes(x=factor(trt.full, levels=trt.names), y=N.g.m2, color=type)) + 
  geom_boxplot() + 
  labs(y="N stock (g/m2)") +
  theme(axis.title.x=element_blank())
p1/p2/p3

# bar plots
bm.cn.stat = bm.cn.fill %>% 
             group_by(Trt, trt.full, Type) %>% 
             summarize(mean.mass = mean(mass.g.m2), mean.C = mean(C.g.m2), mean.N = mean(N.g.m2), 
                       sd.mass = sd(mass.g.m2), sd.C = sd(C.g.m2), sd.N = sd(N.g.m2), .groups="drop")
p1 = ggplot(bm.cn.stat, aes(x=factor(trt.full, levels=trt.names), y=mean.mass, fill=Type)) + 
            geom_bar(position = position_dodge(), stat="identity") + 
            theme(axis.text.x=element_blank(), axis.title.x=element_blank()) + labs(y="Biomass (g/m2)") +
            geom_errorbar(aes(ymin=mean.mass-sd.mass, ymax=mean.mass+sd.mass), width=0.3, position=position_dodge(0.9))
p2 = ggplot(bm.cn.stat, aes(x=factor(trt.full, levels=trt.names), y=mean.C, fill=Type)) + 
            geom_bar(position = position_dodge(), stat="identity") + 
            theme(axis.text.x=element_blank(), axis.title.x=element_blank()) + labs(y="C stock (g/m2)") +
            geom_errorbar(aes(ymin = mean.C-sd.C, ymax = mean.C+sd.C), width = 0.3, position=position_dodge(0.9))
p3 = ggplot(bm.cn.stat, aes(x=factor(trt.full, levels=trt.names), y=mean.N, fill=Type)) + 
            geom_bar(position=position_dodge(), stat="identity") + 
            theme(axis.title.x=element_blank()) + labs(y="N stock (g/m2)") +
            geom_errorbar(aes(ymin = mean.N-sd.N, ymax = mean.N+sd.N), width = 0.3, position=position_dodge(0.9))
p1/p2/p3

## estimate total C & N stocks in biomass
bm.cn.sum = bm.cn.data %>% group_by(Trt, trt.full, Num) %>% summarize(tot.mass = sum(mass.g.m2), 
                                                                      tot.C = sum(C.g.m2),
                                                                      tot.N = sum(N.g.m2), .groups="drop")
