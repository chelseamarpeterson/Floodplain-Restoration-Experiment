path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo"
setwd(path_to_repo)

library(tidyr)
library(dplyr)
library(patchwork)
library(tidyverse)
library(reshape2)

################################################################################
## step 1: clean C & N concentration data

# read in data
cn.data = read.csv("Understory_Analysis/Raw_Data/Vegetation_CN_Concentrations_060523.csv")
colnames(cn.data) = c("well","id","N.percent","C.percent")

# separate Id column
cn.data = cn.data %>% separate(id, c("plot","type"), remove=T)
cn.data = cn.data %>% separate(plot, into = c('trt','num'), sep = 1, remove=T)

# calculate CN ratio
cn.data$CN = as.numeric(cn.data$C)/as.numeric(cn.data$N)

# add column for full treatment name
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")
cn.data$trt.full = rep(0, nrow(cn.data))
for (i in 1:6) { cn.data$trt.full[which(cn.data$trt == trts[i])] = trt.names[i] }

# types to match stock data
type.abvs = c("FWD","HL","HJ","PA","mix")
type.long = c("Fine.woody.debris","Herbaceous.litter","Humulus.japonicus","Phalaris.arundinacea","Mixed.species")
for (i in 1:5) { cn.data$type[which(cn.data$type == type.abvs[i])] = type.long[i] }

# separate CN data by type
cn.data$num = as.numeric(cn.data$num)
fwd.cn.data = cn.data[which(cn.data$type == "Fine.woody.debris"),]
hl.cn.data = cn.data[which(cn.data$type == "Herbaceous.litter"),]
bm.cn.data = cn.data[-which(cn.data$type %in% c("Fine.woody.debris","Herbaceous.litter")),]

################################################################################
## step 2: average organize biomass data by plot

# read in biomass data 
bm.all = read.csv("Understory_Analysis/Clean_Data/Biomass_By_Quadrat_and_Year.csv")

# isolate 2022 data
bm.2022 = bm.all[which(bm.all$year == 2022),]

# melt 2022 data
bm.melt = melt(bm.2022, id.vars=c("year","trt","trt.full","num","quad"),
               value.name="mass", variable.name="type")

# average litter, biomass, and fwd data over quadrats
bm.ave = bm.melt %>% 
         group_by(year, trt, trt.full, num, type) %>% 
         summarise(ave.mass = mean(mass))

# estimate C and N stocks for fwd
fwd.mass.data = bm.ave[which(bm.ave$type == "Fine.woody.debris"),]
fwd.stocks = left_join(fwd.mass.data, fwd.cn.data, by=c("type","trt","trt.full","num"))
fwd.stocks$C.g.m2 = fwd.stocks$ave.mass * fwd.stocks$C.percent / 100
fwd.stocks$N.g.m2 = fwd.stocks$ave.mass * fwd.stocks$N.percent / 100

# estimate C and N stocks for hl
hl.mass.data = bm.ave[which(bm.ave$type == "Herbaceous.litter"),]
hl.stocks = left_join(hl.mass.data, hl.cn.data, by=c("type","trt","trt.full","num"))
hl.stocks$C.g.m2 = hl.stocks$ave.mass * hl.stocks$C.percent / 100
hl.stocks$N.g.m2 = hl.stocks$ave.mass * hl.stocks$N.percent / 100

# estimate C and N stocks for biomass by type
bm.mass.data = bm.ave[-which(bm.ave$type %in% c("Fine.woody.debris","Herbaceous.litter","Herbaceous.biomass")),]
bm.stocks = left_join(bm.mass.data, bm.cn.data, by=c("type","trt","trt.full","num"))
na.id = which(is.na(bm.stocks$C.percent))
bm.stocks$C.g.m2 = rep(0, nrow(bm.stocks))
bm.stocks$N.g.m2 = rep(0, nrow(bm.stocks))
bm.stocks$C.g.m2[-na.id] = bm.stocks$ave.mass[-na.id] * bm.stocks$C.percent[-na.id] / 100
bm.stocks$N.g.m2[-na.id] = bm.stocks$ave.mass[-na.id] * bm.stocks$N.percent[-na.id] / 100

# fill in zeros for total herbaceous biomass
tot.bm.mass.data = bm.ave[which(bm.ave$type %in% c("Herbaceous.biomass")),]
tot.bm.stock.data = bm.stocks[,c("trt","trt.full","num","type","C.g.m2","N.g.m2")] %>% 
                    group_by(trt, trt.full, num) %>% 
                    summarise(C.g.m2 = sum(C.g.m2),
                              N.g.m2 = sum(N.g.m2))
tot.bm.stocks = left_join(tot.bm.mass.data, tot.bm.stock.data, by=c("trt","trt.full","num"))

# write all C and N stock data to csv files
cols = c("trt","trt.full","num","type","ave.mass","C.g.m2","N.g.m2")
all.stock.data = rbind(tot.bm.stocks[,cols], bm.stocks[,cols], hl.stocks[,cols], fwd.stocks[,cols])
write.csv(all.stock.data, "Understory_Analysis/Clean_Data/CN_Stock_Summary_2022.csv", row.names=F)


