path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch2_Floodplain_Experiment/Floodplain-Experiment-Repo"
setwd(path_to_repo)

library(tidyr)
library(dplyr)
library(tidyverse)
library(reshape2)

### scripts that combines C/N concentration data with dry masses to estimate
### C/N stocks of fine woody debris, herbaceous litter, and herbaceous biomass
### in aggregate and by species group

################################################################################
## step 1: clean C/N concentration data

# read in data
cn.data = read.csv("Understory_Analysis/Raw_Data/Vegetation_CN_Concentrations_060523.csv")
colnames(cn.data) = c("well","id","n.percent","c.percent")

# separate ID column
cn.data = cn.data %>% separate(id, c("treatment_plot","type"), remove=T)
cn.data = cn.data %>% separate(treatment_plot, into = c('treatment','plot'), sep=1, remove=T)

# calculate CN ratio
cn.data$cn.ratio = as.numeric(cn.data$c.percent)/as.numeric(cn.data$n.percent)

# add column for full treatment name
trt.letters = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")
n.t = length(trt.letters)
cn.data$full.treatment.name = rep(0, nrow(cn.data))
for (i in 1:n.t) { cn.data$full.treatment.name[which(cn.data$treatment == trt.letters[i])] = trt.names[i] }

# types to match stock data
type.abvs = c("FWD","HL","HJ","PA","mix")
type.long = c("Fine.woody.debris","Herbaceous.litter","Humulus.japonicus","Phalaris.arundinacea","Mixed.species")
n.types = length(type.abvs)
for (i in 1:n.types) { cn.data$type[which(cn.data$type == type.abvs[i])] = type.long[i] }

# separate CN data by type
cn.data$plot = as.numeric(cn.data$plot)
fwd.cn.data = cn.data[which(cn.data$type == "Fine.woody.debris"),]
hl.cn.data = cn.data[which(cn.data$type == "Herbaceous.litter"),]
bm.cn.data = cn.data[-which(cn.data$type %in% c("Fine.woody.debris","Herbaceous.litter")),]

################################################################################
## step 2: average biomass data by plot

# read in biomass data 
bm.df = read.csv("Understory_Analysis/Clean_Data/Biomass_By_Quadrat_Sep2022.csv")

# melt data
bm.melt = melt(bm.df, id.vars=c("treatment","full.treatment.name","plot","quadrat"),
               value.name="mass", variable.name="type")

# average litter, biomass, and fwd data over quadrats
bm.ave = bm.melt %>% 
         group_by(treatment, full.treatment.name, plot, type) %>% 
         summarise(ave.mass = mean(mass))

# estimate C and N stocks for fwd
fwd.mass.data = bm.ave[which(bm.ave$type == "Fine.woody.debris"),]
fwd.stocks = left_join(fwd.mass.data, fwd.cn.data, by=c("type","treatment","full.treatment.name","plot"))
fwd.stocks$c.g.m2 = fwd.stocks$ave.mass * fwd.stocks$c.percent / 100
fwd.stocks$n.g.m2 = fwd.stocks$ave.mass * fwd.stocks$n.percent / 100

# estimate C and N stocks for hl
hl.mass.data = bm.ave[which(bm.ave$type == "Herbaceous.litter"),]
hl.stocks = left_join(hl.mass.data, hl.cn.data, by=c("type","treatment","full.treatment.name","plot"))
hl.stocks$c.g.m2 = hl.stocks$ave.mass * hl.stocks$c.percent / 100
hl.stocks$n.g.m2 = hl.stocks$ave.mass * hl.stocks$n.percent / 100

# estimate C and N stocks for biomass by type
bm.mass.data = bm.ave[-which(bm.ave$type %in% c("Fine.woody.debris","Herbaceous.litter","Herbaceous.biomass")),]
bm.stocks = left_join(bm.mass.data, bm.cn.data, by=c("type","treatment","full.treatment.name","plot"))
na.id = which(is.na(bm.stocks$c.percent))
bm.stocks$c.g.m2 = rep(0, nrow(bm.stocks))
bm.stocks$n.g.m2 = rep(0, nrow(bm.stocks))
bm.stocks$c.g.m2[-na.id] = bm.stocks$ave.mass[-na.id] * bm.stocks$c.percent[-na.id] / 100
bm.stocks$n.g.m2[-na.id] = bm.stocks$ave.mass[-na.id] * bm.stocks$n.percent[-na.id] / 100

# fill in zeros for total herbaceous biomass
tot.bm.mass.data = bm.ave[which(bm.ave$type %in% c("Herbaceous.biomass")),]
tot.bm.stock.data = bm.stocks[,c("type","treatment","full.treatment.name","plot","c.g.m2","n.g.m2")] %>% 
                    group_by(treatment, full.treatment.name, plot) %>% 
                    summarise(c.g.m2 = sum(c.g.m2),
                              n.g.m2 = sum(n.g.m2))
tot.bm.stocks = left_join(tot.bm.mass.data, tot.bm.stock.data, by=c("treatment","full.treatment.name","plot"))

# write all C and N stock data to csv files
cols = c("treatment","full.treatment.name","plot","type","ave.mass","c.g.m2","n.g.m2")
all.stock.data = rbind(tot.bm.stocks[,cols], bm.stocks[,cols], hl.stocks[,cols], fwd.stocks[,cols])
write.csv(all.stock.data, "Understory_Analysis/Clean_Data/CN_Stock_Summary_Sep2022.csv", row.names=F)


