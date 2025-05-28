path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Public-Repo"
setwd(path_to_repo)

library(ggplot2)
library(readxl)
library(ggfortify)
library(patchwork)
library(reshape2)
library(dplyr)
library(tidyr)

# treatments
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

## read in biomass data for 2022 and 2023
bm.2023 = read_excel("Understory_Analysis/Raw_Data/Joslin_Vegetation_Cover_Biomass_Aug2023.xlsx", sheet="Biomass")
bm.2022 = read_excel("Understory_Analysis/Raw_Data/HerbaceousBiomass_Sep2022.xlsx", sheet="ActualBiomass")

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

# normalize biomass data by quadrat size (g/m^2)
bm.all$mass.g.m2 = bm.all$sum.mass/(30^2)*(100^2)

# separate plot column into treatment and plot number 
bm.all = bm.all %>% separate(plot, into = c("trt", "num"), sep = 1, remove=F)

# isolate biomass data
bm.dat = bm.all[!(bm.all$spp == "Herbaceous litter" | bm.all$spp == "Fine woody debris"),]

# add column for full treatment name
bm.dat$trt.full = rep(0, nrow(bm.dat))
for (i in 1:6) {bm.dat$trt.full[which(bm.dat$trt == trts[i])] = trt.names[i]}
bm.sp.sum = bm.dat %>% group_by(year, plot, trt, trt.full, spp) %>% summarise(mean.mass = sum(mass.g.m2)/3)

# make data frame for total herbaceous biomass, herbaceous litter, and fine woody debris in each quadrat 
n.y = 2; n.t = 6; n.p = 3; n.q = 5
df.quad.sum = data.frame(matrix(nrow = n.y*n.t*n.p*n.q, ncol=11))
colnames(df.quad.sum) = c("year","trt","trt.full","num","quad","Herbaceous biomass","Herbaceous litter","Fine woody debris","Mixed species","Phalaris arundinacea","Humulus japonicus")
k = 1
for (y in c(2022, 2023)) {
  for (l in 1:n.t) {
    t = trts[l]
    for (i in 1:n.p) {
      for (j in 1:n.q) {
        # get data for given year, treatment, plot, and quadrat
        quad.id = which((bm.all$year == y & bm.all$trt == t) & (bm.all$num == as.character(i) & bm.all$quad == paste("Q", j, sep="")))
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
        df.quad.sum[k, c("year","trt","trt.full","num","quad")] = c(y, t, trt.names[l], i, paste("Q", j, sep=""))
        
        # increment counter
        k = k + 1
      }
    }
  }
}

# write data by quadrat to file
write.csv(df.quad.sum, "Understory_Analysis/Clean_Data/Biomass_By_Quadrat_and_Year.csv", row.names=F)


# melt biomass dataframe
df.melt = melt(df.quad.sum, id.vars = c("year","trt","trt.full","num","quad"), variable.name = "type", value.name="mass.g.m2")

# dataframes for herb biomass, herb lit, and fwd in each plot summed over quadrats
bm.plt.total = data.frame(matrix(nrow=n.t*n.p, ncol=3))
hl.plt.total = data.frame(matrix(nrow=n.t*n.p, ncol=3))
fwd.plt.total = data.frame(matrix(nrow=n.t*n.p, ncol=3))
colnames(bm.plt.total) = c("trt","plt","mass.g.m2")
colnames(hl.plt.total) = c("trt","plt","mass.g.m2")
colnames(fwd.plt.total) = c("trt","plt","mass.g.m2")
k = 1
for (t in trts) {
  for (i in 1:n.p) {
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

# sort biomass dataframe
bm.sort = bm.plt.sum[order(bm.plt.sum$plot,-bm.plt.sum$mass.g),]

# data frame of most abundant species in each plot
df.sp.sort = data.frame(matrix(nrow=n.t*n.p, ncol=max.n.sp+2))
colnames(df.sp.sort) = c("trt","plot",seq(1,max.n.sp))
k = 1
for (t in trts) {
  for (i in 1:n.p) {
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
  for (i in 1:n.p) {
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


# sum up litter across quadrats by species
lit.plt.sum = data.frame(matrix(nrow=n.plots, ncol=2))
colnames(lit.plt.sum) = c("trt","mass.g")
k = 1
for (t in trts) {
  for (i in 1:n.p) {
    plt.ind = which(lit.dat$plot == paste(t, i, sep=""))
    lit.plt.sum[k,"trt"] = t
    lit.plt.sum[k,"mass.g"] = sum(lit.dat$dry.mass[plt.ind])
    k = k + 1
  }
}
lit.plt.sum$mass.g.m2 = lit.plt.sum$mass.g/5/(30^2)*(100^2)
