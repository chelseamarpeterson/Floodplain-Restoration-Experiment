setwd('C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1')

library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(readxl)
library(reshape2)

### script that brings together all biomass carbon stock and plant species richness data
### and writes the results to a CSV file

# treatment labels
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

# list of species families
genus.families = list("Acer" = "Sapindaceae",
                      "Betula" = "Betulaceae",
                      "Carya" = "Juglandaceae",
                      "Fraxinus" = "Oleaceae",
                      "Platanus" = "Platanaceae",
                      "Quercus" = "Fagaceae",
                      "Salix" = "Salicaceae",
                      "Populus" = "Salicaceae",
                      "Ulmus" = "Ulmaceae",
                      "Celtis" = "Cannabaceae",
                      "Morus" = "Moraceae",
                      "Gleditsia" = "Fabaceae")
families = c("Sapindaceae","Betulaceae","Juglandaceae","Oleaceae","Moraceae",
             "Fagaceae","Salicaceae","Ulmaceae","Platanaceae","Cannabaceae",
             "Fabaceae")

# read in woody C data 
wood.c.df = read.csv("C:/Users/Chels/OneDrive - University of Illinois - Urbana/PhD/Databases/Doraisami_2021_Wood_C_Database.csv")

# isolate rows for dead material
dead.wood.c.df = wood.c.df[which(wood.c.df$dead.alive == "dead"),]
dead.wood.c.df = dead.wood.c.df[-which(dead.wood.c.df$tissue %in% c("bark","root","branch")),]
dead.wood.c.df = dead.wood.c.df %>% separate(binomial.resolved, c("genus","spp"), sep="_", remove=F)
dead.wood.c.df$species = paste(dead.wood.c.df$genus, dead.wood.c.df$spp, sep= " ")

# isolate rows for live material
live.wood.c.df = wood.c.df[which(wood.c.df$dead.alive == "alive" & wood.c.df$growth.form == "tree"),]
live.wood.c.df = live.wood.c.df[which(live.wood.c.df$position == "standing"),]
live.wood.c.df = live.wood.c.df[-which(live.wood.c.df$tissue %in% c("bark","branch","coarseroot","fineroot","twig","root","sapwood")),]
live.wood.c.df = live.wood.c.df %>% separate(binomial.resolved, c("genus","spp"), sep="_", remove=F)
live.wood.c.df$species = paste(live.wood.c.df$genus, live.wood.c.df$spp, sep= " ")

### clean coarse and fine woody debris data ---- ###############################
### double-checked #############################################################

# load dataframe and change columne names
cwd.data = read.table("TreeData/CWD_June2023.csv", header=T, sep=",")
colnames(cwd.data) = c("plot","position","diameter","species","fwd.count")

# split plot name column
cwd.data = cwd.data %>% separate(plot, c("trt","num"), sep = 1, remove=F)

# add column for full treatment name
cwd.data$trt.full = rep(0, nrow(cwd.data))
for (i in 1:6) { cwd.data$trt.full[which(cwd.data$trt == trts[i])] = trt.names[i] }

## make data frame for fine woody debris (2.5-7.5 cm) count (#/m)
fwd.counts = cwd.data[which(is.na(cwd.data$diameter)), c("trt","trt.full","num","fwd.count")]
row.names(fwd.counts) = seq(1, nrow(fwd.counts))
trt.sort = sort(fwd.counts$trt, index.return=T)
fwd.counts = fwd.counts[trt.sort$ix,]

# estimate FWD volume (m^3/ha) and C storage (Mg/ha)
C1 = 21.1659
ftPerM = 3.28084
cmPerIn = 2.54
C.content.fwd = 0.484 #% (averaged from measured FWD data)
cwd.density.decay.class3 = 0.26 # g/cm^3
fwd.counts$fwd.vol = (1/(48*ftPerM)) * (C1*(pi^2)/8) * fwd.counts$fwd.count * ((5/cmPerIn)^2) #m3/ha
fwd.counts$fwd.carbon = fwd.counts$fwd.vol * C.content.fwd * cwd.density.decay.class3 # m3/ha * g/cm3 * g/g = Mg/ha
fwd.counts$fwd.count = fwd.counts$fwd.count / 48 # count -> count/m

# #make data frame for coarse woody debris (>7.5 cm)
cwd.diams = cwd.data[-which(is.na(cwd.data$diameter)), -which(colnames(cwd.data) == "fwd.count")]
cwd.diams = cwd.diams %>% separate(species, c("genus","spp"), sep=" ", remove=F)

# calculate CWD basal area (m2/km)
mPerKm = 1000
cwd.diams$cwd.area = (pi * (cwd.diams$diameter/100/2)^2) / 48 * mPerKm

# calculate average downed angiosperm C content
C.contents = c(47.2,47.8,47.7,48.1,47.4,47.3)
n.samples = c(8,21,25,32,17,4)
cwd.C.ave = sum((C.contents*n.samples))/sum(n.samples)

# add column for potential C content by taxon
cwd.C.content.decay.class3 = 0.481 # %
cwd.diams$taxon.c.content = rep(0, nrow(cwd.diams))
cwd.diams$family = rep(NA, nrow(cwd.diams))
for (i in 1:nrow(cwd.diams)) {
  if (cwd.diams$species[i] != "") {
    cwd.diams$family[i] = genus.families[[cwd.diams$genus[i]]]
    df.sp.id = which(dead.wood.c.df$species == cwd.diams$species[i])
    df.gen.id = which(dead.wood.c.df$genus.resolved == cwd.diams$genus[i]) 
    df.fam.id = which(dead.wood.c.df$family.resolved == cwd.diams$family[i]) 
    if (length(df.sp.id) > 0) {
      cwd.diams[i,"taxon.c.content"] = mean(dead.wood.c.df[df.sp.id,"tissue.c"])/100
    } else if (length(df.gen.id > 0)) {
      cwd.diams[i,"taxon.c.content"] = mean(dead.wood.c.df[df.gen.id,"tissue.c"])/100
    } else if (length(df.fam.id > 0)) {
      cwd.diams[i,"taxon.c.content"] = mean(dead.wood.c.df[df.fam.id,"tissue.c"])/100
    } else {
      cwd.diams[i,"taxon.c.content"] = cwd.C.content.decay.class3
    }
  } else {
    cwd.diams[i,"taxon.c.content"] = cwd.C.content.decay.class3
  }
}

# estimate CWD volume (m^3/ha) and carbon storage (Mg/ha)
cwd.diams$cwd.vol = (1/(48*ftPerM)) * (C1*(pi^2)/8) * (cwd.diams$diameter/cmPerIn)^2 #m3/ha
cwd.diams$cwd.carbon = cwd.diams$cwd.vol * cwd.diams$taxon.c.content*cwd.density.decay.class3 # Mg/ha

# sum areas (m2/km), volumes (m^3/ha), and carbon stocks (Mg/ha) by plot
cwd.sum = cwd.diams %>% 
          group_by(trt, trt.full, num) %>% 
          summarize(cwd.area.sum = sum(cwd.area),
                    cwd.vol.sum = sum(cwd.vol),
                    cwd.carbon.sum = sum(cwd.carbon))

# add zeros to CWD areas and volumes
for (i in 1:6) {
  for (j in 1:3) {
    trt.plot.id = which(cwd.sum$trt == trts[i] & cwd.sum$num == j)
    if (length(trt.plot.id) == 0) {
      zero.row = data.frame(matrix(nrow=1, ncol=6))
      colnames(zero.row) = colnames(cwd.sum)
      zero.row[1, c("trt","trt.full","num")] = c(trts[i], trt.names[i], j)
      zero.row[1, c("cwd.area.sum","cwd.vol.sum","cwd.carbon.sum")] = rep(0,3)
      cwd.sum = rbind(cwd.sum, zero.row)
    }
  }
}
trt.sort = sort(cwd.sum$trt, index.return=T)
cwd.sum = cwd.sum[trt.sort$ix,]

# add to biomass dataframe
all.bm.data = right_join(fwd.counts[,c("trt","trt.full","num","fwd.count","fwd.vol","fwd.carbon")], 
                         cwd.sum[,c("trt","trt.full","num","cwd.area.sum","cwd.vol.sum","cwd.carbon.sum")], 
                         by=c("trt","trt.full","num"))

# sum areas (m2/km), volumes (m^3/ha), and carbon stocks (Mg/ha) by species
cwd.diams$species[which(cwd.diams$species == "")] = "Unknown"
cwd.sp.sum = cwd.diams %>% 
             group_by(trt, trt.full, num, species) %>% 
             summarize(cwd.area.sum = sum(cwd.area),
                       cwd.vol.sum = sum(cwd.vol),
                       cwd.carbon.sum = sum(cwd.carbon))

# reshape data frame
cwd.sp.wide = pivot_wider(cwd.sp.sum, 
                          id_cols = c(trt, num),
                          names_from = species, 
                          values_from = cwd.carbon.sum)
cwd.sp.wide[is.na(cwd.sp.wide)] = 0

# add zeros to CWD stocks
cwd.spp = colnames(cwd.sp.wide)[3:6]
for (i in 1:6) {
  for (j in 1:3) {
    trt.plot.id = which(cwd.sp.wide$trt == trts[i] & cwd.sp.wide$num == j)
    if (length(trt.plot.id) == 0) {
      zero.row = data.frame(matrix(ncol=length(colnames(cwd.sp.wide)), nrow=1))
      colnames(zero.row) = colnames(cwd.sp.wide)
      zero.row[1, c("trt","num")] = c(trts[i], j)
      zero.row[1, cwd.spp] = rep(0, length(cwd.spp))
      cwd.sp.wide = rbind(cwd.sp.wide, zero.row)
    }
  }
}
trt.sort = sort(cwd.sp.wide$trt, index.return=T)
cwd.sp.wide = cwd.sp.wide[trt.sort$ix,]

# reshape dataframe
cwd.sp.melt = melt(cwd.sp.wide, 
                   id.vars=c("trt","num"), 
                   variable.name="species",
                   value.name="cwd.carbon.sum")

# write cleaned cwd data to csv
write.csv(cwd.sp.melt, "TreeData/CWD_Carbon_Stocks_By_Species.csv", row.names=F)

### clean snag data ---- #######################################################
### double-checked #############################################################

# load dataframe and change column names
snag.dbh.data = read.table("TreeData/DBH_Snags_June2023.csv", header=T, sep=",")
colnames(snag.dbh.data) = c("redo","plot","live","species","dbh","stem.count","notes")

# split plot name column
snag.dbh.data = snag.dbh.data %>% separate(plot, c("trt","num"), sep = 1, remove=F)

# add full treatment name
snag.dbh.data$trt.full = rep(0, nrow(snag.dbh.data))
for (i in 1:6) { snag.dbh.data$trt.full[which(snag.dbh.data$trt == trts[i])] = trt.names[i] }

# isolate all dead trees
snag.data = snag.dbh.data[which(snag.dbh.data$live == "D"),]

# make dataframe for diameter measurements
snag.diams = snag.data[-which(snag.data$dbh == "<2.5"),]

# calculate snag basal area (m2/ha) [cm2/m2 == m2/ha] and volume (m3/ha)
snag.diams$dbh = as.numeric(snag.diams$dbh) # cm
snag.diams$basal.area = (pi * (snag.diams$dbh/2)^2) / (12 * 48) #m2/ha
snag.diams$vol.min = snag.diams$basal.area * (4.5 / ftPerM) #m3/ha

# calculate snag C storage
snag.C.content.decay.class3 = 0.48 #% for standing angiosperms
snag.density.decay.class3 = 0.35 #g/cm3 for standing species
snag.diams$taxon.c.content = rep(0, nrow(snag.diams))
snag.diams = snag.diams %>% separate(species, c("genus","spp"), sep = " ", remove=F)
snag.diams$family = rep(NA, nrow(snag.diams))
for (i in 1:nrow(snag.diams)) {
  if (snag.diams$species[i] != "Unknown") {
    snag.diams$family[i] = genus.families[[snag.diams$genus[i]]]
    df.sp.id = which(dead.wood.c.df$species == snag.diams$species[i])
    df.gen.id = which(dead.wood.c.df$genus.resolved == snag.diams$genus[i]) 
    df.fam.id = which(dead.wood.c.df$family.resolved == snag.diams$family[i]) 
    if (length(df.sp.id) > 0) {
      snag.diams[i,"taxon.c.content"] = mean(dead.wood.c.df[df.sp.id,"tissue.c"])/100
    } else if (length(df.gen.id > 0)) {
      snag.diams[i,"taxon.c.content"] = mean(dead.wood.c.df[df.gen.id,"tissue.c"])/100
    } else if (length(df.fam.id > 0)) {
      snag.diams[i,"taxon.c.content"] = mean(dead.wood.c.df[df.fam.id,"tissue.c"])/100
    } else {
      snag.diams[i,"taxon.c.content"] = snag.C.content.decay.class3
    }
  } else {
    snag.diams[i,"taxon.c.content"] = snag.C.content.decay.class3
  }
}
snag.diams$carbon.min = snag.diams$vol.min * snag.diams$taxon.c.content * snag.density.decay.class3 # Mg/ha

# calculate total basal area within each plot
snag.sum = snag.diams %>% 
           group_by(trt, trt.full, num) %>% 
           summarize(snag.area = sum(basal.area),
                     snag.minvol = sum(vol.min),
                     snag.carbonmin = sum(carbon.min))

# add zeros to diameter sums
for (i in 1:6) {
  for (j in 1:3) {
    trt.plot.id = which(snag.sum$trt == trts[i] & snag.sum$num == j)
    if (length(trt.plot.id) == 0) {
      zero.row = data.frame(matrix(ncol=6, nrow=1))
      colnames(zero.row) = colnames(snag.sum)
      zero.row[1, c("trt","trt.full","num")] = c(trts[i], trt.names[i], j)
      zero.row[1, c("snag.area","snag.minvol","snag.carbonmin")] = rep(0,3)
      snag.sum = rbind(snag.sum, zero.row)
    }
  }
}
trt.sort = sort(snag.sum$trt, index.return=T)
snag.sum = snag.sum[trt.sort$ix,]

# estimate volume and C stocks of dead stems
snag.counts = snag.data[which(snag.data$dbh == "<2.5"),]
snag.counts = snag.counts[,-which(colnames(snag.counts) %in% c("plot","redo","live","dbh","notes"))]
colnames(snag.counts)[4] = "dead.stem.count"
snag.counts$dead.stem.count = snag.counts$dead.stem.count / (12 * 48) * (10^4) # count/ha
snag.counts$dead.stem.volmin = snag.counts$dead.stem.count * pi * ((1.25/100)^2) * (4.5/ftPerM) #m3/ha

# calculate dead stem C storage
snag.counts$taxon.c.content = rep(0, nrow(snag.counts))
snag.counts = snag.counts %>% separate(species, c("genus","spp"), sep = " ", remove=F)
snag.counts$family = rep(NA, nrow(snag.counts))
for (i in 1:nrow(snag.counts)) {
  if (snag.counts$species[i] != "Unknown") {
    snag.counts$family[i] = genus.families[[snag.counts$genus[i]]]
    df.sp.id = which(dead.wood.c.df$species == snag.counts$species[i])
    df.gen.id = which(dead.wood.c.df$genus.resolved == snag.counts$genus[i]) 
    df.fam.id = which(dead.wood.c.df$family.resolved == snag.counts$family[i]) 
    if (length(df.sp.id) > 0) {
      snag.counts[i,"taxon.c.content"] = mean(dead.wood.c.df[df.sp.id,"tissue.c"])/100
    } else if (length(df.gen.id > 0)) {
      snag.counts[i,"taxon.c.content"] = mean(dead.wood.c.df[df.gen.id,"tissue.c"])/100
    } else if (length(df.fam.id > 0)) {
      snag.counts[i,"taxon.c.content"] = mean(dead.wood.c.df[df.fam.id,"tissue.c"])/100
    } else {
      snag.counts[i,"taxon.c.content"] = snag.C.content.decay.class3
    }
  } else {
    snag.counts[i,"taxon.c.content"] = snag.C.content.decay.class3
  }
}
snag.counts$dead.stem.carbonmin = snag.counts$dead.stem.volmin * snag.counts$taxon.c.content * snag.density.decay.class3 #Mg/ha

# sum snag counts over species
snag.count.sum = snag.counts %>% 
                 group_by(trt, trt.full, num) %>% 
                 summarize(dead.stem.count = sum(dead.stem.count),
                           dead.stem.volmin = sum(dead.stem.volmin),
                           dead.stem.carbonmin = sum(dead.stem.carbonmin))

# add zeros to small snag count data
for (i in 1:6) {
  for (j in 1:3) {
    trt.plot.id = which(snag.count.sum$trt == trts[i] & snag.count.sum$num == j)
    if (length(trt.plot.id) == 0) {
      zero.row = data.frame(matrix(ncol=6, nrow=1))
      colnames(zero.row) = c("trt","trt.full","num","dead.stem.count","dead.stem.volmin","dead.stem.carbonmin")
      zero.row[1, c("trt","trt.full","num")] = c(trts[i], trt.names[i], j)
      zero.row[1, c("dead.stem.count","dead.stem.volmin","dead.stem.carbonmin")] = rep(0,3)
      snag.count.sum = rbind(snag.count.sum, zero.row)
    }
  }
}
trt.sort = sort(snag.count.sum$trt, index.return=T)
snag.count.sum = snag.count.sum[trt.sort$ix,]

# add to biomass dataframe
all.bm.data = right_join(all.bm.data,
                         snag.sum[,c("trt","trt.full","num","snag.area","snag.minvol","snag.carbonmin")], 
                         by=c("trt","trt.full","num"))
all.bm.data = right_join(all.bm.data,
                         snag.count.sum[,c("trt","trt.full","num","dead.stem.count","dead.stem.volmin","dead.stem.carbonmin")], 
                         by=c("trt","trt.full","num"))

# species-level estimates of snag C stocks
snag.sp.sum = snag.diams %>% 
              group_by(trt, trt.full, num, species) %>% 
              summarize(snag.carbonmin = sum(carbon.min))

# reshape data frame
snag.wide = pivot_wider(snag.sp.sum, 
                        id_cols = c(trt, num),
                        names_from = species, 
                        values_from = snag.carbonmin)
snag.wide[is.na(snag.wide)] = 0

# add zeros to diameter sums
snag.spp = colnames(snag.wide)[3:13]
for (i in 1:6) {
  for (j in 1:3) {
    trt.plot.id = which(snag.wide$trt == trts[i] & snag.wide$num == j)
    if (length(trt.plot.id) == 0) {
      zero.row = data.frame(matrix(ncol=length(colnames(snag.wide)), nrow=1))
      colnames(zero.row) = colnames(snag.wide)
      zero.row[1, c("trt","num")] = c(trts[i], j)
      zero.row[1, snag.spp] = rep(0, length(snag.spp))
      snag.wide = rbind(snag.wide, zero.row)
    }
  }
}
trt.sort = sort(snag.wide$trt, index.return=T)
snag.wide = snag.wide[trt.sort$ix,]

# reshape dataframe
snag.melt = melt(snag.wide, 
                 id.vars=c("trt","num"), 
                 variable.name="species",
                 value.name="snag.carbonmin")
trt.sort = sort(snag.melt$trt, index.return=T)
snag.melt = snag.melt[trt.sort$ix,]

# write cleaned snag data to csv
write.csv(snag.melt, "TreeData/Snag_Carbon_Stocks_By_Species.csv", row.names=F)

### clean live stem count data ---- ############################################
### double-checked #############################################################

# read in 2022 dbh data
dbh.data.2022 = read_excel('TreeData/DBH_08022022.xlsx')
colnames(dbh.data.2022) = c("plot","spp","dbh","stem.count")
stem.data.2022 = dbh.data.2022[which(dbh.data.2022$dbh == "<3"), 
                               c("plot","spp","stem.count")]

# get full species names
allo.df = read.csv("TreeData/Databases/Jenkins2004.csv")
stem.data.2022$species = rep(0, nrow(stem.data.2022))
for (i in 1:nrow(stem.data.2022)) {
  spp.id = which(allo.df$spp == stem.data.2022$spp[i])
  stem.data.2022$species[i] = paste(allo.df[spp.id,"genus"], allo.df[spp.id,"species"], " ")
}

# split plot name column
stem.data.2022 = stem.data.2022 %>% separate(plot, c("trt", "num"), sep = 1, remove=F)

# read in 2023 dbh data
stem.data.2023 = snag.dbh.data[which(snag.dbh.data$dbh == "<2.5" & snag.dbh.data$live == "L"),]
stem.data.2023 = stem.data.2023[, c("plot","trt","num","species","stem.count")]

# combine 2022 and 2023 data
redo.plots = unique(snag.dbh.data$plot[which(snag.dbh.data$redo == "Y")])
redo.ind = which(stem.data.2022$plot %in% redo.plots) # no overlap, so safe to do simple combination
all.stem.data = rbind(stem.data.2022[,c("plot","trt","num","species","stem.count")], 
                      stem.data.2023[,c("plot","trt","num","species","stem.count")])
all.stem.data$species = trimws(all.stem.data$species)


# used allodb to estimate aboveground biomass for each species
all.stem.data = all.stem.data %>% separate(species, c("genus","spp"), sep=" ", remove=F)
all.stem.data$agb = get_biomass(dbh = rep(1.25, nrow(all.stem.data)),
                                genus = all.stem.data$genus, 
                                species = all.stem.data$spp, 
                                coords = c(-90.1835, 41.5542))/1000 # [kg -> Mg]
all.stem.data$agb = all.stem.data$agb * all.stem.data$stem.count # Mg
all.stem.data$agb.ha = all.stem.data$agb / (12 * 48) * (10^4) # Mg/ha

# identify corresponding C density for each species
all.stem.data$c.content = rep(0,nrow(all.stem.data))
all.stem.data$level.c = rep(0,nrow(all.stem.data))
for (i in 1:nrow(all.stem.data)) {
  df.sp.id = which(live.wood.c.df$species == all.stem.data$species[i])
  df.gen.id = which(live.wood.c.df$genus.resolved == all.stem.data$genus[i])
  if (length(df.sp.id) == 0) {
    all.stem.data$c.content[i] = mean(live.wood.c.df[df.gen.id,"tissue.c"])/100
    all.stem.data$level.c[i] = "genus"
  } else {
    all.stem.data$c.content[i] = mean(live.wood.c.df[df.sp.id,"tissue.c"])/100
    all.stem.data$level.c[i] = "species"
  }
}

# estimate stem C content
all.stem.data$live.stem.count = all.stem.data$stem.count / (12 * 48) * (10^4) # count/ha
all.stem.data$live.stem.carbon = all.stem.data$agb.ha * all.stem.data$c.content # Mg/ha

# add full treatment name
all.stem.data$trt.full = rep(0, nrow(all.stem.data))
for (i in 1:6) { all.stem.data$trt.full[which(all.stem.data$trt == trts[i])] = trt.names[i] }

# sum data by plot
live.stem.sum = all.stem.data %>% 
                group_by(trt, trt.full, num) %>% 
                summarize(live.stem.count = sum(live.stem.count),
                          live.stem.carbon = sum(live.stem.carbon))

# add zeros
for (i in 1:6) {
  for (j in 1:3) {
    trt.plot.id = which(live.stem.sum$trt == trts[i] & live.stem.sum$num == j)
    if (length(trt.plot.id) == 0) {
      zero.row = data.frame(matrix(ncol=5, nrow=1))
      colnames(zero.row) = colnames(live.stem.sum)
      zero.row[1, c("trt","trt.full","num")] = c(trts[i], trt.names[i], j)
      zero.row[1, c("live.stem.count","live.stem.carbon")] = rep(0,2)
      live.stem.sum = rbind(live.stem.sum, zero.row)
    }
  }
}
trt.sort = sort(live.stem.sum$trt, index.return=T)
live.stem.sum = live.stem.sum[trt.sort$ix,]

# add to biomass dataframe
all.bm.data = right_join(all.bm.data,
                         live.stem.sum[,c("trt","trt.full","num","live.stem.count","live.stem.carbon")], 
                         by=c("trt","trt.full","num"))

# sum data by species
live.stem.sp.sum = all.stem.data %>% 
                   group_by(trt, trt.full, num, species) %>% 
                   summarize(live.stem.carbon = sum(live.stem.carbon))

# reshape data frame
stem.sp.wide = pivot_wider(live.stem.sp.sum, 
                           id_cols = c(trt, num),
                           names_from = species, 
                           values_from = live.stem.carbon)
stem.sp.wide[is.na(stem.sp.wide)] = 0

# add zeros to diameter sums
live.stem.spp = colnames(stem.sp.wide)[3:8]
for (i in 1:6) {
  for (j in 1:3) {
    trt.plot.id = which(stem.sp.wide$trt == trts[i] & stem.sp.wide$num == j)
    if (length(trt.plot.id) == 0) {
      zero.row = data.frame(matrix(ncol=length(colnames(stem.sp.wide)), nrow=1))
      colnames(zero.row) = colnames(stem.sp.wide)
      zero.row[1, c("trt","num")] = c(trts[i], j)
      zero.row[1, live.stem.spp] = rep(0, length(live.stem.spp))
      stem.sp.wide = rbind(stem.sp.wide, zero.row)
    }
  }
}
trt.sort = sort(stem.sp.wide$trt, index.return=T)
stem.sp.wide = stem.sp.wide[trt.sort$ix,]

# reshape dataframe
stem.sp.melt = melt(stem.sp.wide, 
                    id.vars=c("trt","num"), 
                    variable.name="species",
                    value.name="live.stem.carbon")
trt.sort = sort(stem.sp.melt$trt, index.return=T)
stem.sp.melt = stem.sp.melt[trt.sort$ix,]

# write cleaned live stem data to csv
write.csv(stem.sp.melt, "TreeData/LiveStem_Carbon_Stocks_By_Species.csv", row.names=F)

### read in calculated tree and soil C stock data ---- #########################
### double-checked #############################################################

# tree data
C.stock.data = read.csv("TreeData/Biomass_C_stocks_by_plot.csv", header=T)
C.stock.2022 = C.stock.data[which(C.stock.data$year==2022), c("trt","trt.full","num","abC_ha3")]
C.stock.2022$num = as.character(C.stock.2022$num)
colnames(C.stock.2022)[4] = c("live.woody.C.stock")
trt.sort = sort(C.stock.2022$trt, index.return=T)
C.stock.2022 = C.stock.2022[trt.sort$ix,]
all.bm.data = right_join(all.bm.data, C.stock.2022, by=c("trt","trt.full","num"))

# soil data
ave.soil.data = read.csv("SoilData/Soil_Data_Averaged_by_Plot_June2023.csv", header=T)
ave.soil.data$num = as.character(ave.soil.data$num)
all.bm.data = right_join(all.bm.data, 
                         ave.soil.data[,c("trt","trt.full","num","SOC","TC")], 
                         by=c("trt","trt.full","num"))

#### read in C stocks for herbaceous litter, FWD, and biomass --- ##############
### double-checked #############################################################

# read in understory data
understory.data = read.csv("VegetationData/CN_Summary_2022.csv", header=T)

# convert g/m2 to Mg/ha [(g/m2) * (1 Mg/10^6 g) * (10^4 m2/ha)]
understory.data$C.Mg.ha = understory.data$C.g.m2 / 100

# add herbaceous/fwd C stock data to woody dataframe
types = c("FWD","HL","PA","HJ","mix")
all.bm.data$num = as.integer(all.bm.data$num)
for (t in types) {
  t.ind = which(understory.data$type == t)
  if (t == "FWD" || t == "HL") {
    t.bm.data = understory.data[t.ind, c("trt","trt.full","num","C.Mg.ha","CN")] 
    colnames(t.bm.data)[4:5] = c(paste(t, ".C.Mg.ha", sep=""), paste(t, ".CN.ratio", sep=""))
  } else {
    t.bm.data = understory.data[t.ind, c("trt","trt.full","num","C.Mg.ha")]
    colnames(t.bm.data)[4] = paste(t, ".C.Mg.ha", sep="")
  }
  all.bm.data = right_join(all.bm.data, t.bm.data)
}

# sum together herbaceous C stocks
all.bm.data = all.bm.data %>% 
              mutate(total.C.Mg.ha = PA.C.Mg.ha + HJ.C.Mg.ha + mix.C.Mg.ha)

## calculate weighted biomass C & N concentrations

# fill in zeros of C/N data
understory.data$C[which(is.na(understory.data$C))] = 0
understory.data$N[which(is.na(understory.data$N))] = 0
all.bm.data$bm.CN.ratio = rep(0, nrow(all.bm.data))
for (i in 1:6) {
  for (j in 1:3) {
    t = trts[i]
    t.df = understory.data[which(understory.data$trt == t & understory.data$num == j),] 
    t.C.conc = t.df[which(t.df$type %in% c("mix","PA","HJ")),"C"]
    t.N.conc = t.df[which(t.df$type %in% c("mix","PA","HJ")),"N"]
    t.bm = t.df[which(t.df$type %in% c("mix","PA","HJ")),"sum.mass"]
    t.C.weighted = as.numeric((t.C.conc %*% t.bm)/sum(t.bm))
    t.N.weighted = as.numeric((t.N.conc %*% t.bm)/sum(t.bm))
    t.num.id = which(all.bm.data$trt == t & all.bm.data$num == j)
    all.bm.data[t.num.id,"bm.CN.ratio"] = t.C.weighted/t.N.weighted
  }
}

#### read in species richness data --- #########################################
### double-checked #############################################################

# herbaceous species data
sp.richness.df = read.csv("VegetationData/Species_Richness_by_Plot_and_Year.csv")

# average richness by year
sp.aves = sp.richness.df %>% 
          group_by(Trt, Num) %>%
          summarize(N = mean(N))
colnames(sp.aves) = c("trt","num","N.herb")

# tree species data
C.stock.sp.df = read.csv("TreeData/Biomass_C_stocks_by_species.csv", header=T)
C.stock.sp.df.2022 = C.stock.sp.df[which(C.stock.sp.df$year == 2022),]

# make dataframe for unique tree species
n.t = 6; n.p = 3; nums = seq(1,3)
n.tree.sp.df = data.frame(matrix(nrow=n.t*n.p, ncol=3))
colnames(n.tree.sp.df) = c("trt","num","N.tree")
k = 1
for (t in 1:n.t) {
  for (n in 1:n.p) {
    tn.id = which(C.stock.sp.df.2022$trt == trts[t] & C.stock.sp.df.2022$num == nums[n])
    tn.df = C.stock.sp.df.2022[tn.id,]
    tn.sp = length(sort(unique(tn.df$spp)))
    n.tree.sp.df[k, c("trt","num","N.tree")] = c(trts[t], as.numeric(nums[n]), as.numeric(tn.sp))
    k = k + 1
  }
}

# add richness values to biomass data frame
all.bm.data$num = as.numeric(all.bm.data$num)
n.tree.sp.df$num = as.numeric(n.tree.sp.df$num)
n.tree.sp.df$N.tree = as.numeric(n.tree.sp.df$N.tree)
all.bm.data = right_join(all.bm.data, sp.aves, by=c("trt","num"))
all.bm.data = right_join(all.bm.data, n.tree.sp.df, by=c("trt","num"))

# write data to file
write.csv(all.bm.data, "TreeData/Clean_Cstocks_Richness.csv", row.names=F)
