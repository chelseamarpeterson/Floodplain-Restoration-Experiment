path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo"
setwd(path_to_repo)

library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(reshape2)

################################################################################
# script that brings together all biomass carbon stocks and writes the results to a CSV file

# conversions
ft_per_m = 3.28084
cm_per_in = 2.54
m_per_km = 10^3
cm_per_m = 100
m2_per_ha = 10^4
kg_per_Mg = 10^3

# constants
C1 = 21.1659
fwd.C.content = 0.484 #% (averaged from measured FWD data)
cwd.C.content.decay.class3 = 0.481 # %
cwd.density.decay.class3 = 0.26 # g/cm^3
snag.C.content.decay.class3 = 0.48 #% for standing angiosperms
snag.density.decay.class3 = 0.35 #g/cm3 for standing species
min.snag.height = 4.5 # ft
med.stem.width = 1.25 # cm
med.fwd.diam = 5 # cm

# plot dimensions 
plot.length = 48 #m 
plot.width = 12 #m 
plot.area.m2 = plot.length*plot.width #m2
plot.area.ha = plot.area/m2_per_ha #ha
plot.length.ft = plot.length * ft_per_m #ft

# treatment labels
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

# list of species families
family.df = read.csv("Tree_Analysis/Tree_Databases/Tree_Families.csv")
n.f = nrow(family.df)
genus.families = list()
for (i in 1:n.f) { genus.families[[family.df$Genus[i]]] = family.df$Family[i] }
families = family.df$Family

## read in woody C data 
wood.c.df = read.csv("Tree_Analysis/Tree_Databases/Doraisami_2021_Wood_C_Database.csv")

# isolate rows for dead material
dead.wood.c.df = wood.c.df[which(wood.c.df$dead.alive == "dead"),]
dead.wood.c.df = dead.wood.c.df[which(dead.wood.c.df$tissue == "stem"),]
dead.wood.c.df = dead.wood.c.df %>% separate(binomial.resolved, c("genus","spp"), sep="_", remove=F)
dead.wood.c.df$species = paste(dead.wood.c.df$genus, dead.wood.c.df$spp, sep= " ")

# isolate rows for live material
live.wood.c.df = wood.c.df[which(wood.c.df$dead.alive == "alive" & wood.c.df$growth.form == "tree"),]
live.wood.c.df = live.wood.c.df[which(live.wood.c.df$position == "standing" & live.wood.c.df$tissue == "stem"),]
live.wood.c.df = live.wood.c.df %>% separate(binomial.resolved, c("genus","spp"), sep="_", remove=F)
live.wood.c.df$species = paste(live.wood.c.df$genus, live.wood.c.df$spp, sep= " ")

################################################################################
# step: clean coarse and fine woody debris data

# load dataframe and change columne names
cwd.data = read.csv("Tree_Analysis/Raw_Data/CWD_June2023.csv")
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
fwd.counts$fwd.vol = (1/plot.length.ft) * (C1*(pi^2)/8) * fwd.counts$fwd.count * ((med.fwd.diam/cm_per_in)^2) #m3/ha
fwd.counts$int.fwd.carbon = fwd.counts$fwd.vol * fwd.C.content * cwd.density.decay.class3 # m3/ha * g/cm3 * g/g = Mg/ha

# make data frame for coarse woody debris (>7.5 cm)
cwd.diams = cwd.data[-which(is.na(cwd.data$diameter)), -which(colnames(cwd.data) == "fwd.count")]
cwd.diams = cwd.diams %>% separate(species, c("genus","spp"), sep=" ", remove=F)

# calculate average downed angiosperm C content [Table 3, Harmon et al. (2013)]
C.contents = c(47.2,47.8,47.7,48.1,47.4,47.3)
n.samples = c(8,21,25,32,17,4)
cwd.C.ave = sum((C.contents*n.samples))/sum(n.samples)

# add column for potential C content by taxon
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
cwd.diams$cwd.vol = (1/plot.length.ft) * (C1*(pi^2)/8) * (cwd.diams$diameter/cm_per_in)^2 #m3/ha
cwd.diams$cwd.carbon = cwd.diams$cwd.vol * cwd.diams$taxon.c.content * cwd.density.decay.class3 # Mg/ha

# sum carbon stocks (Mg/ha) by plot
cwd.sum = cwd.diams %>% 
          group_by(trt, trt.full, num) %>% 
          summarize(cwd.carbon = sum(cwd.carbon))

# add zeros to CWD areas and volumes
for (i in 1:6) {
  for (j in 1:3) {
    trt.plot.id = which(cwd.sum$trt == trts[i] & cwd.sum$num == j)
    if (length(trt.plot.id) == 0) {
      zero.row = data.frame(matrix(nrow=1, ncol=4))
      colnames(zero.row) = colnames(cwd.sum)
      zero.row[1, c("trt","trt.full","num")] = c(trts[i], trt.names[i], j)
      zero.row[1, "cwd.carbon"] = 0
      cwd.sum = rbind(cwd.sum, zero.row)
    }
  }
}
trt.sort = sort(cwd.sum$trt, index.return=T)
cwd.sum = cwd.sum[trt.sort$ix,]

# sum carbon stocks (Mg/ha) by species
cwd.diams$species[which(cwd.diams$species == "")] = "Unknown"
cwd.sp.sum = cwd.diams %>% 
             group_by(trt, trt.full, num, species) %>% 
             summarize(cwd.carbon = sum(cwd.carbon))

# reshape data frame
cwd.sp.wide = pivot_wider(cwd.sp.sum, 
                          id_cols = c(trt, num),
                          names_from = species, 
                          values_from = cwd.carbon)
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
                   value.name="cwd.carbon")

# write cleaned cwd data to csv
write.csv(cwd.sp.melt, "Tree_Analysis/Clean_Data/CWD_Carbon_Stocks_By_Species.csv", row.names=F)

################################################################################
# step: clean snag data

# load dataframe and change column names
snag.dbh.data = read.table("Tree_Analysis/Raw_Data/DBH_Snags_June2023.csv", header=T, sep=",")
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

# split species column
snag.diams = snag.diams %>% separate(species, c("genus","spp"), sep = " ", remove=F)

# calculate snag volume (m3/ha)
snag.diams$dbh = as.numeric(snag.diams$dbh) # cm
snag.diams$basal.area = (pi * (snag.diams$dbh/2)^2) / plot.area.m2 #m2/ha = cm2/m2
snag.diams$vol.min = snag.diams$basal.area * (min.snag.height / ft_per_m) #m3/ha

# calculate snag C storage
snag.diams$taxon.c.content = rep(0, nrow(snag.diams))
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
           summarize(snag.carbonmin = sum(carbon.min))

# add zeros to diameter sums
for (i in 1:6) {
  for (j in 1:3) {
    trt.plot.id = which(snag.sum$trt == trts[i] & snag.sum$num == j)
    if (length(trt.plot.id) == 0) {
      zero.row = data.frame(matrix(nrow=1, ncol=4))
      colnames(zero.row) = colnames(snag.sum)
      zero.row[1, c("trt","trt.full","num")] = c(trts[i], trt.names[i], j)
      zero.row[1, "snag.carbonmin"] = 0
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
snag.counts$dead.stem.count = snag.counts$dead.stem.count / plot.area.ha # count/ha
snag.counts$dead.stem.volmin = snag.counts$dead.stem.count * pi * ((med.stem.width/cm_per_m)^2) * (min.snag.height/ft_per_m) # m3/ha

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
                 summarize(dead.stem.carbonmin = sum(dead.stem.carbonmin))

# add zeros to small snag count data
for (i in 1:6) {
  for (j in 1:3) {
    trt.plot.id = which(snag.count.sum$trt == trts[i] & snag.count.sum$num == j)
    if (length(trt.plot.id) == 0) {
      zero.row = data.frame(matrix(nrow=1, ncol=4))
      colnames(zero.row) = c("trt","trt.full","num","dead.stem.carbonmin")
      zero.row[1, c("trt","trt.full","num")] = c(trts[i], trt.names[i], j)
      zero.row[1, "dead.stem.carbonmin"] = 0
      snag.count.sum = rbind(snag.count.sum, zero.row)
    }
  }
}
trt.sort = sort(snag.count.sum$trt, index.return=T)
snag.count.sum = snag.count.sum[trt.sort$ix,]

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
write.csv(snag.melt, "Tree_Analysis/Clean_Data/Snag_Carbon_Stocks_By_Species.csv", row.names=F)

################################################################################
# step: clean live stem count data

# read in 2022 dbh data
dbh.data.2022 = read.csv('Tree_Analysis/Raw_Data/DBH_August2022.csv')
colnames(dbh.data.2022) = c("plot","spp","dbh","stem.count")
stem.data.2022 = dbh.data.2022[which(dbh.data.2022$dbh == "<3"), c("plot","spp","stem.count")]

# get full species names
allo.df = read.csv("Tree_Analysis/Tree_Databases/Chojnacky2014.csv")
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
colnames(all.stem.data)[4] = "spp"
all.stem.data = all.stem.data %>% separate(spp, c("genus","species"), sep=" ", remove=F)
all.stem.data$abg = get_biomass(dbh = rep(med.stem.width, nrow(all.stem.data)),
                               genus = all.stem.data$genus, 
                               species = all.stem.data$spp, 
                               coords = c(-90.1835, 41.5542)) * all.stem.data$stem.count / kg_per_Mg # [kg -> Mg]
all.stem.data$abg.ha = all.stem.data$abg / plot.area.ha # Mg/ha

# get root ratios from Chojnacky and estimate belowground biomass
all.stem.data$r = exp(allo.df[1,"cr.b0"] + allo.df[1,"cr.b1"]*log(med.stem.width)) # ratio = exp(b0 + b1*log(dbh [cm]))
all.stem.data$bg.ha = all.stem.data$abg.ha * all.stem.data$r # Mg/ha

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
all.stem.data$abg.live.stem.carbon = all.stem.data$abg.ha * all.stem.data$c.content # Mg/ha
all.stem.data$bg.live.stem.carbon = all.stem.data$bg.ha * all.stem.data$c.content # Mg/ha
all.stem.data$tot.live.stem.carbon = all.stem.data$abg.live.stem.carbon + all.stem.data$bg.live.stem.carbon

# add full treatment name
all.stem.data$trt.full = rep(0, nrow(all.stem.data))
for (i in 1:6) { all.stem.data$trt.full[which(all.stem.data$trt == trts[i])] = trt.names[i] }

# sum data by plot
live.stem.sum = all.stem.data %>% 
                group_by(trt, trt.full, num) %>% 
                summarize(abg.live.stem.carbon = sum(abg.live.stem.carbon),
                          bg.live.stem.carbon = sum(bg.live.stem.carbon),
                          tot.live.stem.carbon = sum(tot.live.stem.carbon))

# add zeros
for (i in 1:6) {
  for (j in 1:3) {
    trt.plot.id = which(live.stem.sum$trt == trts[i] & live.stem.sum$num == j)
    if (length(trt.plot.id) == 0) {
      zero.row = data.frame(matrix(nrow=1, ncol=6))
      colnames(zero.row) = colnames(live.stem.sum)
      zero.row[1, c("trt","trt.full","num")] = c(trts[i], trt.names[i], j)
      zero.row[1, c("abg.live.stem.carbon","bg.live.stem.carbon","tot.live.stem.carbon")] = 0
      live.stem.sum = rbind(live.stem.sum, zero.row)
    }
  }
}
trt.sort = sort(live.stem.sum$trt, index.return=T)
live.stem.sum = live.stem.sum[trt.sort$ix,]

# sum data by species
live.stem.sp.sum = all.stem.data %>% 
                   group_by(trt, trt.full, num, spp) %>% 
                   summarize(tot.live.stem.carbon = sum(tot.live.stem.carbon))

# reshape data frame
stem.sp.wide = pivot_wider(live.stem.sp.sum, 
                           id_cols = c(trt, num),
                           names_from = spp, 
                           values_from = tot.live.stem.carbon)
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
                    value.name="tot.live.stem.carbon")
trt.sort = sort(stem.sp.melt$trt, index.return=T)
stem.sp.melt = stem.sp.melt[trt.sort$ix,]

# write cleaned live stem data to csv
write.csv(stem.sp.melt, "Tree_Analysis/Clean_Data/LiveStem_Carbon_Stocks_By_Species.csv", row.names=F)

################################################################################
# step: calculate hypothetical C stock if frax. pen. were alive

snag.frax.pen = snag.diams[which(snag.diams$species == "Fraxinus pennsylvanica"),]
snag.frax.pen$hyp.biomass.live = get_biomass(dbh = snag.frax.pen$dbh, 
                                             genus = snag.frax.pen$genus, 
                                             species = snag.frax.pen$spp, 
                                             coords = c(-90.1835, 41.5542))/kg_per_Mg
snag.frax.pen$hyp.biomass.live.ha = snag.frax.pen$hyp.biomass.live / plot.area.ha # Mg/ha
snag.frax.pen$c.content = rep(0,nrow(snag.frax.pen))
snag.frax.pen$level.c = rep(0,nrow(snag.frax.pen))
for (i in 1:nrow(snag.frax.pen)) {
  df.sp.id = which(live.wood.c.df$species == snag.frax.pen$species[i])
  df.gen.id = which(live.wood.c.df$genus.resolved == snag.frax.pen$genus[i])
  if (length(df.sp.id) == 0) {
    snag.frax.pen$c.content[i] = mean(live.wood.c.df[df.gen.id,"tissue.c"])/100
    snag.frax.pen$level.c[i] = "genus"
  } else {
    snag.frax.pen$c.content[i] = mean(live.wood.c.df[df.sp.id,"tissue.c"])/100
    snag.frax.pen$level.c[i] = "species"
  }
}
snag.frax.pen$hyp.Cstock.live.ha = snag.frax.pen$hyp.biomass.live.ha*snag.frax.pen$c.content
snag.frax.sum = snag.frax.pen %>% 
  group_by(trt, trt.full, num) %>% 
  summarize(snag.frax.liveC = sum(hyp.Cstock.live.ha)) 
for (i in 1:6) {
  for (j in 1:3) {
    trt.plot.id = which(snag.frax.sum$trt == trts[i] & snag.frax.sum$num == j)
    if (length(trt.plot.id) == 0) {
      zero.row = data.frame(matrix(nrow=1, ncol=4))
      colnames(zero.row) = colnames(snag.frax.sum)
      zero.row[1, c("trt","trt.full","num")] = c(trts[i], trt.names[i], j)
      zero.row[1, "snag.frax.liveC"] = 0
      snag.frax.sum = rbind(snag.frax.sum, zero.row)
    }
  }
}
trt.sort = sort(snag.frax.sum$trt, index.return=T)
snag.frax.live = snag.frax.sum[trt.sort$ix,]

snag.frax.dead = snag.wide[, c("trt","num","Fraxinus pennsylvanica")]
colnames(snag.frax.dead) = c("trt","num","snag.frax.deadC")
snag.frax.sum = left_join(snag.frax.live, snag.frax.dead,
                          by = c("trt","num"))

################################################################################
# step: read in calculated tree C stock data ---- ##############################

# tree data
C.stock.data = read.csv("Tree_Analysis/Clean_Data/WoodyBiomass_C_Stocks_By_Plot.csv", header=T)
C.stock.2022 = C.stock.data[which(C.stock.data$year==2022), c("trt","trt.full","num","abC_ha3","bgC_ha3","bmC_ha3")]
C.stock.2022$num = as.character(C.stock.2022$num)
colnames(C.stock.2022)[4:6] = c("aboveground.woody.c.stock","belowground.woody.c.stock","total.woody.c.stock")
trt.sort = sort(C.stock.2022$trt, index.return=T)
C.stock.2022 = C.stock.2022[trt.sort$ix,]

# combine CWD/FWD/snag/woody biomass dataframes 
all.bm.data = right_join(fwd.counts[,c("trt","trt.full","num","int.fwd.carbon")], 
                         cwd.sum[,c("trt","trt.full","num","cwd.carbon")], 
                         by=c("trt","trt.full","num"))
all.bm.data = right_join(all.bm.data,
                         snag.sum[,c("trt","trt.full","num","snag.carbonmin")], 
                         by=c("trt","trt.full","num"))
all.bm.data = right_join(all.bm.data,
                         snag.frax.sum[,c("trt","trt.full","num","snag.frax.liveC","snag.frax.deadC")], 
                         by=c("trt","trt.full","num"))
all.bm.data = right_join(all.bm.data,
                         snag.count.sum[,c("trt","trt.full","num","dead.stem.carbonmin")], 
                         by=c("trt","trt.full","num"))
all.bm.data = right_join(all.bm.data,
                         live.stem.sum[,c("trt","trt.full","num","abg.live.stem.carbon","bg.live.stem.carbon","tot.live.stem.carbon")], 
                         by=c("trt","trt.full","num"))
all.bm.data = right_join(all.bm.data, C.stock.2022, by=c("trt","trt.full","num"))

################################################################################
# step: read in C stocks for herbaceous litter, FWD, and biomass

# read in understory data
understory.data = read.csv("Understory_Analysis/Clean_Data/CN_Stock_Summary_2022.csv", header=T)

# convert g/m2 to Mg/ha [(g/m2) * (1 Mg/10^6 g) * (10^4 m2/ha)]
understory.data$C.Mg.ha = understory.data$C.g.m2 / 100
understory.data$N.Mg.ha = understory.data$N.g.m2 / 100

# reshape data frame
understory.wide = pivot_wider(understory.data, 
                              id_cols = c(trt, num),
                              names_from = type, 
                              values_from = C.Mg.ha)

# update column names then join
colnames(understory.wide)[3:8] = c("HB.C.Mg.ha","MB.C.Mg.ha","PAB.C.Mg.ha",
                                   "HJB.C.Mg.ha","HL.C.Mg.ha","FWD.C.Mg.ha")

###############################################################################
# step: calculate C:N ratios of FWD, HL, and HB

# separate data by type
understory.fwd = understory.data[which(understory.data$type == "Fine.woody.debris"),]
understory.hl = understory.data[which(understory.data$type == "Herbaceous.litter"),]
understory.hb = understory.data[which(understory.data$type == "Herbaceous.biomass"),]

# FWD C:N ratio
understory.wide$FWD.CN.ratio = understory.fwd$C.Mg.ha / understory.fwd$N.Mg.ha
understory.wide$FWD.CN.ratio[which(is.na(understory.wide$FWD.CN.ratio))] = 0

# HL C:N ratio
understory.wide$HL.CN.ratio = understory.hl$C.Mg.ha / understory.hl$N.Mg.ha
understory.wide$HL.CN.ratio[which(is.na(understory.hl$HL.CN.ratio))] = 0

# HN C:N ratio
understory.wide$HB.CN.ratio = understory.hb$C.Mg.ha / understory.hb$N.Mg.ha
understory.wide$HB.CN.ratio[which(is.na(understory.hb$HB.CN.ratio))] = 0

## join understory data to tree and woody debris data
understory.wide$num = as.character(understory.wide$num)
all.bm.data = right_join(all.bm.data, understory.wide, by=c("trt","num"))

# write all C stock and C:N ratio to file
write.csv(all.bm.data, "Tree_Analysis/Clean_Data/All_Vegetation_C_Stocks_By_Plot.csv", row.names=F)