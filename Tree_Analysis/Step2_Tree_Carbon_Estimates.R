path_to_tree_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo/Tree_Analysis"
setwd(path_to_tree_folder)

library(allodb)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)

# conversions
m2_per_ha = 10^4
kg_per_Mg = 10^3

# constants
tree.C.content = 0.48 # %

# plot dimensions 
plot.length = 48 #m 
plot.width = 12 #m 
plot.area.m2 = plot.length*plot.width #m2
plot.area.ha = plot.area/m2_per_ha #ha

################################################################################
# step 1: estimate above- and belowground biomass from diameter-at-breast height (DBH) data (cm) 

# define data-related variables
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")
years = c(2022, 2023)
n.y = length(years)

# list of genus families
family.df = read.csv("Tree_Databases/Tree_Families.csv")
n.f = nrow(family.df)
genus.families = list()
for (i in 1:n.f) { genus.families[[family.df$Genus[i]]] = family.df$Family[i] }

# read in all DBH data
dbh.all = read.csv("Clean_Data/DBH_Data_Clean_2023.csv", header=T)

# make list for unique species
uni.spp = sort(unique(dbh.all$spp))
n.sp = length(uni.spp)

# read in allometric equation matrices
allo.df1 = read.csv("Tree_Databases/Jenkins2004.csv", header=T); rownames(allo.df1) = allo.df1$spp
allo.df2 = read.csv("Tree_Databases/Chojnacky2014.csv", header=T); rownames(allo.df2) = allo.df2$spp

# estimate biomass with allometric eqns in 1: Jenkins (2004) & 2: Chojnacky (2014)
new.mat = data.frame(matrix(nrow = length(dbh.all$plot), ncol = 6))
colnames(new.mat) = c("ab1","ab2","r1","r2","bg1","bg2")
dbh.all = cbind(dbh.all,new.mat)
for (sp in uni.spp) {
  sp.ind = which(dbh.all$spp == sp)
  
  # Jenkins
  dbh.all[sp.ind,"ab1"] = exp(allo.df1[sp,"ab.b0"] + allo.df1[sp,"ab.b1"]*log(dbh.all[sp.ind,"dbh"]))/kg_per_Mg # bm = exp(b0 + b1*log(dbh [cm])) [kg -> Mg]
  dbh.all[sp.ind,"r1"] = exp(allo.df1[sp,"cr.b0"] + allo.df1[sp,"cr.b1"]/dbh.all[sp.ind,"dbh"]) # ratio = exp(b0 + b1/(dbh [cm]))
  dbh.all[sp.ind,"bg1"] = dbh.all[sp.ind,"r1"]*dbh.all[sp.ind,"ab1"] 
  
  # Chojnacky
  dbh.all[sp.ind,"ab2"] = exp(allo.df2[sp,"ab.b0"] + allo.df2[sp,"ab.b1"]*log(dbh.all[sp.ind,"dbh"]))/kg_per_Mg # bm = exp(b0 + b1*log(dbh [cm])) [kg -> Mg]
  dbh.all[sp.ind,"r2"] = exp(allo.df2[sp,"cr.b0"] + allo.df2[sp,"cr.b1"]*log(dbh.all[sp.ind,"dbh"])) # ratio = exp(b0 + b1*log(dbh [cm]))
  dbh.all[sp.ind,"bg2"] = dbh.all[sp.ind,"r2"]*dbh.all[sp.ind,"ab2"]
}

# estimate biomass with allodb
dbh.all$ab3 = get_biomass(dbh = dbh.all$dbh, 
                          genus = dbh.all$genus, 
                          species = dbh.all$species, 
                          coords = c(-90.1835, 41.5542))/kg_per_Mg # [kg -> Mg]
dbh.all$bg3 = dbh.all$ab3 * dbh.all$r2
mean(dbh.all$r2)

# check for NAs
sum(is.na(dbh.all[,c("ab1","ab2","r1","r2","bg1","bg2","ab3","bg3")]))

# sum above & belowground biomass by species within each plot
biomass.by.species = dbh.all %>% 
  group_by(redo, year, plot, trt, trt.full, num, spp, genus, species) %>% 
  summarise(ab1=sum(ab1), ab2=sum(ab2), ab3=sum(ab3), 
            bg1=sum(bg1), bg2=sum(bg2), bg3=sum(bg3))

# divide biomass by total plot area (Mg/ha)
bm.vars = c("ab1","ab2","ab3","bg1","bg2","bg3")
bm.vars.ha = c("ab_ha1","ab_ha2","ab_ha3","bg_ha1","bg_ha2","bg_ha3")
biomass.by.species[, bm.vars] = biomass.by.species[, bm.vars]/plot.area.ha
colnames(biomass.by.species)[which(colnames(biomass.by.species) %in% bm.vars)] = bm.vars.ha

################################################################################
# step 2: estimate carbon stocks from biomass (Mg C/ha)

## estimate C density by species

# clean wood database
wood.c.df = read.csv("Tree_Databases/Doraisami_2021_Wood_C_Database.csv")
live.wood.c.df = wood.c.df[which(wood.c.df$dead.alive == "alive" & wood.c.df$growth.form == "tree"),]
live.stem.c.df = live.wood.c.df[which(live.wood.c.df$tissue == "stem"),]
live.stem.c.df = live.stem.c.df %>% separate(binomial.resolved, c("genus","spp"), sep="_", remove=F)
live.stem.c.df$species = paste(live.stem.c.df$genus, live.stem.c.df$spp, sep= " ")

# get best-available density for each species
biomass.by.species$taxon.c.content = rep(0, nrow(biomass.by.species))
biomass.by.species$binomial = paste(biomass.by.species$genus, biomass.by.species$species, sep=" ")
biomass.by.species$family = rep(NA, nrow(biomass.by.species))
for (i in 1:nrow(biomass.by.species)) {
  if (biomass.by.species$binomial[i] != "") {
    biomass.by.species$family[i] = genus.families[[biomass.by.species$genus[i]]]
    df.sp.id = which(live.stem.c.df$species == biomass.by.species$binomial[i])
    df.gen.id = which(live.stem.c.df$genus.resolved == biomass.by.species$genus[i]) 
    df.fam.id = which(live.stem.c.df$family.resolved == biomass.by.species$family[i]) 
    if (length(df.sp.id) > 0) {
      biomass.by.species[i,"taxon.c.content"] = mean(live.stem.c.df[df.sp.id,"tissue.c"])/100
    } else if (length(df.gen.id > 0)) {
      biomass.by.species[i,"taxon.c.content"] = mean(live.stem.c.df[df.gen.id,"tissue.c"])/100
    } else if (length(df.fam.id > 0)) {
      biomass.by.species[i,"taxon.c.content"] = mean(live.stem.c.df[df.fam.id,"tissue.c"])/100
    } else {
      biomass.by.species[i,"taxon.c.content"] = tree.C.content
    }
  } else {
    biomass.by.species[i,"taxon.c.content"] = tree.C.content
  }
}

# estimate C stock by species
id.vars = c("redo","year","plot","trt","trt.full","num","spp","genus","species")
C.vars.ha = c("abC_ha1","abC_ha2","abC_ha3","bgC_ha1","bgC_ha2","bgC_ha3")
C.stocks.by.species = biomass.by.species[, c(id.vars, bm.vars.ha)]
C.stocks.by.species[,bm.vars.ha] = biomass.by.species[,bm.vars.ha] * biomass.by.species$taxon.c.content
colnames(C.stocks.by.species)[which(colnames(C.stocks.by.species) %in% bm.vars.ha)] = C.vars.ha

# sum C stocks by plot (Mg/ha)
C.stocks.by.plot <- C.stocks.by.species %>% 
                    group_by(redo, year, plot, trt, trt.full, num) %>% 
                    summarise(abC_ha1=sum(abC_ha1), abC_ha2=sum(abC_ha2), abC_ha3=sum(abC_ha3), 
                              bgC_ha1=sum(bgC_ha1), bgC_ha2=sum(bgC_ha2), bgC_ha3=sum(bgC_ha3))

# estimate total biomass C stocks per ha
C.stocks.by.plot = C.stocks.by.plot %>% 
                      group_by(redo, year, plot, trt, trt.full, num) %>%
                      summarise(abC_ha1 = abC_ha1, abC_ha2 = abC_ha2, abC_ha3 = abC_ha3,
                                bgC_ha1 = bgC_ha1, bgC_ha2 = bgC_ha2, bgC_ha3 = bgC_ha3,
                                bmC_ha1 = abC_ha1 + bgC_ha1,
                                bmC_ha2 = abC_ha2 + bgC_ha2,
                                bmC_ha3 = abC_ha3 + bgC_ha3)

################################################################################
# step 3: compare dbh and C stock results for all redo v. original

## make dataframe to compare C stock estimates before and after resampling

# redo plots from 2023
redo.ind = which(C.stocks.by.plot$redo == "Y")
redo.plot.df = C.stocks.by.plot[redo.ind,]
redo.plots = unique(redo.plot.df$plot)

# original plots from 2022
orig.ind = which((C.stocks.by.plot$redo == "N" & C.stocks.by.plot$year == 2022) & (C.stocks.by.plot$plot %in% redo.plots))
orig.plot.df = C.stocks.by.plot[orig.ind,]

# combine redo and original plots
redo.df.comp.plot = rbind(redo.plot.df, orig.plot.df)

# compare biomass results for each re-sample plot
ggplot(redo.df.comp.plot, aes(x=plot, y=bmC_ha3, fill=redo)) + 
       geom_bar(stat="identity", position=position_dodge()) + 
       labs(x="Plot", y="Woody Biomass C storage (Mg/ha)")

# compare distributions with and without redo
df.no.redo = C.stocks.by.plot[which(C.stocks.by.plot$redo =="N" & C.stocks.by.plot$year==2022),]
df.w.redo = C.stocks.by.plot[-c(orig.ind, which(C.stocks.by.plot$year==2013)),]
df.no.redo$with.redo = "F"
df.w.redo$with.redo = "T"
df.plot.comp = rbind(df.no.redo, df.w.redo)
ggplot(df.plot.comp, aes(x=factor(trt.full, levels=trt.names), y=bmC_ha3, fill=with.redo)) + 
       geom_boxplot() +
       labs(x="", y="Woody Biomass C storage (Mg/ha)")

################################################################################
# step 4: replace 2022 data with 2023 redo data

# for C stocks by plot 
C.stocks.by.plot.redo = C.stocks.by.plot[-orig.ind,]
C.stocks.by.plot.redo[which(C.stocks.by.plot.redo$year == 2023),"year"] = 2022

orig.ind.spp = which((C.stocks.by.species$redo == "N" & C.stocks.by.species$year == 2022) & (C.stocks.by.species$plot %in% redo.plots))
C.stocks.by.species.redo = C.stocks.by.species[-orig.ind.spp,]
C.stocks.by.species.redo[which(C.stocks.by.species.redo$year == 2023),"year"] = 2022

# estimate total biomass C stocks per ha
C.stocks.by.species.redo = C.stocks.by.species.redo %>% 
                            group_by(redo, year, plot, trt, trt.full, num, spp, genus, species) %>%
                            summarise(abC_ha1 = abC_ha1, abC_ha2 = abC_ha2, abC_ha3 = abC_ha3,
                                      bgC_ha1 = bgC_ha1, bgC_ha2 = bgC_ha2, bgC_ha3 = bgC_ha3,
                                      bmC_ha1 = abC_ha1 + bgC_ha1,
                                      bmC_ha2 = abC_ha2 + bgC_ha2,
                                      bmC_ha3 = abC_ha3 + bgC_ha3)

# write C stocks results to file
write.csv(C.stocks.by.plot.redo, "Clean_Data/WoodyBiomass_C_Stocks_By_Plot.csv", row.names=F)
write.csv(C.stocks.by.species.redo, "Clean_Data/WoodyBiomass_C_Stocks_By_Species.csv", row.names=F)
