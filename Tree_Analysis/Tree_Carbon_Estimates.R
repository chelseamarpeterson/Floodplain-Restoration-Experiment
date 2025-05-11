setwd('C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/TreeData')

library(allodb)
library(ggplot2)
library(patchwork)
library(dplyr)

### estimate above- and belowground biomass from dbh data (cm) 

# define data-related variables
years = c(2013, 2022, 2023)
n.y = length(years)

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

# read in all DBH data
dbh.all = read.csv("DBH_Data_Clean_2023.csv", header=T)

# make list for unique species
uni.spp = sort(unique(dbh.all$spp))
n.sp = length(uni.spp)

# read in allometric equation matrices
allo.df1 = read.csv("Databases/Jenkins2004.csv", header=T); rownames(allo.df1) = allo.df1$spp
allo.df2 = read.csv("Databases/Chojnacky2014.csv", header=T); rownames(allo.df2) = allo.df2$spp

# estimate biomass with allometric eqns in Jenkins (2004) & Chojnacky (2014)
new.mat = data.frame(matrix(nrow = length(dbh.all$plot), ncol = 6))
colnames(new.mat) = c("ab1","ab2","r1","r2","bg1","bg2")
dbh.all = cbind(dbh.all,new.mat)
for (sp in uni.spp) {
  sp.ind = which(dbh.all$spp == sp)
  dbh.all[sp.ind,"ab1"] = exp(allo.df1[sp,"ab.b0"] + allo.df1[sp,"ab.b1"]*log(dbh.all[sp.ind,"dbh"]))/1000 # bm = exp(b0 + b1*log(dbh [cm])) [kg -> Mg]
  dbh.all[sp.ind,"ab2"] = exp(allo.df2[sp,"ab.b0"] + allo.df2[sp,"ab.b1"]*log(dbh.all[sp.ind,"dbh"]))/1000 # bm = exp(b0 + b1*log(dbh [cm])) [kg -> Mg]
  dbh.all[sp.ind,"r1"] = exp(allo.df1[sp,"cr.b0"] + allo.df1[sp,"cr.b1"]/dbh.all[sp.ind,"dbh"]) # ratio = exp(b0 + b1/(dbh [cm]))
  dbh.all[sp.ind,"r2"] = exp(allo.df2[sp,"cr.b0"] + allo.df2[sp,"cr.b1"]*log(dbh.all[sp.ind,"dbh"])) # ratio = exp(b0 + b1*log(dbh [cm]))
  dbh.all[sp.ind,"bg1"] = dbh.all[sp.ind,"r1"]*dbh.all[sp.ind,"ab1"] 
  dbh.all[sp.ind,"bg2"] = dbh.all[sp.ind,"r2"]*dbh.all[sp.ind,"ab2"]
}

# estimate biomass with allodb
dbh.all$ab3 = get_biomass(dbh = dbh.all$dbh, 
                          genus = dbh.all$genus, 
                          species = dbh.all$species, 
                          coords = c(-90.1835, 41.5542))/1000 # [kg -> Mg]
dbh.all$bg3 = 0.20 * dbh.all$ab3

# check for NAs
sum(is.na(dbh.all[,c("ab1","ab2","r1","r2","bg1","bg2","ab3","bg3")]))

### estimate carbon stocks from biomass (cm) ----

# sum above & belowground biomass by species within each plot
biomass.by.species <- dbh.all %>% 
                      group_by(redo, year, plot, trt, trt.full, num, spp, genus, species) %>% 
                      summarise(ab1=sum(ab1), ab2=sum(ab2), ab3=sum(ab3), 
                                bg1=sum(bg1), bg2=sum(bg2), bg3=sum(bg3))

# divide biomass by total plot area (Mg/ha)
areas = c(5*50, 12*48, 12*48)/10^4 # 10^4 m2 = ha
bm.vars = c("ab1","ab2","ab3","bg1","bg2","bg3")
bm.vars.ha = c("ab_ha1","ab_ha2","ab_ha3","bg_ha1","bg_ha2","bg_ha3")
for (i in 1:n.y) {
  year.ind = which(biomass.by.species$year == years[i])
  biomass.by.species[year.ind, bm.vars] = biomass.by.species[year.ind, bm.vars]/areas[i]
}
colnames(biomass.by.species)[which(colnames(biomass.by.species) %in% bm.vars)] = bm.vars.ha

# estimate C density by species
wood.c.df = read.csv("Databases/Doraisami_2021_Wood_C_Database.csv")
live.wood.c.df = wood.c.df[which(wood.c.df$dead.alive == "alive" & wood.c.df$growth.form == "tree"),]
live.stem.c.df = live.wood.c.df[-which(live.wood.c.df$tissue %in% c("bark","twig","branch","coarseroot","fineroot","sapwood")),]
live.stem.c.df = live.stem.c.df %>% separate(binomial.resolved, c("genus","spp"), sep="_", remove=F)
live.stem.c.df$species = paste(live.stem.c.df$genus, live.stem.c.df$spp, sep= " ")

tree.C.content = 0.48 #%
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

### compare dbh and C stock results for all redo v. original ----

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

## replace 2022 data with 2023 redo data ----

# for C stocks by plot 
C.stocks.by.plot.redo = C.stocks.by.plot[-orig.ind,]
C.stocks.by.plot.redo[which(C.stocks.by.plot.redo$year == 2023),"year"] = 2022

orig.ind.spp = which((C.stocks.by.species$redo == "N" & C.stocks.by.species$year == 2022) & (C.stocks.by.species$plot %in% redo.plots))
C.stocks.by.species.redo = C.stocks.by.species[-orig.ind.spp,]
C.stocks.by.species.redo[which(C.stocks.by.species.redo$year == 2023),"year"] = 2022

# write C stocks results to file
write.csv(C.stocks.by.plot.redo, "Biomass_C_stocks_by_plot.csv", row.names=F)
write.csv(C.stocks.by.species.redo, "Biomass_C_stocks_by_species.csv", row.names=F)

################################################################################
### compare results across allometric equations
library(reshape2)

# replace redo plot biomass data
biomass.redo.ind = which(biomass.by.plot$redo == "N" & (biomass.by.plot$year == "2022" & biomass.by.plot$plot %in% redo.plots))
biomass.by.plot.redo = biomass.by.plot[-biomass.redo.ind,]

# aboveground biomass
trt.method.abg = biomass.by.plot.redo[,which(colnames(biomass.by.plot.redo) %in% c("year","trt","trt.full","ab_ha1","ab_ha2","ab_ha3"))]
trt.method.abg.melt = melt(trt.method.abg, 
                           id.vars = c("year","trt","trt.full"), 
                           variable.name = "method", value.name="ab_ha")

names = c("Jenkins (2003)","Chojnacky (2014)","allodb (2022)")
ggplot(trt.method.abg.melt, aes(x=factor(trt.full,level=trt.names), y=ab_ha, fill=factor(method))) +
       geom_boxplot() + 
       scale_fill_discrete(labels = names) +
       labs(x="Treatment", y="Aboveground biomass (Mg/ha)") +
       theme(legend.title = element_blank(), legend.text = element_text(size=14), 
             axis.title = element_text(size=16), axis.text = element_text(size=12))

trt.method.abg.melt$log_ab_ha = log(trt.method.abg.melt$ab_ha + 1)
aov.method.log.abg = aov(log_ab_ha ~ method*trt*year, data=trt.method.abg.melt)
summary(aov.method.log.abg)
plot(aov.method.log.abg)
hist(residuals(aov.method.log.abg))

aov.method.abg = aov(ab_ha ~ method*trt*year, data=trt.method.abg.melt)
summary(aov.method.abg)
plot(aov.method.abg)
hist(residuals(aov.method.abg))

# belowground biomass
trt.method.bg = biomass.by.plot[,which(colnames(biomass.by.plot) %in% c("year","trt","trt.full","bg_ha1","bg_ha2","bg_ha3"))]
trt.method.bg.melt = melt(trt.method.bg, id.vars = c("year","trt","trt.full"), 
                     variable.name = "method", value.name="bg_ha")

ggplot(trt.method.bg.melt, aes(x=factor(trt.full,level=trt.names), y=bg_ha, fill=factor(method))) +
        geom_boxplot() + theme(legend.title = element_blank()) +
        ylab("Belowground Biomass (Mg/ha)") + xlab("Treatment") + 
        scale_fill_discrete(labels = names) 

trt.method.bg.melt$log_bg_ha = log(trt.method.bg.melt$bg_ha + 1)
aov.method.log.bg = aov(log_bg_ha~method*year*trt, data=trt.method.bg.melt)
summary(aov.method.log.bg)
plot(aov.method.log.bg)

################################################################################
## do PCA on biomass carbon stocks by species

uni.spp = unique(C.stocks.by.species.redo$spp)
n.sp = length(uni.spp)
  
# make df
colors = c("coral1","darkgoldenrod","chartreuse3","cornflowerblue","darkorchid1","burlywood4")
years = c(2013, 2022)
n.yrs = 2
df.c = data.frame(matrix(nrow = n.trt * n.plot * n.yrs, ncol = n.sp + 5))
colnames(df.c) = c("year","trt","trt.full","num","color", uni.spp)
k = 1
for (y in years) {
  for (i in 1:n.trt) {
    t = trts[i]
    n = trt.names[i]
    c = colors[i]
    for (j in 1:n.plot) {
      df.c[k, c("year","trt","trt.full","num","color")] = c(y, t, n, j, c)
      p.ind = which((C.stocks.by.species.redo$year == y & C.stocks.by.species.redo$trt == t) & C.stocks.by.species.redo$num == j)
      p.uni.sp = unique(C.stocks.by.species.redo$spp[p.ind])
      col.ind = which(colnames(df.c) %in% p.uni.sp)
      if (length(p.ind) > 0) {
        df.c[k, col.ind] = C.stocks.by.species.redo[p.ind, "bmC_ha3"]
        df.c[k, which(is.na(df.c[k,]))] = 0
      } else {
        df.c[k, uni.spp] = 0
      }
      k = k + 1
    }
  }
}
            
## run PCA
library(ggfortify)
pca.C = prcomp(df.c[,-seq(1,5)], center=T, scale=T)
summary(pca.C)

# numerical PCA results
pca.C$rotation[,1:3]
pcs = data.frame(as.matrix(df.c[,-seq(1,5)]) %*% as.matrix(pca.C$rotation))
pcs$trt = df.c$trt
df.c$year = factor(df.c$year)

## auto plot
df.c$trt.full = factor(df.c$trt.full, levels=trt.names)
autoplot(pca.C, data=df.c, loadings=T, loadings.label=T, frame=T, 
         col="trt.full", shape="year", x=1, y=3)

## plot with better labels
CAloadings <- data.frame(Variables = rownames(pca.C$rotation), pca.C$rotation)

# add PCA scores to the dataset
df.c[, c('PC1', 'PC2')] = pca.C$x[, 1:2]

# save variable loadings in a separate dataset
rot = as.data.frame(pca.C$rotation[, 1:2])
rot$var = rownames(pca.C$rotation)

# rescale the loadings to fit nicely within the scatterplot of our data
mult = max(abs(df.c[, c('PC1', 'PC2')])) / max(abs(rot[, 1:2])) / 2
rot[, 1:2] = rot[, 1:2] * mult

# if there are many variables to plot, you can play with ggrepel 
library(ggrepel)  
library(ggalt)
library(plyr)

ggplot(data = rot, aes(x=0, y=0, xend = PC1, yend = PC2, label = var)) +
        geom_point(data = df.c, aes(x=PC1, y=PC2, color=factor(trt.full, levels=trt.names), shape=year), inherit.aes = FALSE, size=3) +
        geom_segment(color = 'red', arrow = arrow(length = unit(0.03, "npc"))) +
        geom_label_repel(aes(PC1 * 1, PC2 * 1)) +
        labs(x = "PC1 (13.0%)", y = "PC2 (11.5%)", color="") + 
        theme(legend.title = element_blank(), legend.text = element_text(size=10), 
              axis.title = element_text(size=14), axis.text = element_text(size=10))

################################################################################
## compare allometric equations

# calculate difference between two sets of allo eqns
ab_diff = trt.df$ab_C_ha1-trt.df$ab_C_ha2
bg_diff = trt.df$bg_C_ha1-trt.df$bg_C_ha2
par(mfrow=c(1,2))
hist(ab_diff)
hist(bg_diff)
sum(ab_diff)
sum(bg_diff)

# plot two sets of allo eqns
sp.colors = rep(0,10)
colors = c('cyan','orange','black','red','blue','green')
sp.order = c('Acersai','Aceneg','Quercus','Salix/Populus','Fraxinus','Carya')
for (i in 1:10) {
  sp = uni.sp[i]
  if (sp %in% oaks) {
    sp.colors[i] = 'black'
  } else if (sp %in% willows | sp=="Popdel") {
    sp.colors[i] = 'red'
  } else if (sp == "Fralan") {
    sp.colors[i] = 'blue'
  } else if (sp == "Carill") {
    sp.colors[i] = 'green'
  } else if (sp == "Aceneg") {
    sp.colors[i] = 'orange'
  } else if (sp == "Acesai") {
    sp.colors[i] = 'cyan'
  }
}
dhb.range = seq(2.5,55,2.5)
par(mfrow=c(1,2))
plot(dhb.range,allo.df1[uni.sp[1],"ab.b0"]+allo.df1[uni.sp[1],"ab.b1"]*dhb.range,
     type='n',ylim=c(0,2.6),xlab="DBH (cm)",ylab="Aboveground Biomass (Mg)",
     main="Jenkins et al. (2003)")
for (i in 1:10) {
  sp = uni.sp[i]
  if (sp %in% oaks) {
    lines(dhb.range,exp(allo.df1[sp,"ab.b0"]+allo.df1[sp,"ab.b1"]*log(dhb.range))/1000,
          col=sp.colors[i],lwd=2,lty=1)
  } else if (sp == "Carill") {
    lines(dhb.range,exp(allo.df1[sp,"ab.b0"]+allo.df1[sp,"ab.b1"]*log(dhb.range))/1000,
          col=sp.colors[i],lwd=3,lty=1)
  } else {
    lines(dhb.range,exp(allo.df1[sp,"ab.b0"]+allo.df1[sp,"ab.b1"]*log(dhb.range))/1000,
          col=sp.colors[i],lwd=2,lty=2)
  }
}
legend(legend=sp.order,col=colors,x="topleft",lwd=2,lty=2)

plot(dhb.range,allo.df2[uni.sp[1],"ab.b0"]+allo.df2[uni.sp[1],"ab.b1"]*dhb.range,
     type='n',ylim=c(0,2.6),xlab="DBH (cm)",ylab="",
     main="Chojnacky et al. (2014)")
for (i in 1:10) {
  sp = uni.sp[i]
  if (sp %in% oaks | sp == "Carill") {
    lines(dhb.range,exp(allo.df2[sp,"ab.b0"]+allo.df2[sp,"ab.b1"]*log(dhb.range))/1000,
          col=sp.colors[i],lwd=2,lty=1)
  } else {
    lines(dhb.range,exp(allo.df2[sp,"ab.b0"]+allo.df2[sp,"ab.b1"]*log(dhb.range))/1000,
          col=sp.colors[i],lwd=2,lty=2)
  }
}

