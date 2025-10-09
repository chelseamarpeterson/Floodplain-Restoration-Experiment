path_to_soil_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo/Soil_Analysis"
setwd(path_to_soil_folder)

library(tidyr)
library(dplyr)
library(reshape2)
library(tidyverse)

### script that brings together all files with quadrat-level soil data into
### one clean dataframe and averages quadrat level SOC stock data

################################################################################
## data upload and preparation

# treatments
trt.letters = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")
n.t = length(trt.letters)

# read in all soil data
temp.data = read.csv("Raw_Data/Soil_Temperature_June2023.csv")[,seq(1,6)]
bd.data = read.csv("Raw_Data/Bulk_Density_June2023.csv")
chem.data = read.csv("Raw_Data/Waypoint_Results_June2023.csv")
cn.data = read.csv("Raw_Data/TC_TOC_Results_2025.csv")
pom.data = read.csv("Raw_Data/POM_Results_2025.csv")
quad.data = read.csv("Raw_Data/Quadrat_GIS_Data_PolygonMean_WGS1984Aux_2020Lidar.csv")

# split plot name column
bd.data = bd.data %>% separate(Treatment_Plot, c("Treatment","Plot"), sep = 1, remove = T)
chem.data = chem.data %>% separate(Treatment_Plot, c("Treatment","Plot"), sep = 1, remove = T)
cn.data = cn.data %>% separate(Treatment_Plot_Quadrat, c("Treatment_Plot","Quadrat"), sep = " ", remove = T)
cn.data = cn.data %>% separate(Treatment_Plot, c("Treatment","Plot"), sep = 1, remove = T)
pom.data = pom.data %>% separate(Treatment_Plot_Quadrat, c("Treatment_Plot","Quadrat"), sep = " ", remove = T)
pom.data = pom.data %>% separate(Treatment_Plot, c("Treatment","Plot"), sep = 1, remove = T)

# update column names
colnames(temp.data)[5] = "temperature"
colnames(chem.data)[seq(8,26)] = c("texture.class","ph","p",
                                   "k","ca","mg","som.percent",
                                   "n.release","no3.n","nh4.n","cec",
                                   "k.sat","ca.sat","mg.sat","h.sat",
                                   "k.meq","ca.meq","mg.meq","h.meq")
colnames(bd.data)[c(5,20,21,22)] = c("Depth","gravimetric.moisture",
                                     "volumetric.moisture","bulk.density")
colnames(cn.data)[seq(4,8)] = c("bulk.n.percent","bulk.c.percent","bulk.cn.ratio","toc.n.percent","toc.percent")
colnames(pom.data)[seq(5,7)] = c("pom.mass","poc.n.percent","poc.percent")
colnames(quad.data)[5:6] = c("elevation.mean","vegetation.height.mean")

# assume second R1 Q1 is equal to R1 Q3
bd.ind1 = which((bd.data$Treatment == "R" & bd.data$Plot == "1") & (bd.data$Quadrat == "Q1" & bd.data$Depth == "0-15"))[2]
bd.ind2 = which((bd.data$Treatment == "R" & bd.data$Plot == "1") & (bd.data$Quadrat == "Q3" & bd.data$Depth == "0-15"))
bd.data[bd.ind2, colnames(bd.data)[6:23]] = bd.data[bd.ind1, colnames(bd.data)[6:23]]
bd.data = bd.data[-bd.ind1,]

# compute moisture and bulk density for 0-30 cm
bd.sum = bd.data %>% 
         group_by(Treatment, Plot, Quadrat) %>%
         summarise(mass.water = sum(`Mass.of.water..g.`),
                   sum.soil = sum(`Mass.of.dry.soil..g.`),
                   sum.vol = sum(`Corrected.volume..cm3.`),
                   mass.coarse = sum(`Mass.of.coarse.roots.rocks..g.`))
bd.sum$gravimetric.moisture = bd.sum$mass.water/bd.sum$sum.soil
bd.sum$volumetric.moisture = bd.sum$mass.water/bd.sum$sum.vol
bd.sum$bulk.density = bd.sum$sum.soil/bd.sum$sum.vol
bd.sum$coarse.material = bd.sum$mass.coarse
bd.sum = bd.sum[,-which(colnames(bd.sum) %in% c("mass.water","sum.soil","sum.vol","mass.coarse"))]

# average the POM reps
pom.ave = pom.data[,c("Treatment","Plot","Quadrat","Rep","poc.percent","pom.mass")] %>%
          group_by(Treatment, Plot, Quadrat) %>%
          summarize(poc.percent = mean(poc.percent*pom.mass/10)) # POM-C (%) = (g POM-C/g POM)*(g POM/g soil)

# add full treatment name
temp.data$full.treatment.name = rep(0, nrow(temp.data))
for (i in 1:n.t) { temp.data$full.treatment.name[which(temp.data$Treatment == trt.letters[i])] = trt.names[i] }

# combine all bulk density, moisture, temperature, C/N, vegetation height, elevation, and other chemical data
temp.data$Plot = as.character(temp.data$Plot)
quad.data$Plot = as.character(quad.data$Plot)
soil.data = left_join(temp.data[,c("full.treatment.name","Treatment","Plot","Quadrat","temperature")], 
                      bd.sum, 
                      by=c("Treatment","Plot","Quadrat"))
soil.data = left_join(soil.data, 
                      cn.data[,-which(colnames(cn.data) %in% c("toc.n.percent"))],
                      by=c("Treatment","Plot","Quadrat"))
soil.data = left_join(soil.data, 
                      chem.data[,-which(colnames(chem.data) %in% c("Number","Note","n.release","k.sat","ca.sat","mg.sat","h.sat"))], 
                      by=c("Treatment","Plot","Quadrat"))
soil.data = left_join(soil.data, 
                      pom.ave,
                      by=c("Treatment","Plot","Quadrat"))
soil.data = left_join(soil.data,
                      quad.data[,-which(colnames(quad.data) %in% c("Trt_Plot_Quad","Shape_Length","Shape_Area"))], 
                      by=c("Treatment","Plot","Quadrat"))

################################################################################
## clean up aggregate data 

# load mass data
ag.data = read.csv("Raw_Data/Aggregate_Masses_2023.csv")

# update columns
key.cols = c("Size.class","Treatment_Plot_Quadrat","Mass.of.sample.w.o.fragments..g.")
ag.data = ag.data[,colnames(ag.data) %in% key.cols]
colnames(ag.data) = c("size","treatment_plot_quadrat","soil.mass")

# split sample plotber into plot and quadrat, plot into trt and plot
ag.data = ag.data %>% separate(treatment_plot_quadrat, c("treatment_plot","quadrat"), sep=" ", remove=T)
ag.data = ag.data %>% separate(treatment_plot, c("treatment","plot"), sep=1, remove=F)

# size classes, plots, and quadrats
sizes = c("fPOM",">4.75 mm",">2 mm",">250 um",">53 um","<=53 um")         
trt_plots = sort(unique(ag.data$treatment_plot))
quads = sort(unique(ag.data$quadrat))
n.s = length(sizes)
n.tp = length(trt_plots)

# make new dataframe to estimate mass of fraction <53 um
ag.df = data.frame(matrix(nrow=0, ncol=6))
ag.cols = c("size","treatment","plot","quadrat","soil.mass","rel.mass")
colnames(ag.df) = ag.cols
ag.data$rel.mass = rep(0, nrow(ag.data))
for (p in trt_plots) {
  for (q in quads) {
    # plot and quadrat id
    pq.id = which(ag.data$treatment_plot == p & ag.data$quadrat == q)
    
    # get total mass and estimate <53 mass
    pq.full.mass = ag.data[which(ag.data[pq.id,"size"] == "Full"),"soil.mass"]
    pq.53.mass = pq.full.mass - sum(ag.data[pq.id[2:length(pq.id)],"soil.mass"])
    
    # get subset of dataframe for plot p and quad q
    pq.df = ag.data[pq.id[2:length(pq.id)], ag.cols]
    
    # add row to p and q dataframe
    pq.row = data.frame(matrix(nrow=1, ncol=6))
    colnames(pq.row) = ag.cols
    pq.row[1, "size"] = "<=53 um"
    pq.row[1, ag.cols[2:4]] = pq.df[1, ag.cols[2:4]]
    pq.row[1, c("soil.mass","rel.mass")] = c(pq.53.mass, 0)
    pq.df = rbind(pq.df, pq.row)
    
    # estimate relative masses
    pq.df[,"rel.mass"] = pq.df[,"soil.mass"]/as.numeric(pq.full.mass)*100
    
    # append new dataframe
    ag.df = rbind(ag.df, pq.df)   
  }
}

# make all categorical variables into factors
ag.df$quadrat = as.factor(ag.df$quadrat)
ag.df$size = factor(ag.df$size, levels=sizes)

# convert to long form
ag.df$size.lab = rep("", nrow(ag.df))
size.labs = c("fPOM","g475mm","g2mm","g250um","g53um","l53um")
for (i in 1:n.s) {
  size.id = which(ag.df$size == sizes[i])
  ag.df[size.id,"size.lab"] = size.labs[i]
}
ag.new = ag.df[,c("treatment","plot","quadrat","rel.mass","size.lab")]
ag.wide = pivot_wider(ag.new,id_cols = c(treatment, plot, quadrat),
                      names_from = size.lab, values_from = rel.mass)

# compute mean weight diameter from aggregate size distribution
sizes = c(0.053, 0.250, 2, 4.75, 8)
mean.diams = (sizes[2:5]-sizes[1:4])/2 + sizes[1:4]
ag.wide$MWD = (mean.diams[1]*ag.wide$g53um + mean.diams[2]*ag.wide$g250um + mean.diams[3]*ag.wide$g2mm + mean.diams[4]*ag.wide$g475mm)/100

# add aggregate data to soil dataframe
colnames(soil.data) = tolower(colnames(soil.data))
colnames(ag.wide) = tolower(colnames(ag.wide))
soil.data = right_join(soil.data, ag.wide, by=c("treatment","plot","quadrat"))

################################################################################
## final quadrat-level data cleaning and writing to file

# calculate inorganic and mineral-associated organic c
soil.data$bulk.c.percent = pmax(soil.data$bulk.c.percent, soil.data$toc.percent)
soil.data$tic.percent = soil.data$bulk.c.percent - soil.data$toc.percent
soil.data$moc.percent = pmax(soil.data$toc.percent - soil.data$poc.percent, 0)

# calculate POC:MOC ratio
soil.data$poc.moc.ratio = soil.data$poc.percent/soil.data$moc.percent

# write all soil data at quadrat level to csv
write.csv(soil.data, "Clean_Data/Soil_Data_by_Quadrat_June2023.csv", row.names=F)

