path_to_soil_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo/Soil_Analysis"
setwd(path_to_soil_folder)

library(ggplot2)
library(tidyr)
library(patchwork)
library(dplyr)
library(reshape2)
library(viridis)
library(tidyverse)

################################################################################
## Step 1: Data upload and preparation

# treatment names
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

# read in all soil data
temp.data = read.csv("Raw_Data/Soil_Temperature_June2023.csv")[,seq(1,4)]
bd.data = read.csv("Raw_Data/Bulk_Density_June2023.csv")
chem.data = read.csv("Raw_Data/Waypoint_Results_June2023.csv")
cn.data = read.csv("Raw_Data/TC_TOC_Results_2025.csv")
pom.data = read.csv("Raw_Data/POM_Results_2025.csv")

# split plot name column
temp.data = temp.data %>% separate(Plot, c("Trt", "Num"), sep = 1, remove = F)
bd.data = bd.data %>% separate(Plot, c("Trt", "Num"), sep = 1, remove = F)
chem.data = chem.data %>% separate(Plot, c("Trt", "Num"), sep = 1, remove = F)
cn.data = cn.data %>% separate(Plot, c("Plot", "Quad"), sep = " ", remove = T)
cn.data = cn.data %>% separate(Plot, c("Trt", "Num"), sep = 1, remove = F)
pom.data = pom.data %>% separate(Plot, c("Plot", "Quad"), sep = " ", remove = T)
pom.data = pom.data %>% separate(Plot, c("Trt", "Num"), sep = 1, remove = F)

# update column names
colnames(temp.data)[6] = "Temperature"
colnames(chem.data)[seq(9,27)] = c("Texture.Class","pH","P",
                                   "K","Ca","Mg","SOM",
                                   "N.rel","NO3","NH4","CEC",
                                   "K.sat","Ca.sat","Mg.sat","H.sat",
                                   "K.meq","Ca.meq","Mg.meq","H.meq")
colnames(bd.data)[c(5,6,21,22,23)] = c("Quad","Depth","Gravimetric.Moisture",
                                       "Volumetric.Moisture","Bulk.Density")
colnames(cn.data)[seq(5,9)] = c("bulk.N","bulk.C","bulk.CN","ton","toc")
colnames(pom.data)[seq(6,8)] = c("pom.mass","pon","poc")

# assume second R1 Q1 is equal to R1 Q3
bd.ind1 = which((bd.data$Plot == "R1" & bd.data$Quad == "Q1") & bd.data$Depth == "0-15")[2]
bd.ind2 = which((bd.data$Plot == "R1" & bd.data$Quad == "Q3") & bd.data$Depth == "0-15")
bd.data[bd.ind2, colnames(bd.data)[7:24]] = bd.data[bd.ind1, colnames(bd.data)[7:24]]
bd.data = bd.data[-bd.ind1,]

# compute moisture and bulk density for 0-30 cm
bd.sum = bd.data %>% 
         group_by(Plot,Trt, Num, Quad) %>%
         summarise(mass.water = sum(`Mass.of.water..g.`),
                   sum.soil = sum(`Mass.of.dry.soil..g.`),
                   sum.vol = sum(`Corrected.volume..cm3.`))
bd.sum$Gravimetric.Moisture = bd.sum$mass.water/bd.sum$sum.soil
bd.sum$Volumetric.Moisture = bd.sum$mass.water/bd.sum$sum.vol
bd.sum$Bulk.Density = bd.sum$sum.soil/bd.sum$sum.vol
bd.sum = bd.sum[,-which(colnames(bd.sum) %in% c("mass.water","sum.soil","sum.vol"))]

# average the POM reps
pom.ave = pom.data[,c("Plot","Trt","Num","Quad","Rep","pon","poc","pom.mass")] %>%
          group_by(Plot, Trt, Num, Quad) %>%
          summarize(pon = mean(pon*pom.mass/10),
                    poc = mean(poc*pom.mass/10)) # POM-C (%) = (g POM-C/g POM)*(g POM/g soil)

# sort chem data by treatment and quadrat
chem.sort = sort(chem.data$Plot, index.return=T)
chem.data = chem.data[chem.sort$ix,]

# combine all bulk density, moisture, temperature, C/N, and other chemical data
soil.data = right_join(temp.data[,-which(colnames(temp.data) %in% c("Date","BD15.Starting.Depth"))], bd.sum, 
                       by=c("Plot","Trt","Num","Quad"))
soil.data = right_join(soil.data, cn.data, by=c("Plot","Trt","Num","Quad"))
soil.data = right_join(soil.data, chem.data[,-which(colnames(chem.data) %in% c("Number","Note"))], 
                       by=c("Plot","Trt","Num","Quad"))
soil.data = right_join(soil.data, pom.ave, by=c("Plot","Trt","Num","Quad"))

# remove miscellaneous or redundant columns
leave.out = c("Texture.Class","H.sat","K.sat","Ca.sat","Mg.sat","N.rel")
soil.data = soil.data[,-which(colnames(soil.data) %in% leave.out)]

################################################################################
## Step 2: Clean up aggregate data 

# load mass data
ag.data = read.csv("Raw_Data/Aggregate_Masses_2023.csv")

# update columns
key.cols = c("Boat.number","Size.class","Sample.number","Mass.of.sample.w.o.fragments..g.")
ag.data = ag.data[,colnames(ag.data) %in% key.cols]
colnames(ag.data) = c("boat","size","sample","soil.mass")

# histogram of fPOM values
#hist(ag.data[which(ag.data$size == "fPOM"),"soil.mass"]$soil.mass)

# split sample number into plot and quadrat, plot into trt and num
ag.data = ag.data %>% separate(sample, c("plot","quad"), sep = " ", remove = F)
ag.data = ag.data %>% separate(plot, c("trt","num"), sep = 1, remove = F)

# size classes, plots, and quadrats
sizes = c("fPOM",">4.75 mm",">2 mm",">250 um",">53 um","<=53 um")         
plots = sort(unique(ag.data$plot))
quads = sort(unique(ag.data$quad))
n.s = length(sizes)
n.p = length(plots)

# make new dataframe to estimate mass of fraction <53 um
ag.df = data.frame(matrix(nrow=0, ncol=7))
ag.cols = c("size","plot","trt","num","quad","soil.mass","rel.mass")
colnames(ag.df) = ag.cols
ag.data$rel.mass = rep(0, nrow(ag.data))
for (p in plots) {
  for (q in quads) {
    # plot and quadrat id
    pq.id = which(ag.data$plot == p & ag.data$quad == q)
    
    # get total mass and estimate <53 mass
    pq.full.mass = ag.data[which(ag.data[pq.id,"size"] == "Full"),"soil.mass"]
    pq.53.mass = pq.full.mass - sum(ag.data[pq.id[2:length(pq.id)],"soil.mass"])
    
    # get subset of dataframe for plot p and quad q
    pq.df = ag.data[pq.id[2:length(pq.id)], ag.cols]
    
    # add row to p and q dataframe
    pq.row = data.frame(matrix(nrow=1, ncol=7))
    colnames(pq.row) = ag.cols
    pq.row[1, "size"] = "<=53 um"
    pq.row[1, ag.cols[2:6]] = pq.df[1, ag.cols[2:6]]
    pq.row[1, c("soil.mass","rel.mass")] = c(pq.53.mass, 0)
    pq.df = rbind(pq.df, pq.row)
    
    # estimate relative masses
    pq.df[,"rel.mass"] = pq.df[,"soil.mass"]/as.numeric(pq.full.mass)*100
    
    # append new dataframe
    ag.df = rbind(ag.df, pq.df)   
  }
}

# make all categorical variables into factors
ag.df$quad = as.factor(ag.df$quad)
ag.df$size = factor(ag.df$size, levels=sizes)

# convert to long form
ag.df$size.lab = rep("", nrow(ag.df))
size.labs = c("fPOM","g475mm","g2mm","g250um","g53um","l53um")
for (i in 1:n.s) {
  size.id = which(ag.df$size == sizes[i])
  ag.df[size.id,"size.lab"] = size.labs[i]
}
ag.new = ag.df[,c("trt","num","quad","rel.mass","size.lab")]
ag.wide = pivot_wider(ag.new,id_cols = c(trt, num, quad),
                      names_from = size.lab, values_from = rel.mass)

# compute mean weight diameter from aggregate size distribution
sizes = c(0.053, 0.250, 2, 4.75, 8)
mean.diams = (sizes[2:5]-sizes[1:4])/2 + sizes[1:4]
ag.wide$MWD = (mean.diams[1]*ag.wide$g53um + mean.diams[2]*ag.wide$g250um + mean.diams[3]*ag.wide$g2mm + mean.diams[4]*ag.wide$g475mm)/100

# add aggregate data to soil dataframe
colnames(soil.data) = tolower(colnames(soil.data))
colnames(ag.wide) = tolower(colnames(ag.wide))
soil.data = right_join(soil.data, ag.wide, by=c("trt","num","quad"))

################################################################################
## Step 3: Final data cleaning and saving to file

# add full treatment name
soil.data$trt.full = rep(0, nrow(soil.data))
for (i in 1:6) { soil.data$trt.full[which(soil.data$trt == trts[i])] = trt.names[i] }

# calculate inorganic and mineral-associated organic c
soil.data$bulk.c = pmax(soil.data$bulk.c, soil.data$toc)
soil.data$tic = soil.data$bulk.c - soil.data$toc
soil.data$maoc = pmax(soil.data$toc - soil.data$poc, 0)

# calculate POC:MAOC ratio
soil.data$poc_maoc = soil.data$poc/soil.data$maoc

# write all soil data at quadrat level to csv
id.vars = c("plot","trt","trt.full","num","quad")
soil.vars = colnames(soil.data)[-which(colnames(soil.data) %in% id.vars)]
write.csv(soil.data[,c(id.vars,soil.vars)], "Clean_Data/All_Soil_Data_2023.csv", row.names=F)

################################################################################
## Step 4: calculate averages at the plot level

# estimate SOC stocks (Mg/ha)
soil.data$tic.stock = soil.data$tic * soil.data$bulk.density * 30
soil.data$maoc.stock = soil.data$maoc * soil.data$bulk.density * 30
soil.data$poc.stock = soil.data$poc * soil.data$bulk.density * 30
soil.data$soc.stock = soil.data$toc * soil.data$bulk.density * 30
soil.data$tc.stock = soil.data$bulk.c * soil.data$bulk.density * 30

# average results across plots and write to csv
soil.aves = soil.data %>% 
            group_by(plot, trt, trt.full, num) %>%
            summarize(maoc.stock = mean(maoc.stock),
                      poc.stock = mean(poc.stock),
                      soc.stock = mean(soc.stock),
                      tic.stock = mean(tic.stock),
                      tc.stock = mean(tc.stock),
                      cn.ratio = mean(bulk.cn),
                      sand = mean(sand),
                      silt = mean(silt),
                      clay = mean(clay),
                      bd = mean(bulk.density),
                      moisture = mean(gravimetric.moisture),
                      temperature = mean(temperature),
                      ph = mean(ph),
                      p = mean(p),
                      k = mean(k),
                      ca = mean(ca),
                      mg = mean(mg),
                      no3 = mean(no3),
                      nh4 = mean(nh4),
                      cec = mean(cec))
write.csv(soil.aves, "Clean_Data/Soil_Data_Averaged_by_Plot_June2023.csv", row.names=F)
