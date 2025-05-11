setwd("C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/SoilData")

library(ggplot2)
library(readxl)
library(tidyr)

# treatment names
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

# read in POM masses
pom.dat = read_excel("Size_Fractionation_Masses.xlsx", sheet="Samples")
cn.dat = read_excel("CN_Analysis_Results_Clean_May2024.xlsx")

# split plot column
pom.dat = pom.dat %>% separate(Plot, c("Trt", "Num"), sep=1, remove=F)
cn.dat = cn.dat %>% separate(Sample, c("Plot", "Quad"), sep=" ", remove=F)
cn.dat = cn.dat %>% separate(Plot, c("Trt", "Num"), sep=1, remove=F)

# update column names
colnames(pom.dat) = c("plot","trt","num","quad","rep","batch",
                      "bottle","start.mass","mesh.pom","mesh","pom","notes")
colnames(cn.dat) = c("well","sample","plot","trt","num","quad","N","C","CN")
pom.dat$maom = pom.dat$start.mass - pom.dat$pom

# add treatment name column
pom.dat$trt.full = pom.dat$trt
cn.dat$trt.full = cn.dat$trt
for (i in 1:6) {pom.dat$trt.full[which(pom.dat$trt == trts[i])] = trt.names[i]}
for (i in 1:6) {cn.dat$trt.full[which(cn.dat$trt == trts[i])] = trt.names[i]}

# make barplot
pom.dat$trt.full = factor(pom.dat$trt.full, levels=trt.names)
ggplot(pom.dat, aes(y=trt.full, x=pom)) + 
       geom_boxplot() + labs(y="", x="POM mass (g)")


cn.dat$trt.full = factor(cn.dat$trt.full, levels=trt.names)
ggplot(cn.dat, aes(y=trt.full, x=C)) + 
       geom_boxplot() + labs(y="", x="Carbon (%)")
       
ggplot(cn.dat, aes(y=trt.full, x=N)) + 
       geom_boxplot() + labs(y="", x="Nitrogen (%)")    

ggplot(cn.dat, aes(y=trt.full, x=CN)) + 
       geom_boxplot() + labs(y="", x="C/N ratio")  

