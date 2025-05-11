path_to_tree_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Public-Repo/Tree_Analysis"

library(tidyr)

### script that loads and cleans diameter-at-breast-height (DBH) data from 
### Joslin wetland mitigation site in Henry County

## define important variables
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")
years = c(2013, 2022)
n.trt = length(trts)
n.plot = 3
n.yrs = 2

## read in data from 2013, 2020, and 2023
setwd(path_to_tree_data)
data.2013 = read.csv('Raw_Data/DBH_2013.csv')
data.2022 = read.csv('Raw_Data/DBH_August2022.csv')
data.2023 = read.csv('Raw_Data/DBH_Snags_June2023.csv', header=T)

# update column names
colnames(data.2013) = c("plot","spp","dbh","basal.area","stem.count")
colnames(data.2022) = c("plot","spp","dbh","stem.count")
colnames(data.2023) = c("redo","plot","live","spp","dbh","stem.count","notes")

# clean up 2013 data
data.2013 = data.2013[-which(data.2013['dbh'] == "<50"),-which(colnames(data.2013) %in% c("stem.count"))] # remove rows with dbh <5 cm
data.2013$dbh = as.numeric(unlist(data.2013[,'dbh']))/10 # convert mm to cm
data.2013$basal.area = as.numeric(unlist(data.2013[,'basal.area']))
data.2013$year = 2013

# clean up 2022 data
data.2022 = data.2022[-which(data.2022['dbh'] == "<3"),-which(colnames(data.2022) %in% c("stem.count"))] # remove rows with dbh <3 cm
data.2022$dbh = as.numeric(unlist(data.2022[,'dbh']))
data.2022$basal.area = pi*(data.2022$dbh/100/2)^2  # m^2
data.2022$year = 2022

# remove dead tree and stem information 
data.2023 = data.2023[which(data.2023$live == "L"),]
data.2023 = data.2023[-which(data.2023$dbh == "<2.5"), -which(colnames(data.2023) %in% c("live","stem.count","notes"))]
data.2023$dbh = as.numeric(data.2023$dbh)
data.2023$basal.area = pi*(data.2023$dbh/100/2)^2
data.2023$year = 2023
data.2023 = data.2023 %>% separate(spp, into = c("genus", "species"), sep = " ", remove=F)

# combine 2013 & 2022 DBH data into one dataframe
new.cols = c("year","plot","spp","dbh","basal.area")
dbh.all = rbind(data.2013[,new.cols], data.2022[,new.cols])
dbh.all$redo = "N"

## read in allometric equation matrices
allo.df1 = read.csv("Allometric_Equations/Jenkins2004.csv", header=T); rownames(allo.df1) = allo.df1$spp
allo.df2 = read.csv("Allometric_Equations/Chojnacky2014.csv", header=T); rownames(allo.df2) = allo.df2$spp

# add columns for genus and species
uni.spp = sort(unique(dbh.all$spp))
n.sp = length(uni.spp)
dbh.all$genus = dbh.all$spp
dbh.all$species = dbh.all$spp
for (sp in uni.spp) {
  sp.ind1 = which(dbh.all$spp == sp)
  sp.ind2 = which(allo.df1$spp == sp)
  dbh.all[sp.ind1, c("genus","species")] = allo.df1[sp.ind2, c("genus","species")]
}

# add treatment and plot number columns to 2023 data and combined 2013/2022 data
dbh.all = dbh.all %>% separate(plot, into = c("trt", "num"), sep = 1, remove=F)
data.2023 = data.2023 %>% separate(plot, into = c("trt", "num"), sep = 1, remove=F)

# add columns for full treatment name
dbh.all$trt.full = rep(0, nrow(dbh.all))
data.2023$trt.full = rep(0, nrow(data.2023))
for (i in 1:6) { dbh.all$trt.full[which(dbh.all$trt == trts[i])] = trt.names[i] }
for (i in 1:6) { data.2023$trt.full[which(data.2023$trt == trts[i])] = trt.names[i] }

# combine 2023 data with 2013/2022 data
order.cols = c("redo","year","plot","trt","trt.full","num","spp","genus","species","dbh","basal.area")
dbh.all = rbind(dbh.all[,order.cols], data.2023[,order.cols])

# replace full names with abbreviations
names = c("Acer saccharinum","Fraxinus pennsylvanica","Morus alba","Platanus occidentalis","Quercus bicolor")
abbrs = c("Acesai","Fraxpen","Moralb","Plaocc","Quebic")
for (i in 1:length(names)) {
  name.id = which(dbh.all$spp == names[i])
  dbh.all$spp[name.id] = abbrs[i]
}

# write DBH data to file
write.csv(dbh.all, "Clean_Data/DBH_Data_Clean_2023.csv", row.names=F)
