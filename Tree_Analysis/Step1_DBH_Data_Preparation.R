path_to_tree_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo/Tree_Analysis"
setwd(path_to_tree_folder)

library(tidyr)

### script that loads and cleans diameter-at-breast-height (dbh.cm) data from 
### Joslin wetland mitigation site in Henry County

# define important variables
trt.letters = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")
n.t = length(trt.letters)
n.p = 3

# conversions
cm.per.m = 100

# read in data from 2022 and re-sample data from 2023
data.2022 = read.csv('Raw_Data/DBH_August2022.csv')
data.2023 = read.csv('Raw_Data/DBH_Snags_June2023.csv', header=T)

# update column names
colnames(data.2022) = c("treatment_plot","spp","dbh.cm","stem.count")
colnames(data.2023) = c("redo","treatment_plot","live","spp","dbh.cm","stem.count","notes")

# clean up 2022 data
data.2022 = data.2022[-which(data.2022['dbh.cm'] == "<3"),
                      -which(colnames(data.2022) %in% c("stem.count"))] # remove rows with dbh.cm <3 cm
data.2022$dbh.cm = as.numeric(unlist(data.2022[,'dbh.cm']))
data.2022$basal.area.m2 = pi*(data.2022$dbh.cm/cm.per.m/2)^2  
data.2022$year = 2022
data.2022$redo = "N"

# remove dead tree and stem information 
data.2023 = data.2023[which(data.2023$live == "L"),]
data.2023 = data.2023[-which(data.2023$dbh.cm == "<2.5"), 
                      -which(colnames(data.2023) %in% c("live","stem.count","notes"))]
data.2023$dbh.cm = as.numeric(data.2023$dbh.cm)
data.2023$basal.area.m2 = pi*(data.2023$dbh.cm/cm.per.m/2)^2 # m^2
data.2023$year = 2023
data.2023 = data.2023 %>% separate(spp, into = c("genus", "species"), sep = " ", remove=F)

# read in allometric equation matrices
allo.df1 = read.csv("Tree_Databases/Jenkins2004.csv", header=T); rownames(allo.df1) = allo.df1$spp
allo.df2 = read.csv("Tree_Databases/Chojnacky2014.csv", header=T); rownames(allo.df2) = allo.df2$spp

# add columns for genus and species
uni.spp = sort(unique(data.2022$spp))
n.sp = length(uni.spp)
data.2022$genus = data.2022$spp
data.2022$species = data.2022$spp
for (sp in uni.spp) {
  sp.ind1 = which(data.2022$spp == sp)
  sp.ind2 = which(allo.df1$spp == sp)
  data.2022[sp.ind1, c("genus","species")] = allo.df1[sp.ind2, c("genus","species")]
}

# split treatment_plot columns
data.2022 = data.2022 %>% separate(treatment_plot, into = c("treatment","plot"), sep = 1, remove=T)
data.2023 = data.2023 %>% separate(treatment_plot, into = c("treatment","plot"), sep = 1, remove=T)

# add columns for full treatment name
data.2022$full.treatment.name = rep(0, nrow(data.2022))
data.2023$full.treatment.name = rep(0, nrow(data.2023))
for (i in 1:n.t) { data.2022$full.treatment.name[which(data.2022$treatment == trt.letters[i])] = trt.names[i] }
for (i in 1:n.t) { data.2023$full.treatment.name[which(data.2023$treatment == trt.letters[i])] = trt.names[i] }

# combine 2022 and 2023 data
order.cols = c("year","redo","treatment","full.treatment.name","plot",
               "spp","genus","species","dbh.cm","basal.area.m2")
dbh.cm.all = rbind(data.2022[,order.cols], data.2023[,order.cols])

# replace full names with abbreviations
names = c("Acer saccharinum","Fraxinus pennsylvanica","Morus alba","Platanus occidentalis","Quercus bicolor")
abbrs = c("Acesai","Fraxpen","Moralb","Plaocc","Quebic")
for (i in 1:length(names)) {
  name.id = which(dbh.cm.all$spp == names[i])
  dbh.cm.all$spp[name.id] = abbrs[i]
}

# write dbh.cm data to file
write.csv(dbh.cm.all, "Clean_Data_By_Species/DBH_Data_Clean_2022_2023.csv", row.names=F)
