path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo"
setwd(path_to_repo)

library(readxl)
library(vegan)
library(dplyr)
library(ggvegan)
library(ggrepel)
library(tidyr)

################################################################################
## script estimates understory species richness by plot and year

# treatment names
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

# read in species cover data from 2022 and 2023
cover2022 = read.csv("Understory_Analysis/Raw_Data/VegetationCover_Sep2022.csv")
cover2023 = read.csv("Understory_Analysis/Raw_Data/VegetationCover_Aug2023.csv")
colnames(cover2022)[which(colnames(cover2022) == "Cover class")] = "Cover"
colnames(cover2023)[which(colnames(cover2023) == "Cover class")] = "Cover"

# convert cover categories to percentages
cover.meds = c(0.5, 3, 15, 37.5, 62.5, 85, 97.5)
for (i in 1:7) {cover2022$Cover[which(cover2022$Cover == i-1)] = cover.meds[i]}
for (i in 1:7) {cover2023$Cover[which(cover2023$Cover == i-1)] = cover.meds[i]}

# combine species cover dfs
cover2022$Year = 2022; cover2023$Year = 2023
all.cov = rbind(cover2022[,-which(colnames(cover2022) %in% c("Date","Notes"))], 
                cover2023[,-which(colnames(cover2022) %in% c("Date","Notes"))])

# replace grass spp. with poaceae
all.cov$Spp[which(all.cov$Spp == "Grass spp.")] = "Poaceae spp."

# trim white space
all.cov$Spp = trimws(all.cov$Spp)

# get unique plots and quads
n.y = 2; n.t = 6; n.p = 3; n.q = 5
years = c(2022, 2023)
plots = unique(sort(all.cov$Plot))
quads = unique(sort(all.cov$Quadrat))
nums = c(1,2,3)

# remove bareground and woody debris
sp.only = all.cov[-which(all.cov$Spp %in% c("Bareground","Woody debris")),]

# split plot column
sp.only = sp.only %>% separate(Plot, c("Trt","Num"), sep = 1, remove=T)

# make dataframes for richness by quadrat and plot
n.sp.quad.df = data.frame(matrix(nrow=n.y*n.t*n.p*n.q, ncol=5))
n.sp.plot.df = data.frame(matrix(nrow=n.y*n.t*n.p, ncol=4))
colnames(n.sp.quad.df) = c("Year","Trt","Num","Quad","N")
colnames(n.sp.plot.df) = c("Year","Trt","Num","N")

# loop through each year and treatment to count species
k = 1
h = 1
for (y in 1:n.y) {
  for (t in 1:n.t) {
    for (n in 1:n.p) {
      ytn.id = which((sp.only$Year == years[y] & sp.only$Trt == trts[t]) & (sp.only$Num == nums[n]))
      ytn.df = sp.only[ytn.id,]
      ytn.sp = length(sort(unique(ytn.df$Spp)))
      n.sp.plot.df[h, c("Year","Trt","Num","N")] = c(years[y], trts[t], nums[n], ytn.sp)
      h = h + 1
      for (q in 1:n.q) {
        ytnq.id = which((sp.only$Year == years[y] & sp.only$Trt == trts[t]) & (sp.only$Num == nums[n] & sp.only$Quadrat == quads[q]))
        ytnq.df = sp.only[ytnq.id,]
        ytnq.sp = length(sort(unique(ytnq.df$Spp)))
        n.sp.quad.df[k, c("Year","Trt","Num","Quad","N")] = c(years[y], trts[t], nums[n], quads[q], ytnq.sp)
        k = k + 1
      }
    }
  }
}

# save species richness dataframe
write.csv(all.cov, "Understory_Analysis/Clean_Data/Clean_Cover_Data_2022_2023.csv", row.names=F)
write.csv(n.sp.quad.df, "Understory_Analysis/Clean_Data/Understory_Species_Richness_by_Quadrat_and_Year.csv", row.names=F)
write.csv(n.sp.plot.df, "Understory_Analysis/Clean_Data/Understory_Species_Richness_by_Plot_and_Year.csv", row.names=F)
