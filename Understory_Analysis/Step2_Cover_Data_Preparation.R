path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo"
setwd(path_to_repo)

library(dplyr)
library(tidyr)

### scripts that cleans cover data and estimates species richness for each treatment
### by plot and quadrat

# treatment names
treatment.letters = c("A","B","C","D","E","R")
treatment.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")
n.t = length(treatment.letters)

# read in species cover data from 2022 and 2023
cover.df = read.csv("Understory_Analysis/Raw_Data/VegetationCover_Sep2022.csv")
colnames(cover.df) = c("date","treatment_plot","quadrat","cover_class","spp","notes")

# convert cover categories to percentages
cover.meds = c(0.5, 3, 15, 37.5, 62.5, 85, 97.5); n.cc = length (cover.meds)
cover.df$cover = rep(0, nrow(cover.df))
for (i in 1:n.cc) { cover.df$cover[which(cover.df$cover == i-1)] = cover.meds[i] }

# trim white space
cover.df$spp = trimws(cover.df$spp)

# replace grass spp. with poaceae
cover.df$spp[which(cover.df$spp == "Grass spp.")] = "Poaceae spp."

# separate plot column into treatment and plot plotber 
cover.df = cover.df %>% separate(treatment_plot, into=c("treatment","plot"), sep=1, remove=T)

# get unique plots and quads
plots = unique(sort(cover.df$plot))
quads = unique(sort(cover.df$quadrat))
n.p = length(plots); n.q = length(quads)

# remove bareground and woody debris
sp.only = cover.df[-which(cover.df$spp %in% c("Bareground","Woody debris")),]

# make dataframes for richness by quadrat and plot
n.sp.quad.df = data.frame(matrix(nrow=n.t*n.p*n.q, ncol=4))
n.sp.plot.df = data.frame(matrix(nrow=n.t*n.p, ncol=3))
colnames(n.sp.quad.df) = c("treatment","plot","quadrat","n")
colnames(n.sp.plot.df) = c("treatment","plot","n")

# loop through each treatment and plot to count species
k = 1
h = 1
for (t in 1:n.t) {
  for (n in 1:n.p) {
    tn.id = which((sp.only$treatment == trt.letters[t]) & (sp.only$plot == plots[n]))
    tn.df = sp.only[tn.id,]
    tn.sp = length(sort(unique(tn.df$spp)))
    n.sp.plot.df[h, c("treatment","plot","n")] = c(trt.letters[t], plots[n], tn.sp)
    h = h + 1
    for (q in 1:n.q) {
      tnq.id = which((sp.only$treatment == trt.letters[t]) & (sp.only$plot == plots[n] & sp.only$quadrat == quads[q]))
      tnq.df = sp.only[tnq.id,]
      tnq.sp = length(sort(unique(tnq.df$spp)))
      n.sp.quad.df[k, c("treatment","plot","quadrat","n")] = c(trt.letters[t], plots[n], quads[q], tnq.sp)
      k = k + 1
    }
  }
}

# save species richness dataframe
write.csv(cover.df, "Understory_Analysis/Clean_Data/Clean_Cover_Data_Sep2022.csv", row.names=F)
write.csv(n.sp.quad.df, "Understory_Analysis/Clean_Data/Understory_Species_Richness_by_Quadrat_Sep2022.csv", row.names=F)
write.csv(n.sp.plot.df, "Understory_Analysis/Clean_Data/Understory_Species_Richness_by_Plot_Sep2022.csv", row.names=F)
