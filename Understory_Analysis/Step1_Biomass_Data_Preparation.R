path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo"
setwd(path_to_repo)

library(ggfortify)
library(patchwork)
library(reshape2)
library(dplyr)
library(tidyr)

### script that creates cleans dataframe for dry mass of litter, fine woody debris, 
### and biomass by species group

# treatments
trt.letters = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")
n.t = length(trt.letters)
n.p = 3
n.q = 5

# dimensions and conversions
subplot.length = 30 # cm
cm.per.m = 100 

# read in biomass data for 2022
bm.df = read.csv("Understory_Analysis/Raw_Data/HerbaceousBiomass_Sep2022.csv")

# remove boat mass columns
bm.df = bm.df[, -which(colnames(bm.df) %in% c("Container.mass..g.","Actual.dry.mass...container.w.o.bag..g."))]

# update column names
colnames(bm.df) = c("treatment_plot","quadrat","spp","mass")

# combine duplicate herbaceous litter and fine woody debris data within same quadrat
bm.df = bm.df %>% group_by(treatment_plot, quadrat, spp) %>% summarise(sum.mass = sum(mass))

# normalize biomass data by quadrat size (g/m^2)
bm.df$mass.g.m2 = bm.df$sum.mass/(subplot.length^2)*(cm.per.m^2)

# separate plot column into treatment and plot number 
bm.df = bm.df %>% separate(treatment_plot, into=c("treatment","plot"), sep=1, remove=T)

# add column for full treatment name
bm.df$full.treatment.name = rep(0, nrow(bm.df))
for (i in 1:n.t) {bm.df$full.treatment.name[which(bm.df$treatment == trt.letters[i])] = trt.names[i]}

# make data frame for total herbaceous biomass, herbaceous litter, and fine woody debris in each quadrat 
df.quad.sum = data.frame(matrix(nrow = n.t*n.p*n.q, ncol=10))
id.cols = c("treatment","full.treatment.name","plot","quadrat")
colnames(df.quad.sum) = c(id.cols,"Herbaceous biomass","Herbaceous litter","Fine woody debris",
                          "Mixed species","Phalaris arundinacea","Humulus japonicus")
k = 1
for (l in 1:n.t) {
  t = trt.letters[l]
  for (i in 1:n.p) {
    for (j in 1:n.q) {
      # get data for given year, treatment, plot, and quadrat
      quad.id = which((bm.df$treatment == t) & (bm.df$plot == as.character(i) & bm.df$quadrat == paste("Q", j, sep="")))
      quad.df = bm.df[quad.id,]
      
      # get HL and FWD id
      hl.id = which(quad.df$spp == "Herbaceous litter")
      fwd.id = which(quad.df$spp == "Fine woody debris")
      pa.id = which(quad.df$spp == "Phalaris arundinacea")
      hj.id = which(quad.df$spp == "Humulus japonicus")
      bm.id = which(!(quad.df$spp %in% c("Herbaceous litter","Fine woody debris")))
      m.id = which(!(quad.df$spp %in% c("Herbaceous litter","Fine woody debris",
                                        "Phalaris arundinacea","Humulus japonicus")))
                      
      # sum mass for each category
      df.quad.sum[k, "Herbaceous litter"] =  sum(quad.df[hl.id, "mass.g.m2"])
      df.quad.sum[k, "Fine woody debris"] = sum(quad.df[fwd.id, "mass.g.m2"])
      df.quad.sum[k, "Herbaceous biomass"] = sum(quad.df[bm.id, "mass.g.m2"])
      
      if (length(pa.id) > 0) { 
        df.quad.sum[k, "Phalaris arundinacea"] = sum(quad.df[pa.id, "mass.g.m2"]) 
      } else {
        df.quad.sum[k, "Phalaris arundinacea"] = 0
      }
      
      if (length(hj.id) > 0) { 
        df.quad.sum[k, "Humulus japonicus"] = sum(quad.df[hj.id, "mass.g.m2"]) 
      } else {
        df.quad.sum[k, "Humulus japonicus"] = 0
      }
      
      if (length(m.id) > 0) { 
        df.quad.sum[k, "Mixed species"] = sum(quad.df[m.id, "mass.g.m2"]) 
      } else {
        df.quad.sum[k, "Mixed species"] = 0
      }
      
      # fill in identification columns
      df.quad.sum[k, id.cols] = c(t, trt.names[l], i, paste("Q", j, sep=""))
      
      # increment counter
      k = k + 1
    }
  }
}

# write data by quadrat to file
write.csv(df.quad.sum, "Understory_Analysis/Clean_Data/Biomass_By_Quadrat_Sep2022.csv", row.names=F)
