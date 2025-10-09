path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo"
setwd(path_to_repo)

library(dplyr)

# treatment labels
trt.letters = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")
n.t = length(trt.letters)
plots = seq(1,3)
n.p = length(plots)

################################################################################
## read in all C stock estimates

# read in woody biomass/debris and understory C stock data
c.data = read.csv("Tree_Analysis/Clean_Data_By_Plot/All_Vegetation_C_Stocks_By_Plot.csv")

# read in soil data
soil.data = read.csv("Soil_Analysis/Clean_Data/Soil_Data_by_Quadrat_June2023.csv", header=T)

# estimate SOC stocks (Mg/ha)
soil.data$tic.stock = soil.data$tic.percent * soil.data$bulk.density * 30 # g/cm2 = Mg/ha
soil.data$moc.stock = soil.data$moc.percent * soil.data$bulk.density * 30 # g/cm2 = Mg/ha
soil.data$poc.stock = soil.data$poc.percent * soil.data$bulk.density * 30 # g/cm2 = Mg/ha
soil.data$soc.stock = soil.data$toc.percent * soil.data$bulk.density * 30 # g/cm2 = Mg/ha
soil.data$tc.stock = soil.data$bulk.c.percent * soil.data$bulk.density * 30 # g/cm2 = Mg/ha

# average results across plots
soil.aves = soil.data %>% 
            group_by(full.treatment.name, treatment, plot) %>%
            summarize(moc.stock = mean(moc.stock),
                      poc.stock = mean(poc.stock),
                      soc.stock = mean(soc.stock),
                      tic.stock = mean(tic.stock),
                      tc.stock = mean(tc.stock))

# update column names
colnames(soil.aves)[4:8] = c("MOC","POC","SOC","TIC","TC")

# join vegetation and soil c stocks
c.data = right_join(c.data, 
                    soil.aves, 
                    by=c("treatment","full.treatment.name","plot"))

# estimate aggregate C stocks in different pools
c.data = c.data %>%
         mutate(total.abg.wood.carbon = abg.live.stem.carbon + aboveground.woody.c.stock,
                total.bg.wood.carbon = bg.live.stem.carbon + belowground.woody.c.stock,
                total.live.carbon = total.woody.c.stock + tot.live.stem.carbon + herbaceous.biomass.c.stock,
                total.dead.carbon = snag.carbonmin + dead.stem.carbonmin + cwd.carbon + int.fwd.carbon + fine.woody.debris.c.stock + herbaceous.litter.c.stock)
c.data = c.data %>%
         mutate(total.veg.carbon = total.live.carbon + total.dead.carbon,
                total.ecosystem.carbon = TC + total.veg.carbon)

################################################################################
## estimating richness by in tree and understory layer

# read in tree species data
tree.C.df = read.csv("Tree_Analysis/Clean_Data_By_Species/WoodyBiomass_C_Stocks_By_Species.csv", header=T)

# combine genus and specific epithet
tree.C.df$spp = paste(tree.C.df$genus, tree.C.df$species, sep=" ")

# read in species cover data
cover.df = read.csv("Understory_Analysis/Clean_Data/Clean_Cover_Data_Sep2022.csv")
colnames(cover.df) = tolower(colnames(cover.df))

# make dataframe for total unique species
total.sp.df = data.frame(matrix(nrow=n.t*n.p, ncol=5))
colnames(total.sp.df) = c("treatment","plot","n.herb","n.tree","n.total")
k = 1
for (t in 1:n.t) {
  for (n in 1:n.p) {
    tree.id = which(tree.C.df$treatment == trt.letters[t] & tree.C.df$plot == plots[n])
    tree.df = tree.C.df[tree.id,]
    tree.sp = sort(unique(tree.df$spp))
    n.tree.sp = length(tree.sp)
    
    herb.id = which(cover.df$treatment == trt.letters[t] & cover.df$plot == plots[n])
    herb.df = cover.df[herb.id,]
    herb.sp = sort(unique(herb.df$spp))
    n.herb.sp = length(herb.sp)
    
    n.total.sp = length(unique(c(tree.sp,herb.sp)))
    total.sp.df[k, c("treatment","plot","n.herb","n.tree","n.total")] = c(trt.letters[t], 
                                                                          as.integer(plots[n]),
                                                                          as.integer(n.herb.sp),
                                                                          as.integer(n.tree.sp),
                                                                          as.integer(n.total.sp))
    k = k + 1
  }
}

# add richness values to biomass data frame
c.data$plot = as.numeric(c.data$plot)
total.sp.df$plot = as.numeric(total.sp.df$plot)
total.sp.df$n.herb = as.integer(total.sp.df$n.herb)
total.sp.df$n.tree = as.integer(total.sp.df$n.tree)
total.sp.df$n.total = as.integer(total.sp.df$n.total)
c.sp.data = right_join(c.data, total.sp.df, by=c("treatment","plot"))

# write data to file
write.csv(c.sp.data, "Tree_Analysis/Clean_Data_By_Plot/Clean_Veg_Soil_C_Stocks_Richness_by_Plot.csv", row.names=F)