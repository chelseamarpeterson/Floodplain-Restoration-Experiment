path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo"
setwd(path_to_repo)

# treatment labels
trts = c("A","B","C","D","E","R")

## read in woody biomass/debris and understory C stock data
c.data = read.csv("Tree_Analysis/Clean_Data/All_Vegetation_C_Stocks_By_Plot.csv")

## read in soil data
soil.data = read.csv("Soil_Analysis/Clean_Data/Soil_Data_Averaged_by_Plot_June2023.csv", header=T)

# update soil data column names
soil.c.data = soil.data[,c("trt","trt.full","num","maoc.stock","poc.stock","soc.stock","tic.stock","tc.stock")]
colnames(soil.c.data)[4:8] = c("MAOC","POC","SOC","TIC","TC")

## join vegetation and soil c stocks
c.data = right_join(c.data, 
                    soil.c.data[,c("trt","trt.full","num","MAOC","POC","SOC","TIC","TC")], 
                    by=c("trt","trt.full","num"))

## read in tree species data
tree.C.df = read.csv("Tree_Analysis/Clean_Data/WoodyBiomass_C_Stocks_By_Species.csv", header=T)

# combine genus and specific epithet
tree.C.df$spp = paste(tree.C.df$genus, tree.C.df$species, sep=" ")

# get data for 2022
tree.C.2022 = tree.C.df[which(tree.C.df$year == 2022),]

## read in species cover data
all.cov.df = read.csv("Understory_Analysis/Clean_Data/Clean_Cover_Data_2022_2023.csv")
colnames(all.cov.df) = tolower(colnames(all.cov.df))

# 2022 cover deta
cov.2022 = all.cov.df[which(all.cov.df$year == 2022),]

# split plot column
cov.2022 = cov.2022 %>% separate(plot, c("trt","num"), sep=1, remove=F)
cov.2022$num = as.numeric(cov.2022$num)

# make dataframe for total unique species
n.t = 6; n.p = 3; nums = seq(1,3)
total.sp.2022 = data.frame(matrix(nrow=n.t*n.p, ncol=5))
colnames(total.sp.2022) = c("trt","num","n.herb","n.tree","n.total")
k = 1
for (t in 1:n.t) {
  for (n in 1:n.p) {
    tree.id = which(tree.C.2022$trt == trts[t] & tree.C.2022$num == nums[n])
    tree.df = tree.C.2022[tree.id,]
    tree.sp = sort(unique(tree.df$spp))
    n.tree.sp = length(tree.sp)
    
    herb.id = which(cov.2022$trt == trts[t] & cov.2022$num == nums[n])
    herb.df = cov.2022[herb.id,]
    herb.sp = sort(unique(herb.df$spp))
    n.herb.sp = length(herb.sp)
    
    n.total.sp = length(unique(c(tree.sp,herb.sp)))
    total.sp.2022[k, c("trt","num","n.herb","n.tree","n.total")] = c(trts[t], 
                                                                     as.numeric(nums[n]),
                                                                     as.numeric(n.herb.sp),
                                                                     as.numeric(n.tree.sp),
                                                                     as.numeric(n.total.sp))
    k = k + 1
  }
}

# add richness values to biomass data frame
c.data$num = as.numeric(c.data$num)
total.sp.2022$num = as.numeric(total.sp.2022$num)
c.data = right_join(c.data, total.sp.2022, by=c("trt","num"))


# write data to file
write.csv(c.data, "Tree_Analysis/Clean_Data/Clean_Veg_Soil_Cstocks_Richness.csv", row.names=F)