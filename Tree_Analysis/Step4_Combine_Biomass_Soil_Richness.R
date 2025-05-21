path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Public-Repo"
setwd(path_to_repo)

# treatment labels
trts = c("A","B","C","D","E","R")

## read in woody biomass/debris and understory C stock data
c.data = read.csv("Tree_Analysis/Clean_Data/All_Vegetation_C_Stocks_By_Plot.csv")

## read in soil data
soil.data = read.csv("Soil_Analysis/Clean_Data/Soil_Data_Averaged_by_Plot_June2023.csv", header=T)

# update soil data column names
soil.c.data = soil.data[,c("trt","trt.full","num","maoc.stock","poc.stock","soc.stock","tic.stock","tc.stock")]
colnames(soil.c.data)[4:8] = c("MAOM-C","POM-C","SOC","TIC","TC")

## join vegetation and soil c stocks
c.data = right_join(c.data, 
                    soil.c.data[,c("trt","trt.full","num","MAOM-C","POM-C","SOC","TIC","TC")], 
                    by=c("trt","trt.full","num"))

## read in species richness data

# herbaceous species data
sp.richness.df = read.csv("Understory_Data/Species_Richness_by_Plot_and_Year.csv")

# average richness for 2022 to be consistent with biomass data
sp.aves = sp.richness.df[which(sp.richness.df$Year == 2022),] %>% 
          group_by(Trt, Num) %>% 
          summarize(N = mean(N))
colnames(sp.aves) = c("trt","num","N.herb")

## read in tree species data
C.stock.sp.df = read.csv("Tree_Analysis/Clean_Data/WoodyBiomass_C_Stocks_By_Species.csv", header=T)

# get data for 2022
C.stock.sp.df.2022 = C.stock.sp.df[which(C.stock.sp.df$year == 2022),]

# make dataframe for unique tree species
n.t = 6; n.p = 3; nums = seq(1,3)
n.tree.sp.df = data.frame(matrix(nrow=n.t*n.p, ncol=3))
colnames(n.tree.sp.df) = c("trt","num","N.tree")
k = 1
for (t in 1:n.t) {
  for (n in 1:n.p) {
    tn.id = which(C.stock.sp.df.2022$trt == trts[t] & C.stock.sp.df.2022$num == nums[n])
    tn.df = C.stock.sp.df.2022[tn.id,]
    tn.sp = length(sort(unique(tn.df$spp)))
    n.tree.sp.df[k, c("trt","num","N.tree")] = c(trts[t], as.numeric(nums[n]), as.numeric(tn.sp))
    k = k + 1
  }
}

# add richness values to biomass data frame
c.data$num = as.numeric(c.data$num)
n.tree.sp.df$num = as.numeric(n.tree.sp.df$num)
n.tree.sp.df$N.tree = as.numeric(n.tree.sp.df$N.tree)
c.data = right_join(c.data, sp.aves, by=c("trt","num"))
c.data = right_join(c.data, n.tree.sp.df, by=c("trt","num"))


# write data to file
write.csv(c.data, "Tree_Analysis/Clean_Data/Clean_Veg_Soil_Cstocks_Richness.csv", row.names=F)