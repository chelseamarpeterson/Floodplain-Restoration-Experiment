setwd("C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/VegetationData")

library(readxl)
library(vegan)
library(dplyr)
library(ggvegan)
library(ggrepel)
library(tidyr)

## script that analyzes vegetation cover data

# treatment names
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

# read in species cover data from 2022 and 2023
cover2022 = read_excel("VegetationCover_Sep2022.xlsx")
cover2023 = read_excel("VegetationCover_Aug2023.xlsx")
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

# make column for species abbreviations
all.cov$Spp.abr = all.cov$Spp
all.cov = all.cov %>% separate(Spp, c('Genus', 'Species'), remove=F)
for (i in 1:nrow(all.cov)) {
  if (all.cov$Genus[i] == "Bareground") {
    all.cov[i, "Spp.abr"] = "Bare"
  } else if (all.cov$Genus[i] == "Woody") {
    all.cov[i, "Spp.abr"] = "Debris"
  } else {
    all.cov[i, "Spp.abr"] = paste(substr(all.cov$Genus[i], 1, 3), substr(all.cov$Species[i], 1, 3), sep=".")
  }
}

# get all unique species
uni.sp = sort(unique(all.cov$Spp.abr))
n.sp = length(uni.sp)

# get unique plots and quads
n.y = 2; n.t = 6; n.p = 3; n.q = 5
years = c(2022, 2023)
plots = unique(sort(all.cov$Plot))
quads = unique(sort(all.cov$Quadrat))

################################################################################
# estimate the number of unique species by plot and quadrat

# remove bareground and woody debris
sp.only = all.cov[-which(all.cov$Spp == "Bareground"),]
sp.only = sp.only[-which(sp.only$Spp == "Woody debris"),]

# split plot column
sp.only = sp.only %>% separate(Plot, c("Trt","Num"), sep = 1, remove=T)

# make dataframes
n.sp.quad.df = data.frame(matrix(nrow=n.y*n.t*n.p*n.q, ncol=5))
n.sp.plot.df = data.frame(matrix(nrow=n.y*n.t*n.p, ncol=4))
colnames(n.sp.quad.df) = c("Year","Trt","Num","Quad","N")
colnames(n.sp.plot.df) = c("Year","Trt","Num","N")

# loop through each year and treatment to count species
nums = c(1,2,3)
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
write.csv(n.sp.quad.df, "Species_Richness_by_Quadrat_and_Year.csv", row.names=F)
write.csv(n.sp.plot.df, "Species_Richness_by_Plot_and_Year.csv", row.names=F)

################################################################################
## canonical correlation analysies

# make community data matrix
df.sp = data.frame(matrix(nrow=n.y*n.t*n.p*n.q, ncol=3+n.sp)) 
colnames(df.sp) = c("Year","Plot","Quad", uni.sp)
k = 1
for (y in 2022:2023) {
  for (p in plots) {
    for (q in quads) {
      # fill out year, plot, quad columns
      df.sp[k, c("Year","Plot","Quad")] = c(y, p, q)
      
      # get subset of dataframe
      ypq.id = which((all.cov$Year == y & all.cov$Plot == p) & all.cov$Quadrat == q)
      ypq.sp = all.cov[ypq.id, c("Spp.abr","Cover")]
      
      # get col ids for species 
      sp.id = which(colnames(df.sp) %in% ypq.sp$Spp.abr)
      df.sp[k, sp.id] = ypq.sp$Cover
      
      # assign remainder of columns zero
      df.sp[k, -c(seq(1,3), sp.id)] = 0
      
      # increment counter
      k = k + 1
    }
  }
}
df.sp = df.sp %>% separate(Plot, c("Trt","Num"), sep = 1, remove=F)

# separate species data into years
df.sp.2022 = data.frame(df.sp[which(df.sp$Year == 2022),-which(colnames(df.sp) %in% c("Trt","Plot","Num","Quad","Year"))])
df.sp.2023 = data.frame(df.sp[which(df.sp$Year == 2023),-which(colnames(df.sp) %in% c("Trt","Plot","Num","Quad","Year"))])
row.names(df.sp.2022) = paste(df.sp$Plot[1:90], df.sp$Quad[1:90])
row.names(df.sp.2023) = paste(df.sp$Plot[1:90], df.sp$Quad[1:90])

# count number of plots where each species occurs
quad.count.2022 = sort(apply(df.sp.2022, 2, function(x) sum(x != 0)), decreasing=F)
quad.count.2023 = sort(apply(df.sp.2023, 2, function(x) sum(x != 0)), decreasing=F)
df.count.2022 = data.frame(quad.count.2022)
df.count.2023 = data.frame(quad.count.2023)
colnames(df.count.2022) = "count"
colnames(df.count.2023) = "count"
df.count.2022$spp = factor(row.names(df.count.2022), levels=row.names(df.count.2022))
df.count.2023$spp = factor(row.names(df.count.2023), levels=row.names(df.count.2023))
ggplot(df.count.2023[-which(df.count.2023$count == 0),], aes(y=spp, x=count)) + geom_bar(stat="identity")

# remove species that occur in 1 or fewer quadrats
cols.2022 = rownames(df.count.2022)[-which(df.count.2022$count <= 1)]
cols.2023 = rownames(df.count.2023)[-which(df.count.2023$count <= 1)]
df.sp.2022 = df.sp.2022[,cols.2022]
df.sp.2023 = df.sp.2023[,cols.2023]


## read in soil data matrix
setwd("C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/SoilData")
soil.data = read.csv("All_Soil_Data_2023.csv", header=T)
soil.scale = soil.data[,-which(colnames(soil.data) %in% c("Trt.full"))]

# sort soil data
soil.sort = sort(soil.scale$Plot, index.return=T)
soil.scale = soil.scale[soil.sort$ix,]
soil.scale[,5:ncol(soil.scale)] = scale(soil.scale[,-seq(1,4)], center=T, scale=T)

## run CCA
cca.2022 = cca(df.sp.2022~Temperature+Gravimetric.Moisture+Bulk.Density+Sand+Silt+pH+P+K+Ca+Mg+SOM+NO3+NH4+CEC, data=soil.scale)
cca.2023 = cca(df.sp.2023~Temperature+Gravimetric.Moisture+Bulk.Density+Sand+Silt+pH+P+K+Ca+Mg+SOM+NO3+NH4+CEC, data=soil.scale)
summary(cca.2022)
summary(cca.2023)
plot(cca.2022, display=c("sp","cn"))
plot(cca.2023, display=c("sp","cn"))

## make custom plot

# species dataframe
cca.df.sp = data.frame(cca.2023$CCA$v)

# scale by eigenvalues
eigs = cca.2023$CCA$eig
cca.df.sp.scale = cca.df.sp
for (i in 1:ncol(cca.df.sp)) { cca.df.sp.scale[,i] = sqrt(eigs[i])*cca.df.sp.scale[,i]}
cca.df.sp.scale$sp = row.names(cca.df.sp)

# sites dataframe
cca.df.sites = data.frame(cca.2023$CCA$wa)
cca.df.sites$trt = soil.scale$Trt
cca.df.sites$trt.full = soil.scale$Trt
for (i in 1:6) { cca.df.sites$trt.full[which(cca.df.sites$trt == trts[i])] = trt.names[i]}

# variable dataframe
cca.df.vars = data.frame(cca.2023$CCA$biplot)
cca.df.vars$var = row.names(cca.2023$CCA$biplot)

# plot everything
ggplot(NULL, aes(x=CCA1, y=CCA2)) + labs(title="2023") +
       geom_point(data=cca.df.sites, aes(color=factor(trt.full, levels=trt.names)), size=2.5) +
       geom_point(data=cca.df.sp.scale) + 
       geom_text(data=cca.df.sp.scale, aes(label=sp), size=2) +
       geom_segment(data=cca.df.vars, aes(x=0, y=0, xend=3*CCA1, yend=3*CCA2),
                    arrow=arrow(length=unit(0.2, "cm"), type="open"), color="red") +
       geom_text(data=cca.df.vars, aes(x=3*CCA1, y=3*CCA2, label=var),
                 color="red", size=3) + theme(legend.title=element_blank())

ggplot(NULL, aes(x=CCA1, y=CCA3)) + labs(title="2023") +
      geom_point(data=cca.df.sites, aes(color=factor(trt.full, levels=trt.names)), size=2.5) +
      geom_point(data=cca.df.sp.scale) + 
      geom_text(data=cca.df.sp.scale, aes(label=sp), size=2) +
      geom_segment(data=cca.df.vars, aes(x=0, y=0, xend=3*CCA1, yend=3*CCA3),
                   arrow=arrow(length=unit(0.2, "cm"), type="open"), color="red") +
      geom_text(data=cca.df.vars, aes(x=3*CCA1, y=3*CCA3, label=var),
                color="red", size=3) + theme(legend.title=element_blank())


# species dataframe
cca.df.sp = data.frame(cca.2022$CCA$v)

# scale by eigenvalues
eigs = cca.2022$CCA$eig
cca.df.sp.scale = cca.df.sp
for (i in 1:ncol(cca.df.sp)) { cca.df.sp.scale[,i] = sqrt(eigs[i])*cca.df.sp.scale[,i]}
cca.df.sp.scale$sp = row.names(cca.df.sp)

# sites dataframe
cca.df.sites = data.frame(cca.2022$CCA$wa)
cca.df.sites$trt = soil.scale$Trt
cca.df.sites$trt.full = soil.scale$Trt
for (i in 1:6) { cca.df.sites$trt.full[which(cca.df.sites$trt == trts[i])] = trt.names[i]}

# variable dataframe
cca.df.vars = data.frame(cca.2022$CCA$biplot)
cca.df.vars$var = row.names(cca.2022$CCA$biplot)

# plot everything
ggplot(NULL, aes(x=CCA1, y=CCA2)) + labs(title="2022") +
      geom_point(data=cca.df.sites, aes(color=factor(trt.full, levels=trt.names)), size=2.5) +
      geom_point(data=cca.df.sp.scale) + 
      geom_text(data=cca.df.sp.scale, aes(label=sp), size=2) +
      geom_segment(data=cca.df.vars, aes(x=0, y=0, xend=3*CCA1, yend=3*CCA2),
                   arrow=arrow(length=unit(0.2, "cm"), type="open"), color="red") +
      geom_text(data=cca.df.vars, aes(x=3*CCA1, y=3*CCA2, label=var),
                color="red", size=3) + theme(legend.title=element_blank())


# test significance
anova.cca(cca.2022)
anova.cca(cca.2023)
anova.cca(cca.2022, by="term")
anova.cca(cca.2023, by="term")
anova.cca(cca.2023, by="margin")


