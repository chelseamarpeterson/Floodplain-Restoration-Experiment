setwd('C:/Users/Chels/OneDrive - University of Illinois - Urbana/PhD/VegetationData')


## script that combines 2013 & 2022 vegetation data and adds database information

# read in species database
#sp.db = read_excel('SpeciesDatabase.xlsx')
#WFO.data = read.csv("WFO_Backbone/classification.csv", sep="\t")

# read in cover data 
cover.2013 = read_excel('AllData_2013.xlsx', sheet="quadrat data")
cover.2022 = read_excel('VegetationCover_Sep2022.xlsx')
colnames(cover.2013) = c("plot","quad","cover","spp")
colnames(cover.2022) = c("date","plot","quad","cover","spp","notes")
cover.2013$cover = as.numeric(cover.2013$cover)

# remove litter, woody debris, and bareground rows
cover.2013 = cover.2013[!(cover.2013$spp %in% c("BARE")),]
cover.2022 = cover.2022[!(cover.2022$spp %in% c("Litter","Woody debris","Bareground")),]

# update Carex rows 
cover.2013[which(cover.2013$spp == "Carex sp."),"spp"] = "Carex grisea"
cover.2022[which(cover.2022$spp == "Carex sp."),"spp"] = "Carex grayi"

# change Aster lanceolatus in 2013 data to Symphyotrichum lanceolatum
cover.2013$spp[which(cover.2013$spp == "Aster lanceolatus var. simplex")] = "Symphyotrichum lanceolatum"

# convert cover categories to percentages
cover.meds = c(2.5,15,37.5,62.5,85,97.5)
for (i in 1:6) {cover.2022$cover[which(cover.2022$cover == i)] = cover.meds[i]}

# combine 2013 & 2022 data
cover.2013$year = 2013; cover.2022$year = 2022
all.cov = rbind(cover.2013,cover.2022[,-which(colnames(cover.2022) %in% c("date","notes"))])

# split plot column into treatment and year & quad column into treatment and year
all.cov = all.cov %>% separate(plot, into = c("trt", "plot"), sep = "(?<=[A-Za-z])(?=[0-9])")
all.cov = all.cov %>% separate(quad, into = c("quad.letter","quad"), sep = "(?<=[A-Za-z])(?=[0-9])")

# add C-values, native indicator, wetland indicator, growth form column to cover data
uni.sp = sort(unique(all.cov$spp))
n.pts = length(all.cov$spp)
my.cols = c("MW_family","wis","growth.form","c.value","native","life.dur")
db.cols = c("FAMILY","MW_2012","PHYSIOGNOMY","CC_NATIVE","NATIVE","PERENNIAL")
df = data.frame((matrix(nrow=n.pts, ncol=6)))
colnames(df) = my.cols
all.cov = cbind(all.cov,df)
for (sp in uni.sp) {
  sp.ind1 = which(all.cov$spp == sp)
  sp.ind2 = which(sp.db$PSID == sp | sp.db$Synonym == sp)[1] 
  for (i in 1:3) {all.cov[sp.ind1,my.cols[i]] = sp.db[sp.ind2,db.cols[i]]}
  for (i in 4:6) {all.cov[sp.ind1,my.cols[i]] = as.numeric(sp.db[sp.ind2,db.cols[i]])}
}

# add family and order to cover data
n.pts = length(all.cov$spp)
df.tax = data.frame(matrix(nrow=n.pts, ncol=6))
colnames(df.tax) = c("family","order","node1","node2","node3","node4")
for (i in 1:n.pts) {
  sp = all.cov$spp[i]
  wfo.data = WFO.family(sp, WFO.data=WFO.data, verbose=F)
  print(wfo.data)
  wfo.fam = wfo.data$Family
  wfo.order = wfo.data$Order
  wfo.node1 = wfo.data$Node.1
  wfo.node2 = wfo.data$Node.2
  wfo.node3 = wfo.data$Node.3
  wfo.node4 = wfo.data$Node.4
  if (length(wfo.fam) > 0) {df.tax$family[i] = wfo.fam}
  if (length(wfo.order) > 0) {df.tax$order[i] = wfo.order}
  if (length(wfo.node1) > 0) {df.tax$node1[i] = wfo.node1}
  if (length(wfo.node2) > 0) {df.tax$node2[i] = wfo.node2}
  if (length(wfo.node3) > 0) {df.tax$node3[i] = wfo.node3}
  if (length(wfo.node4) > 0) {df.tax$node4[i] = wfo.node4}
}
all.cov = cbind(all.cov,df.tax)

# check for NAs
apply(all.cov, 2, FUN= function(x) sum(is.na(x)))
#na.ind = which(is.na(all.cov[,c("wis")])); all.cov[na.ind,]

# add numeric column for wetland indicator
all.cov$wis.num = rep(0, length(all.cov$spp))
wis_names = c("OBL","FACW","FAC","FACU","UPL")
for (i in 1:5) {
  wis = wis_names[i]
  wis.id = which(all.cov$wis == wis)
  all.cov[wis.id,"wis.num"] = i
} 

# save cover data with family information
write.csv(all.cov,"PlantCoverClean_2022.csv", row.names=F)
