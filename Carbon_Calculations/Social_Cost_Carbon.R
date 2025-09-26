path_to_repo= "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo"
setwd(path_to_repo)

library(ggplot2)
library(reshape2)
library(patchwork)

# read in social cost of carbon estimates 
scc.df = read.csv("Carbon_Calculations/SCC_Estimates.csv")
colnames(scc.df) = c("stat","per.CO2","per.CO2.C")
scc.df$stat = c("mean","lower","upper")
rownames(scc.df) = c("mean","lower","upper")

# read in establishment costs
est.df = read.csv("Carbon_Calculations/Treatment_Establishment_Costs.csv")
colnames(est.df) = c("trt","cost.2019","cost.2025")

# read in ecosystem carbon estimates
stock.df1 = read.csv("Tree_Analysis/Posteriors/CarbonRichness_Means_Intervals_5Chains.csv")
stock.df2 = read.csv("Tree_Analysis/Clean_Data/Clean_Veg_Soil_Cstocks_Richness.csv")
sort(unique(stock.df1$var))
ecoC.df1 = stock.df1[which(stock.df1$var == "Total ecosystem"),c("trt","mean","X5","X95")]
ecoC.df2 = stock.df2[,c("trt.full","num","total.ecosystem.carbon")]
colnames(ecoC.df1) = c("trt","mean","lower","upper")
colnames(ecoC.df2) = c("trt","num","stock")

# melt ecosystem carbon estimates by statistic
ecoC.df1.melt = melt(ecoC.df1, id.vars=c("trt"), 
                     variable.name="stat",value.name="stock")

# estimate carbon benefit of each treatment
ecoC.df1.melt$carbon.benefit = rep(0, nrow(ecoC.df1.melt))
ecoC.df2.melt = data.frame(matrix(nrow=0,ncol=4))
colnames(ecoC.df2.melt) = c("trt","num","stat","stock")
stats = c("mean","lower","upper")
for (i in 1:3) {
  stat.i = stats[i]
  stat.id = which(ecoC.df1.melt$stat == stat.i)
  ecoC.df1.melt$carbon.benefit[stat.id] = ecoC.df1.melt$stock[stat.id]*scc.df[stat.i,"per.CO2.C"]
  ecoC.df2$stat = stat.i
  ecoC.df2$carbon.benefit = ecoC.df2$stock * scc.df[stat.i,"per.CO2.C"] 
  ecoC.df2.melt = rbind(ecoC.df2.melt, ecoC.df2)
}

# treatment names
trts = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

# estimate net benefit of each treatment by subtracting the establishment cost
ecoC.df1.melt$net.benefit = rep(0, nrow(ecoC.df1.melt))
ecoC.df1.melt$breakeven.scc = rep(0, nrow(ecoC.df1.melt))
ecoC.df2.melt$net.benefit = rep(0, nrow(ecoC.df2.melt))
for (i in 1:6) {
  trt.i = trts[i]
  cost.i = est.df[which(est.df$trt == trt.i),"cost.2025"]
  trt.id1 = which(ecoC.df1.melt$trt == trt.i)
  trt.id2 = which(ecoC.df2.melt$trt == trt.i)
  ecoC.df1.melt[trt.id1,"net.benefit"] = ecoC.df1.melt[trt.id1,"carbon.benefit"] - cost.i
  ecoC.df1.melt[trt.id1,"breakeven.scc"] = cost.i/ecoC.df1.melt[trt.id1,"stock"]* 12.01/44.01
  ecoC.df2.melt[trt.id2,"net.benefit"] = ecoC.df2.melt[trt.id2,"carbon.benefit"] - cost.i
}

# get species richness estimates
n.df1 = stock.df1[which(stock.df1$var == "Total richness"),c("trt","mean","X5","X95")]
n.df2 = stock.df2[,c("trt.full","num","n.total")]
colnames(n.df1) = c("trt","mean","lower","upper")
colnames(n.df2) = c("trt","num","richness")
n.df.melt1 = melt(n.df1, id.vars=c("trt"), variable.name="stat",value.name="richness")

# join richness estimates 
ecoC.n.df.join1 = left_join(ecoC.df1.melt, n.df.melt1, by=c("trt","stat"))
ecoC.n.df.join2 = left_join(ecoC.df2.melt, n.df2, by=c("trt","num"))

# plot net carbon benefit v. richness'
mean.df1 = ecoC.n.df.join1[which(ecoC.n.df.join1$stat == "mean"),]
mean.df2 = ecoC.n.df.join2[which(ecoC.n.df.join2$stat == "mean"),]
scc.df$stat = c("Mean","Lower (5%)","Upper (95%)")
colnames(scc.df)[1] = "Social Cost of Carbon Estimate"
p1 = ggplot(data=mean.df1, 
            aes(y=richness, x=stock, 
                color=factor(trt, levels=trts))) +
            geom_point(size=4) + 
            guides(color="none") + 
            scale_y_continuous(breaks=seq(7,21,by=2),limits=c(7,21)) +
            #scale_x_continuous(breaks=seq(100,200,by=25),
            #                   limits=c(90,200)) +
            labs(x="Ecosystem carbon stock (Mg/ha)",
                 y="Total species richness") +
            theme(text=element_text(size=14))
p2 = ggplot(data=mean.df1, 
            aes(y=richness, x=carbon.benefit/1000, 
                color=factor(trt, levels=trts))) +
            geom_point(size=4) + 
            guides(color="none") + 
            scale_y_continuous(breaks=seq(7,21,by=2),limits=c(7,21)) +
            scale_x_continuous(breaks=seq(100,200,by=25),
                               limits=c(90,200)) +
            labs(x="Total carbon benefit ($1,000/ha)",
                 y="Total species richness") +
            theme(text=element_text(size=14))
p3 = ggplot(data=mean.df1, 
            aes(y=richness, x=net.benefit/1000, 
                color=factor(trt, levels=trts))) +
            geom_point(size=4) +
            guides(color="none") + 
            scale_y_continuous(breaks=seq(7,21,by=2),limits=c(7,21)) +
            labs(x="Net economic benefit ($1,000/ha)",
                 y="") +
            theme(axis.text.y=element_blank(),
                  text=element_text(size=14))
p4 = ggplot() +
     geom_point(data=mean.df1, 
                aes(y=richness, x=breakeven.scc, 
                    color=factor(trt, levels=trts)),
                size=4) +
     guides(color=guide_legend(title="Treatment")) + 
     scale_y_continuous(breaks=seq(7,21,by=2),limits=c(7,21)) +
     labs(x="Breakeven social cost of carbon ($/Mg CO2)",
          y="") +
     theme(axis.text.y=element_blank(),
           text=element_text(size=14)) +
     geom_vline(data=scc.df, 
                 aes(xintercept=per.CO2, 
                     linetype=`Social Cost of Carbon Estimate`))
(p1 + p2)/(p3 + p4)

# print results for table 
cols = colnames(ecoC.n.df.join1[3:7])
for (i in 1:length(cols)) {
  mean.df = ecoC.n.df.join1[which(ecoC.n.df.join1$stat == "mean"),c("trt",cols[i])]
  l.df = ecoC.n.df.join1[which(ecoC.n.df.join1$stat == "lower"),c("trt",cols[i])]
  u.df = ecoC.n.df.join1[which(ecoC.n.df.join1$stat == "upper"),c("trt",cols[i])]
  print(cols[i])
  for (j in 1:6) {
    trt = mean.df$trt[j]
    if (cols[i] == "stock") {
      print(paste(trt, ": ", 
                  round(mean.df[j,cols[i]]), " (", 
                  round(l.df[j,cols[i]]), "-", 
                  round(u.df[j,cols[i]]), ")", sep="")) 
      print(paste(trt, ": ", 
                  round(mean.df[j,cols[i]],1), " (", 
                  round(l.df[j,cols[i]],1), "-", 
                  round(u.df[j,cols[i]],1), ")", sep="")) 
    } else if (cols[i] == "richness") {
      print(paste(trt, ": ", 
                  round(mean.df[j,cols[i]],1), " (", 
                  round(l.df[j,cols[i]],1), "-", 
                  round(u.df[j,cols[i]],1), ")", sep=""))      
      print(paste(trt, ": ", 
                  round(mean.df[j,cols[i]],2), " (", 
                  round(l.df[j,cols[i]],2), "-", 
                  round(u.df[j,cols[i]],2), ")", sep=""))
    } else if (cols[i] %in% c("carbon.benefit","net.benefit")) {
      print(paste(trt, ": ", 
                  round(mean.df[j,cols[i]]/1000), " (", 
                  round(l.df[j,cols[i]]/1000), "-", 
                  round(u.df[j,cols[i]]/1000), ")", sep=""))
      print(paste(trt, ": ", 
                  round(mean.df[j,cols[i]]/1000,1), " (", 
                  round(l.df[j,cols[i]]/1000,1), "-", 
                  round(u.df[j,cols[i]]/1000,1), ")", sep=""))  
    } else {
      print(paste(trt, ": ", 
                  round(mean.df[j,cols[i]]), " (", 
                  round(l.df[j,cols[i]]), "-", 
                  round(u.df[j,cols[i]]), ")", sep=""))
      print(paste(trt, ": ", 
                  round(mean.df[j,cols[i]],1), " (", 
                  round(l.df[j,cols[i]],1), "-", 
                  round(u.df[j,cols[i]],1), ")", sep=""))
    }
  }
}
