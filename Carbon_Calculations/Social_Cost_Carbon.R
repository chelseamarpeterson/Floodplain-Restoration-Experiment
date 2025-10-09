path_to_repo= "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo"
setwd(path_to_repo)

library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)

### script that uses social cost of carbon to estimate absolute and relative 
### carbon benefits of each treatment

# treatment names
trts = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

## read in datasets

# social cost of carbon estimates 
scc.df = read.csv("Carbon_Calculations/SCC_Estimates.csv")
colnames(scc.df) = c("stat","per.CO2","per.CO2.C")
scc.df$stat = c("mean","lower","upper")
rownames(scc.df) = c("mean","lower","upper")

# establishment costs
est.df = read.csv("Carbon_Calculations/Treatment_Establishment_Costs.csv")
colnames(est.df) = c("trt","cost.2019","cost.2025")

# ecosystem carbon estimates
stock.df = read.csv("Tree_Analysis/Posteriors/CarbonRichness_Means_Intervals_5Chains.csv")
sort(unique(stock.df$var))
ecoC.df = stock.df[which(stock.df$var == "Total ecosystem"),c("trt","mean","X5","X95")]
colnames(ecoC.df) = c("trt","mean","lower","upper")

## combine datasets to estimate carbon benefit

# melt ecosystem carbon estimates by statistic
ecoC.df.melt = melt(ecoC.df, id.vars=c("trt"), variable.name="stat",value.name="stock")

# estimate carbon benefit of each treatment
ecoC.df.melt$carbon.benefit = rep(0, nrow(ecoC.df.melt))
stats = c("mean","lower","upper")
for (i in 1:3) {
  stat.i = stats[i]
  stat.id = which(ecoC.df.melt$stat == stat.i)
  ecoC.df.melt$carbon.benefit[stat.id] = ecoC.df.melt$stock[stat.id]*scc.df[stat.i,"per.CO2.C"]
}

# estimate net benefit of each treatment by subtracting the establishment cost
ecoC.df.melt$net.benefit = rep(0, nrow(ecoC.df.melt))
ecoC.df.melt$breakeven.scc = rep(0, nrow(ecoC.df.melt))
for (i in 1:6) {
  trt.i = trts[i]
  cost.i = est.df[which(est.df$trt == trt.i),"cost.2025"]
  trt.id = which(ecoC.df.melt$trt == trt.i)
  ecoC.df.melt[trt.id,"net.benefit"] = ecoC.df.melt[trt.id,"carbon.benefit"] - cost.i
  ecoC.df.melt[trt.id,"breakeven.scc"] = cost.i/ecoC.df.melt[trt.id,"stock"]* 12.01/44.01
}

# get species richness estimates
n.df = stock.df[which(stock.df$var == "Total richness"),c("trt","mean","X5","X95")]
colnames(n.df) = c("trt","mean","lower","upper")
n.df.melt = melt(n.df, id.vars=c("trt"), variable.name="stat",value.name="richness")

# join richness estimates 
ecoC.n.df.join = left_join(ecoC.df.melt, n.df.melt, by=c("trt","stat"))

# plot net carbon benefit v. richness
mean.df = ecoC.n.df.join[which(ecoC.n.df.join$stat == "mean"),]
scc.df$stat = c("Mean","Lower (5%)","Upper (95%)")
colnames(scc.df)[1] = "Social Cost of Carbon Estimate"
p1 = ggplot(data=mean.df, 
            aes(y=richness, x=stock, 
                color=factor(trt, levels=trts))) +
            geom_point(size=4) + 
            guides(color="none") + 
            scale_y_continuous(breaks=seq(7,21,by=2),limits=c(7,21)) +
            labs(x="Ecosystem carbon stock (Mg/ha)",
                 y="Total species richness") +
            theme(text=element_text(size=14))
p2 = ggplot(data=mean.df, 
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
p3 = ggplot(data=mean.df, 
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
     geom_point(data=mean.df, 
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
cols = colnames(ecoC.n.df.join[3:7])
for (i in 1:length(cols)) {
  mean.df = ecoC.n.df.join[which(ecoC.n.df.join$stat == "mean"),c("trt",cols[i])]
  l.df = ecoC.n.df.join[which(ecoC.n.df.join$stat == "lower"),c("trt",cols[i])]
  u.df = ecoC.n.df.join[which(ecoC.n.df.join$stat == "upper"),c("trt",cols[i])]
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
