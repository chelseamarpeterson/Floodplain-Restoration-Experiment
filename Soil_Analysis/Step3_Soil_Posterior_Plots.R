path_to_soil_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Public-Repo/Soil_Analysis"
setwd(path_to_soil_folder)

library(ggplot2)
library(patchwork)
library(RColorBrewer)

################################################################################
### plot Bayesian model results for soil properties

# treatment names
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

# variable labels
var.labels = c("SOM (%)","TOC (%)","TIC (%)","POM-C (%)","MAOM-C (%)","TC (%)","TN (%)","C:N Ratio",
               "Temperature (C)","Gravitational moisture (%)","Bulk density (g/cm3)",
               "Sand (%)","Silt (%)","Clay (%)","NO3 (ppm)","NH4 (ppm)","P (ppm)",
               "K (ppm)","Ca (ppm)","Mg (ppm)","K (meq)","Ca (meq)","Mg (meq)",
               "CEC (meq/100 g)","pH","Floating POM",
               ">4.75 mm","2 - 4.75 mm","250 \u03bcm - 2 mm","53 - 250 \u03bcm",
               "<53 \u03bcm","Mean-weight diameter (mm)")

# read in HPDIs
df.int = read.csv("Posteriors/All_Soil_Posterior_Intervals_ChainCount5.csv")

# print results for tables
cbind(df.int[which(df.int$v == var.labels[6]),c("v","t")],
      round(df.int[which(df.int$v == var.labels[6]),c("mean","X5","X95")],2))

ggplot(df.int, aes(y=factor(t, levels=trt.names), 
                   x=mean)) + 
       geom_errorbarh(aes(xmin=X5, xmax=X95), color="black") +
       geom_errorbarh(aes(xmin=X25, xmax=X75), color="blue") +
       geom_point() + 
       facet_wrap(.~v, scales="free_x",ncol=5)

################################################################################
### make combined plot for CEC, carbon fractions, texture, and aggregates

# stacked plot for aggregate size distribution
t.size = 14
ag.labs = c(">4.75 mm","2 - 4.75 mm","250 \u03bcm - 2 mm","53 - 250 \u03bcm","<53 \u03bcm")
df.stack.ag = df.int[which(df.int$v %in% ag.labs),]
df.stack.ag$position = rep(0, nrow(df.stack.ag))
for (i in 1:6) {
  df.trt = df.stack.ag[which(df.stack.ag$t == trt.names[i]),]
  trt.id = which(df.stack.ag$t == trt.names[i])
  df.trt$position = cumsum(df.trt[seq(5,1,-1),"mean"])[seq(5,1,-1)] - 1.5*df.trt$mean
  df.stack.ag[trt.id,"position"] = df.trt$position
}
ag.labs.new = c(">4.75 mm","2-4.75 mm","0.250-2 mm","0.053-0.250 mm","<0.053 mm")
df.stack.ag$lab = rep(0, nrow(df.stack.ag))
for (i in 1:5) {  df.stack.ag$lab[which(df.stack.ag$v == ag.labs[i])]  = ag.labs.new[i]  }
ag.palette <- tail(brewer.pal(12,"Blues"),5)
p.ag = ggplot(df.stack.ag, aes(y=factor(t, levels=trt.names), 
                               x=mean, 
                               fill=factor(lab, levels=ag.labs.new))) + 
              geom_bar(stat="identity",position="stack") +
              geom_label(aes(label = round(mean,1)), 
                         nudge_x = df.stack.ag$position,
                         color="white",
                         size=2.5, show.legend = FALSE) +
              labs(x="Percent of total mass (%)",
                   y="",fill="Aggregate size class") +
              scale_fill_manual(values=ag.palette) +
              theme(text=element_text(size=t.size)) 
p.ag

# stacked plot for particle size distribution
text.labs = c("Sand (%)","Silt (%)","Clay (%)")
df.stack.text = df.int[which(df.int$v %in% text.labs),]
df.stack.text$position = rep(0, nrow(df.stack.text))
df.stack.text$lower = rep(0, nrow(df.stack.text))
df.stack.text$upper = rep(0, nrow(df.stack.text))
for (i in 1:6) {
  df.trt = df.stack.text[which(df.stack.text$t == trt.names[i]),]
  trt.id = which(df.stack.text$t == trt.names[i])
  df.trt$position = cumsum(df.trt[seq(3,1,-1),"mean"])[seq(3,1,-1)] - 1.5*df.trt$mean
  df.stack.text[trt.id,"position"] = df.trt$position
  df.stack.text[trt.id,"lower"] = (cumsum(df.trt[seq(3,1,-1),"mean"]) - (df.trt[seq(3,1,-1),"mean"]-df.trt[seq(3,1,-1),"5"]))[seq(3,1,-1)]
  df.stack.text[trt.id,"upper"] = (cumsum(df.trt[seq(3,1,-1),"mean"]) + (df.trt[seq(3,1,-1),"95"]-df.trt[seq(3,1,-1),"mean"]))[seq(3,1,-1)] 
}
text.labs.new = c("Sand (0.05-2.0 mm)","Silt (0.002-0.05 mm)","Clay (<0.002 mm)")
df.stack.text$lab = rep(0, nrow(df.stack.text))
for (i in 1:3) { df.stack.text$lab[which(df.stack.text$v == text.labs[i])] = text.labs.new[i] }
text.palette <- tail(brewer.pal(7,"BuPu"),3)
p.text = ggplot(df.stack.text,
                aes(y=factor(t, levels=trt.names), 
                    x=mean, 
                    fill=factor(lab, levels=text.labs.new))) + 
                geom_bar(stat = "identity", position = "stack") +
                geom_label(aes(label = round(mean,1)), 
                           nudge_x = df.stack.text$position,
                           color="white",#fill="white",
                           size=2.5, show.legend = FALSE) +
                labs(x="Percent of total mass (%)",y="", 
                     fill="Particle size class") +
                scale_fill_manual(values=text.palette) +
                theme(text = element_text(size=t.size))
p.text

# stacked plot for carbon fractions
c.labs = c("TIC (%)","POM-C (%)","MAOM-C (%)")
df.stack.c = df.int[which(df.int$v %in% c.labs),]
df.stack.c$position = rep(0, nrow(df.stack.c))
for (i in 1:6) {
  df.trt = df.stack.c[which(df.stack.c$t == trt.names[i]),]
  trt.id = which(df.stack.c$t == trt.names[i])
  df.trt$position = cumsum(df.trt[seq(3,1,-1),"mean"])[seq(3,1,-1)] - 1.5*df.trt$mean
  df.stack.c[trt.id,"position"] = df.trt$position
}
c.labs.new = c("TIC","POM-C","MAOM-C")
df.stack.c$lab = rep(0, nrow(df.stack.c))
for (i in 1:3) { df.stack.c$lab[which(df.stack.c$v == c.labs[i])] = c.labs.new[i] }
c.palette <- tail(brewer.pal(6,"Oranges"),3)
p.c = ggplot(df.stack.c[which(df.stack.c$v %in% c.labs),], 
             aes(y=factor(t, levels=trt.names), x=mean, 
                 fill=factor(lab, levels=c.labs.new))) + 
             geom_bar(stat="identity",position="stack") +
             geom_label(aes(label = round(mean,1)), 
                        nudge_x = df.stack.c$position,
                        color="white",#fill="white",
                        size=2.5, show.legend = FALSE) +
             labs(x="Concentration (% [g C/g soil])",
                  y="",fill="Soil carbon fraction") + 
             scale_fill_manual(values=c.palette) +
             theme(text = element_text(size=t.size))
p.c

# stacked plot for cation exchange
meq.labs = c("K (meq)","Mg (meq)","Ca (meq)")
df.stack.meq = df.int[which(df.int$v %in% meq.labs),]
df.stack.meq$position = rep(0, nrow(df.stack.meq))
for (i in 1:6) {
  df.trt = df.stack.meq[which(df.stack.meq$t == trt.names[i]),]
  trt.id = which(df.stack.meq$t == trt.names[i])
  df.trt$position = cumsum(df.trt[c(2,3,1),"mean"])[c(3,1,2)] - 1.5*df.trt$mean
  df.stack.meq[trt.id,"position"] = df.trt$position
}
meq.labs.new = c("K","Mg","Ca")
df.stack.meq$lab = rep(0, nrow(df.stack.meq))
for (i in 1:3) { df.stack.meq$lab[which(df.stack.meq$v == meq.labs[i])] = meq.labs.new[i] }
meq.palette <- tail(brewer.pal(7,"YlGn"),3)
p.meq = ggplot(df.stack.meq, aes(x=mean,
               y=factor(t, levels=trt.names), 
               fill=factor(lab, levels=meq.labs.new))) + 
               geom_bar(stat="identity",position="stack") +
               geom_label(aes(label = round(mean,1)), 
                          nudge_x = df.stack.meq$position,
                          color="white",
                         size=2.5, show.legend=FALSE) +
               labs(x="Contribution to CEC (meq/100 g)",y="",
                    fill="Cation") +
               scale_fill_manual(values=meq.palette) +
               theme(text = element_text(size=t.size))
p.meq

# combine carbon, meq, texture, & aggregate plots
p.all = (p.c + theme(plot.margin = unit(c(0,60,0,0), "pt")) + p.meq)/(p.text+ theme(plot.margin = unit(c(0,2,0,0), "pt")) +p.ag)
p.all
ggsave("Figures/CEC_And_Carbon_Texture_Aggregates.jpeg", 
       plot = p.all, width = 36, height = 18, units="cm")

################################################################################
# plot all other chemical and physical variables

# points with whiskers
vars = c("Temperature (C)","Gravitational moisture (%)","Bulk density (g/cm3)",
         "pH","C:N Ratio","TN (%)","NO3 (ppm)","NH4 (ppm)","P (ppm)")
var.labels = c("Temperature (C)","Gravitational moisture (%)","Bulk density (g/cm3)",
               "pH","C:N Ratio","Total N (%)","NO3-N (ppm)","NH4-N (ppm)","P (ppm)")
df.plot = df.int[which(df.int$v %in% vars),]
df.plot$lab = rep(0, nrow(df.plot))
for (i in 1:length(vars)) { df.plot$lab[which(df.plot$v == vars[i])] = var.labels[i] }
p.chem.phys = ggplot(df.plot, 
                     aes(y=factor(t, levels=trt.names), x=mean)) + 
                         geom_errorbar(aes(xmin=`X5`, xmax=`X95`), 
                                       width=0.25, color="black", 
                                       position=position_dodge(width=0.5)) +
                         geom_errorbar(aes(xmin=`X25`, xmax=`X75`), 
                                       width=0.25, color="blue", 
                                       position=position_dodge(width=0.5)) +
                         geom_point(position=position_dodge(width=0.5)) +
                         facet_wrap(.~factor(lab, levels=var.labels), 
                                    scales="free_x", ncol=3) +
                         theme(panel.spacing.x = unit(0.4, "cm")) +
                         labs(y="",x="Posterior estimate") + 
                         theme(text = element_text(size=12))
ggsave("Figures/Nutrients_BD_Temp_Moisture_pH.jpeg", 
       plot = p.chem.phys, width = 18, height = 16, units="cm")

## plot continuous distributions
new.df = all.a[["som"]]
colnames(new.df) = trt.names
all.melt = melt(new.df, id.vars=trt.names)
all.melt$var = "SOM (%)"
all.melt = all.melt[,2:4]
colnames(all.melt) = c("trt","value","var")
for (i in 2:n.v.plot) {
  v = var.order[i]
  new.df = all.a[[v]]
  colnames(new.df) = trt.names
  melt.df = melt(new.df, id.vars=trt.names) 
  melt.df$var = var.labels[i]
  melt.df = melt.df[,2:4]
  colnames(melt.df) = c("trt","value","var")
  all.melt = rbind(all.melt, melt.df)
}

all.melt$var = factor(all.melt$var, levels=var.labels)
ggplot(all.melt, aes(x=value, fill=trt)) + 
       geom_density(alpha=0.5, linewidth=0.2) + 
       facet_wrap(.~var, scales = "free") +
       labs(y="Density",x="") + 
       theme(legend.title=element_blank())