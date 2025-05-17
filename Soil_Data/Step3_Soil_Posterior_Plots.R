path_to_soil_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Public-Repo/Soil_Data"
setwd(path_to_soil_folder)

library(ggplot2)
library(patchwork)
library(RColorBrewer)

### plot Bayesian model results for all soil variables

# treatment names
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

# read in HPDIs
df.int = read.csv("Posteriors/All_Soil_Posterior_Intervals_ChainCount5.csv")

# print results for tables
cbind(df.int[which(df.int$v == var.labels[6]),c("v","t")],
      round(df.int[which(df.int$v == var.labels[6]),c("mean","X5","X95")],2))

# stacked plot for aggregate sizes
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

# compare particulate organic, total organic, and total carbon
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
                  y="",fill="Fraction") + 
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
               labs(x="Contribution to base saturation (meq/100 g)",y="",
                    fill="Cation") +
               scale_fill_manual(values=meq.palette) +
               theme(text = element_text(size=t.size))
p.meq

# combine carbon, meq, texture, & aggregate plots
p.all = (p.c + theme(plot.margin = unit(c(0,60,0,0), "pt")) + p.meq)/(p.text+ theme(plot.margin = unit(c(0,2,0,0), "pt")) +p.ag)
p.all
ggsave("Figures/CEC_And_Carbon_Texture_Aggregates.jpeg", 
       plot = p.all, width = 36, height = 20, units="cm")

# plot all other chemical and physical variables
var.set1 = c("Temperature (C)","Gravitational moisture (%)","Bulk density (g/cm3)",
             "pH","TN (%)","C:N Ratio","NH4 (ppm)",
             "NO3 (ppm)","P (ppm)","CEC (meq/100 g)")
ggplot(df.int[which(df.int$v %in% var.set1),], 
       aes(y=factor(t, levels=trt.names), x=mean)) + 
        geom_errorbar(aes(xmin=`X5`, xmax=`X95`), 
                      width=0.25, color="black", 
                      position=position_dodge(width=0.5)) +
        geom_errorbar(aes(xmin=`X25`, xmax=`X75`), 
                      width=0.25, color="blue", 
                      position=position_dodge(width=0.5)) +
        geom_point(position=position_dodge(width=0.5)) +
        facet_wrap(.~factor(v, levels=var.set1), 
                   scales="free_x", ncol=3) +
        theme(panel.spacing.x = unit(0.4, "cm")) +
        labs(y="",x="") + theme(text = element_text(size=12))

leave.out = c("Mean-weight diameter (mm)")
chem.vars = var.labels[-which(var.labels %in% c(text.labs, ag.labs[1:6],leave.out))]
df.chem = df.int[which(df.int$v %in% chem.vars),]
ggplot(df.chem, aes(y=factor(t, levels=trt.names), x=mean)) + 
      geom_errorbar(aes(xmin=`X5`, xmax=`X95`), 
                    width=0.25, color="black", 
                    position=position_dodge(width=0.5)) +
      geom_errorbar(aes(xmin=`X25`, xmax=`X75`), 
                    width=0.25, color="blue", 
                    position=position_dodge(width=0.5)) +
      geom_point(position=position_dodge(width=0.5)) +
      facet_wrap(.~factor(v, levels=chem.vars), scales="free_x", ncol=4) +
      theme(panel.spacing.x = unit(0.4, "cm")) +
      labs(y="",x="")

## plot all posteriors
ggplot(df.int, 
       aes(x=mean,
           y=factor(t, levels=trt.names))) + 
       geom_errorbar(aes(xmin=`5`, xmax=`95`), 
                     width=0.25, color="black", 
                     position=position_dodge(width=0.5)) +
       geom_errorbar(aes(xmin=`25`, xmax=`75`),
                     width=0.25, color="blue", 
                     position=position_dodge(width=0.5)) +
       geom_point(position=position_dodge(width=0.5)) +
       facet_wrap(.~factor(v, levels=var.labels), 
                  scales="free_x", ncol=4) +
       theme(panel.spacing.x = unit(0.4, "cm")) +
       labs(y="",x="")

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