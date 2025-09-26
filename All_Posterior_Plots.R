path_to_repo= "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Chapter2/Floodplain-Experiment-Repo"
setwd(path_to_repo)

library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(rethinking)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggtext)
library(scales)

################################################################################
# Step 1: define treatment and variable names, load datasets 

# treatment names
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

# read posterior distributions for soil properties, carbon stocks, and species richness
df.soil = read.csv("Soil_Analysis/Posteriors/All_Soil_Posterior_Intervals_ChainCount5.csv")
df.stock = read.csv("Tree_Analysis/Posteriors/CarbonRichness_Means_Intervals_5Chains.csv")

# soil variable labels
soil.var.labs = c("SOM (%)","TOC (%)","TIC (%)","MAOC (%)","POC (\u2265 53 \U00B5m)","TC (%)","TN (%)","C:N Ratio",
                  "Temperature (C)","Gravitational moisture (%)","Bulk density (g/cm3)",
                  "Sand (%)","Silt (%)","Clay (%)","NO3 (ppm)","NH4 (ppm)","P (ppm)",
                  "K (ppm)","Ca (ppm)","Mg (ppm)","K (meq)","Ca (meq)","Mg (meq)",
                  "CEC (meq/100 g)","pH","Floating POM",
                  ">4.75 mm","2 - 4.75 mm","250 \u03bcm - 2 mm","53 - 250 \u03bcm",
                  "<53 \u03bcm","Mean-weight diameter (mm)")

# stock and richness variables
stock.var.labels = c("Fine woody debris (2.5-7.5 cm)","Coarse woody debris (\u2265 7.6 cm)",
                     "Standing dead trees (\u2265 2.5 cm)","Live F. pennsylvanica (\u2265 2.5 cm)",
                     "Standing dead F. pennsylvanica (\u2265 2.5 cm)","Diff. F. pennsylvanica (Live-Dead)",
                     "Dead stems (< 2.5 cm)","Live stems (< 2.5 cm)","Live trees (\u2265 2.5 cm)",
                     "Herbaceous biomass","Mixed biomass","P. arundinacea biomass",
                     "H. japonicus biomass","Herbaceous litter","Fine woody debris",
                     "Fine woody debris C:N ratio","Herbaceous litter C:N ratio","Herbaceous biomass C:N ratio",
                     "MAOC","POC","SOC","TIC","TC",
                     "Herbaceous species","Tree species","Total richness",
                     "Total live vegetation","Total dead vegetation","Total vegetation","Total ecosystem",
                     "Live trees w/ F. pennsylvanica","Standing dead trees w/o F. pennsylvanica",
                     "Total live w/ F. pennsylvanica","Total dead w/o F. pennsylvanica",
                     "Total vegetation w/o F. pennsylvanica death","Total ecosystem w/o F. pennsylvanica death")

################################################################################
# Step 2: print results for appendix tables

# soil means and 90% intervals
i = 33
cbind(df.soil[which(df.soil$v == soil.var.labs[i]),c("v","t")],
      round(df.soil[which(df.soil$v == soil.var.labs[i]),c("mean","X5","X95")],1))
cbind(df.soil[which(df.soil$v == soil.var.labs[i]),c("v","t")],
      round(df.soil[which(df.soil$v == soil.var.labs[i]),c("mean","X5","X95")],2))
cbind(df.soil[which(df.soil$v == soil.var.labs[i]),c("v","t")],
      round(df.soil[which(df.soil$v == soil.var.labs[i]),c("mean","X5","X95")],3))

# C stock and richness means and 90% intervals
i = 9
print(cbind(df.stock[which(df.stock$var == stock.var.labels[i]),c("var","trt")],
            round(df.stock[which(df.stock$var == stock.var.labels[i]),c("mean","X5","X95")],1)))
print(cbind(df.stock[which(df.stock$var == stock.var.labels[i]),c("var","trt")],
            round(df.stock[which(df.stock$var == stock.var.labels[i]),c("mean","X5","X95")],2)))

which(stock.var.labels == "Live trees (\u2265 2.5 cm)")
df.plot = df.stock[df.stock$var %in% c("Live trees (\u2265 2.5 cm)","Total live vegetation"),]
ggplot(df.plot, aes(y=trt, x=mean, color=var)) + 
                geom_point(position=position_dodge(0.5)) +
                geom_errorbarh(aes(y=trt, , color=var,
                                   xmin=`X5`,xmax=`X95`),
                               position=position_dodge(0.5),
                               height=0.2)

################################################################################
# Step 3: make combined plot for carbon concentrations and stocks

# IPCC estimates
ipcc.df = read.csv("Carbon_Calculations/IPCC_Carbon_Estimates.csv")

# stacked plot for carbon fractions
c.labs = c("TIC (%)","POC (%)","MAOC (%)")
df.stack.c = df.soil[which(df.soil$v %in% c.labs),]
df.stack.c$position = rep(0, nrow(df.stack.c))
for (i in 1:6) {
  df.trt = df.stack.c[which(df.stack.c$t == trt.names[i]),]
  trt.id = which(df.stack.c$t == trt.names[i])
  df.trt$position = cumsum(df.trt[seq(3,1,-1),"mean"])[seq(3,1,-1)] - 1.5*df.trt$mean
  df.stack.c[trt.id,"position"] = df.trt$position
}
c.labs.new = c("TIC","POC","MAOC")
df.stack.c$lab = rep(0, nrow(df.stack.c))
for (i in 1:3) { df.stack.c$lab[which(df.stack.c$v == c.labs[i])] = c.labs.new[i] }
s.palette <- c(brewer.pal(9,"Greys")[5],brewer.pal(11,"BrBG")[c(2,1)])
p.c.concentrations = ggplot(df.stack.c[which(df.stack.c$v %in% c.labs),], 
                            aes(y=factor(t, levels=trt.names), x=mean, 
                                fill=factor(lab, levels=c.labs.new))) + 
                            geom_bar(stat="identity",position="stack") +
                            labs(x="Concentration (% [g C/g soil])",
                                 y="",fill="Soil carbon fraction") + 
                            scale_fill_manual(values=s.palette) +
                            theme(text = element_text(size=14),
                                  plot.margin=unit(c(1,1,1,1),"lines")) +
                            coord_cartesian(xlim = c(0,4.6), clip="off") +
                            geom_label(x=5.12,y=6.3,label="a",
                                       color="black",fill=alpha("white",0.9),
                                       label.r=unit(0,"pt"),label.size=0,
                                       size=8,fontface="bold") + 
                            guides(fill="none")
stack.vars4 = c("TIC","POC","MAOC")
ipcc.vars4 = c("Annual crops","Revegetated cropland","Natural wetland")
p.c.stocks = ggplot(df.stock[which(df.stock$var %in% stack.vars4),], 
                    aes(y=factor(trt, levels=trt.names), 
                        x=mean, fill=factor(var, levels=stack.vars4))) +
                    geom_bar(stat="identity",position="stack") + 
                    scale_fill_manual(values=s.palette) +
                    labs(fill="",y="",x="Stock (Mg C/ha)",title="") +
                    geom_vline(data=ipcc.df, 
                               aes(xintercept=soc.value, 
                                   color=factor(soc.type,levels=ipcc.vars4),
                                   linetype=factor(soc.type,levels=ipcc.vars4)), 
                               linewidth=1.5) +
                    scale_color_manual(values=c("red","yellow1","royalblue1")) + 
                    scale_linetype_manual(values=c("solid","dotdash","dashed")) +
                    guides(color=guide_legend(title="IPCC soil organic carbon"),
                           linetype=guide_legend(title="IPCC soil organic carbon")) +
                    theme(text=element_text(size=14), 
                          axis.text.y=element_blank(),
                          legend.key.size=unit(0.7,'cm'),
                          plot.margin=unit(c(1,1,1,1),"lines"),
                          legend.key=element_rect(fill="darkgrey")) +
                    coord_cartesian(xlim=c(0,130), clip="off") +
                    scale_x_continuous(breaks=c(0,30,60,90,120))+
                    geom_label(x=145.5,y=6.3,label="b",
                               color="black",fill=alpha("white",0.9),
                               label.r=unit(0,"pt"),label.size=0,
                               size=8,fontface="bold")
p.c.all = p.c.concentrations + p.c.stocks
p.c.all =  p.c.all + theme(plot.margin = margin(0, 0, 0, 0, "cm"))
p.c.all
ggsave("Figures/Figure3_Soil_Carbon_Concentrations_and_Stocks.jpeg", 
       plot=p.c.all, width=30, height=12, units="cm",dpi=600)

################################################################################
# Step 4: make combined plot for CEC, carbon fractions, texture, and aggregates
# (figure not used in paper)

# stacked plot for aggregate size distribution
t.size = 14
ag.labs = c(">4.75 mm","2 - 4.75 mm","250 \u03bcm - 2 mm","53 - 250 \u03bcm","<53 \u03bcm")
df.stack.ag = df.soil[which(df.soil$v %in% ag.labs),]
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
              theme(text=element_text(size=t.size),
                    axis.text.y = element_blank()) +
              coord_cartesian(xlim = c(0,100), clip="off") +
              theme(plot.margin=unit(c(1,1,1,1),"lines")) +
              geom_label(x=111.4,y=6.3,label="d",
                         color="black",fill=alpha("white",0.9),
                         label.r=unit(0,"pt"),label.size=0,
                         size=8,fontface="bold")

# stacked plot for particle size distribution
text.labs = c("Sand (%)","Silt (%)","Clay (%)")
df.stack.text = df.soil[which(df.soil$v %in% text.labs),]
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
                theme(text = element_text(size=t.size)) +
                coord_cartesian(xlim = c(0,100), clip="off") +
                theme(plot.margin=unit(c(1,1,1,1),"lines")) +
                geom_label(x=111.1,y=6.3,label="c",
                           color="black",fill=alpha("white",0.9),
                           label.r=unit(0,"pt"),label.size=0,
                           size=8,fontface="bold")

# stacked plot for cation exchange
meq.labs = c("K (meq)","Mg (meq)","Ca (meq)")
df.stack.meq = df.soil[which(df.soil$v %in% meq.labs),]
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
               theme(text = element_text(size=t.size),
                     axis.text.y = element_blank()) +
               coord_cartesian(xlim = c(0,41), clip="off") +
               theme(plot.margin=unit(c(1,1,1,1),"lines")) +
               geom_label(x=45.7,y=6.3,label="b",
                          color="black",fill=alpha("white",0.9),
                          label.r=unit(0,"pt"),label.size=0,
                          size=8,fontface="bold")

p.upper = p.c.concentrations + p.meq + plot_spacer() + plot_layout(widths=c(0.5,0.5,0))
p.lower = p.text + p.ag + plot_spacer() + plot_layout(widths=c(0.5,0.5,0))
p.all =  p.upper/p.lower + 
         plot_layout(guides = "collect") + 
         theme(legend.position = "right",
               plot.margin = unit(c(0,0,0,0),"mm"))
p.all

################################################################################
# Step 5: plot non-carbon chemical and physical variables

# plot all variables
ggplot(df.soil, aes(y=factor(t, levels=trt.names), 
                   x=mean)) + 
        geom_errorbarh(aes(xmin=X5, xmax=X95), color="black") +
        geom_point() + 
        facet_wrap(.~v, scales="free_x",ncol=5)

# plot select variables 
phys.chem.vars = c("Temperature (C)","Gravitational moisture (%)","Bulk density (g/cm3)",
                   "TN (%)","C:N Ratio","Sand (%)","Silt (%)","Clay (%)",
                   "pH","P (ppm)","NO3 (ppm)","NH4 (ppm)",
                   "K (meq)","Ca (meq)","Mg (meq)",">4.75 mm","2 - 4.75 mm",
                   "250 \u03bcm - 2 mm","53 - 250 \u03bcm","<53 \u03bcm")
phys.chem.labs = c("Temperature (C)","Gravitational moisture (%)","Bulk density (g/cm3)",
                   "TN (%)","C:N ratio","Sand (%)","Silt (%)","Clay (%)",
                   "pH","P (ppm)","NO3-N (ppm)","NH4-N (ppm)",
                   "K (meq/100 g)","Ca (meq/100 g)","Mg (meq/100 g)",
                   ">4.75 mm (%)","2 - 4.75 mm (%)","250 \u03bcm - 2 mm (%)",
                   "53 - 250 \u03bcm (%)","<53 \u03bcm (%)")
df.plot = df.soil[which(df.soil$v %in% phys.chem.vars),]
df.plot$label = rep(0, nrow(df.plot))
for (i in 1:length(phys.chem.vars)) { df.plot$label[which(df.plot$v == phys.chem.vars[i])] = phys.chem.labs[i] }

p.chem.phys = ggplot(df.plot, aes(y=factor(t, levels=trt.names), x=mean)) + 
                    geom_errorbar(aes(xmin=`X5`, xmax=`X95`), 
                                  width=0.25, color="black", 
                                  position=position_dodge(width=0.5)) +
                    geom_point(size=1, position=position_dodge(width=0.5)) +
                    facet_wrap(.~factor(label, levels=phys.chem.labs), 
                               ncol=5, scales="free_x") +
                    theme(panel.spacing.x = unit(0.4, "cm"),
                          text = element_text(size=12)) +
                    labs(y="",x="Posterior estimate") + 
                    scale_x_continuous(breaks = breaks_extended(n = 4)) 
p.chem.phys
ggsave("Figures/Figure4_All_Except_CEC_MWD.jpeg", 
       plot=p.chem.phys,width=27,height=18,units="cm",dpi=600)

################################################################################
# Step 6: make combined plot for C stocks in biomass, debris, and soil + vegetation,
# along with species richness

# define color palletes
v.palette <- brewer.pal(11,"RdYlGn")[c(7,9,11)]
d.palette <- brewer.pal(9,"YlOrBr")[seq(3,8)]
r.palette <- tail(brewer.pal(8,"BuPu"),3)

# update ipcc.df
ipcc.df = ipcc.df[1:2,]

#  live vegetation plot
stack.vars.l = c("Herbaceous biomass",
                 "Belowground woody biomass",
                 "Aboveground woody biomass")
ipcc.vars.l = c("Restored temperate forest","Mature temperate forest")
p.l = ggplot(df.stock[which(df.stock$var %in% stack.vars.l),], 
             aes(y=factor(trt, levels=trt.names), x=mean, 
                 fill=factor(var, levels=stack.vars.l))) +
              geom_bar(stat="identity",position="stack") + 
              labs(fill="",y="",
                   x="Posterior mean stock (Mg C/ha)",title="") + 
              geom_vline(data=ipcc.df,
                         aes(xintercept=woody.value, 
                             color=factor(woody.type,
                                          levels=ipcc.vars.l),
                             linetype=factor(woody.type,
                                             levels=ipcc.vars.l)), 
                         linewidth=1.5) +
              scale_fill_manual(values=c("olivedrab3","olivedrab4",v.palette[3])) +
              scale_color_manual(values=c("red","royalblue1")) +
              scale_linetype_manual(values=c("solid","dashed")) + 
              guides(fill=guide_legend(order=1),
                     color=guide_legend(title="IPCC total woody biomass",order=2),
                     linetype=guide_legend(title="IPCC total woody biomass",order=2)) +
              theme(text=element_text(size=14), 
                    legend.key.size=unit(0.7,'cm'),
                    legend.background=element_rect(fill="transparent", color=NA),
                    plot.margin=unit(c(0,0,1,0),"lines"),
                    plot.background=element_rect(color="black",linewidth=1),
                    legend.key=element_rect(fill="darkgrey")) +
              coord_cartesian(xlim=c(0,180),clip="off") +
              geom_label(x=338,y=6.6,label="a",
                         color="black",fill=alpha("white",0.9),
                         label.r=unit(0,"pt"),label.size=0,
                         size=10,fontface="bold")
p.l
abg.df = df.stock[which(df.stock$var == "Aboveground woody biomass"),]
bg.df = df.stock[which(df.stock$var == "Belowground woody biomass"),]
bg.df$mean/(bg.df$mean+abg.df$mean)*100

# debris plot
stack.vars.d = c("Herbaceous litter","Fine woody debris (< 2.5 cm)",
                "Fine woody debris (2.5-7.5 cm)","Coarse woody debris (\u2265 7.6 cm)",
                "Dead stems (< 2.5 cm)","Standing dead trees (\u2265 2.5 cm)")
ipcc.vars.d = c("Restored temperate forest","Mature temperate forest")
p.d = ggplot(df.stock[which(df.stock$var %in% stack.vars.d),], 
             aes(y=factor(trt, levels=trt.names), x=mean, 
                 fill=factor(var, levels=stack.vars.d))) +
             geom_bar(stat="identity",position="stack") + 
             labs(fill="",y="",
                  x="Posterior mean stock (Mg C/ha)",title="") + 
             geom_vline(data=ipcc.df, 
                        aes(xintercept=debris.value,
                            color=factor(debris.type,
                                         levels=ipcc.vars.d),
                            linetype=factor(debris.type,
                                            levels=ipcc.vars.d)), 
                        linewidth=1.5) +
             scale_color_manual(values=c("red","royalblue1")) +
             scale_linetype_manual(values=c("solid","dashed")) + 
             scale_fill_manual(values=d.palette) + 
             guides(fill=guide_legend(order=1),
                    color=guide_legend(title="IPCC litter\nand woody debris",order=2),
                    linetype=guide_legend(title="IPCC litter\nand woody debris",order=2)) +
             theme(text=element_text(size=14), 
                   legend.key.size=unit(0.7,'cm'),
                   legend.background=element_rect(fill="transparent", color=NA),
                   plot.margin=unit(c(0,0,1,0),"lines"),
                   plot.background=element_rect(color="black",linewidth=1),
                   legend.key=element_rect(fill="darkgrey")) +
             coord_cartesian(xlim = c(0,30.5), clip="off") +
             geom_label(x=58.9,y=6.6,label="b",
                        color="black",fill=alpha("white",0.9),
                        label.r=unit(0,"pt"),label.size=0,
                        size=10,fontface="bold")
p.d

# total ecosystem
stack.vars.e = c("TIC","SOC","Total dead vegetation","Total live vegetation")
ipcc.vars.e = c("Restored forested wetland","Mature forested wetland")
p.e = ggplot(df.stock[which(df.stock$var %in% stack.vars.e),], 
             aes(y=factor(trt, levels=trt.names), x=mean, 
                 fill=factor(var, levels=stack.vars.e))) +
            geom_bar(stat="identity",position="stack") + 
            labs(fill="",y="",
                 x="Posterior mean stock (Mg C/ha)",title="") +
            geom_vline(data=ipcc.df, 
                       aes(xintercept=total.value, 
                           color=factor(total.type,
                                        levels=ipcc.vars.e),
                           linetype=factor(total.type,
                                           levels=ipcc.vars.e)), 
                       linewidth=1.5) +
            scale_fill_manual(values=c(brewer.pal(9,"Greys")[5],
                                       s.palette[3],d.palette[4],
                                       v.palette[3]),
                              labels=c("TIC","SOC",
                                       "Litter and woody debris",
                                       "Living biomass")) +
            scale_color_manual(values=c("red","royalblue1")) +
            scale_linetype_manual(values=c("solid","dashed")) +
            guides(fill=guide_legend(order=1),
                   color=guide_legend(title="IPCC total organic C",order=2),
                   linetype=guide_legend(title="IPCC total organic C",order=2)) +
            theme(text=element_text(size=14), 
                  legend.key.size=unit(0.7,'cm'),
                  legend.background=element_rect(fill="transparent", color=NA),
                  plot.margin=unit(c(0,0,1,0),"lines"),
                  plot.background=element_rect(color="black",linewidth=1),
                  legend.key=element_rect(fill="darkgrey")) +
            coord_cartesian(xlim = c(0,280), clip="off") +
            geom_label(x=525,y=6.6,label="c",
                       color="black",fill=alpha("white",0.9),
                       label.r=unit(0,"pt"),label.size=0,
                       size=10,fontface="bold")
p.e

# richness plot
stack.vars.r = c("Total richness","Herbaceous species","Tree species")
stack.vars.new = c("Herbaceous layer only","Tree layer only","Tree and herbaceous layer")
df.sp = df.stock[which(df.stock$var %in% stack.vars.r),]
df.sp.wide = pivot_wider(df.sp, id_cols=trt, names_from=var, values_from=mean)
df.sp.wide$`Tree and herbaceous layer` = df.sp.wide$`Tree species` + df.sp.wide$`Herbaceous species` - df.sp.wide$`Total richness`
df.sp.new = df.sp.wide[,-which(colnames(df.sp.wide) == "Total richness")]
colnames(df.sp.new) = c("trt",stack.vars.new)
df.sp.melt = melt(df.sp.new, id.vars=c("trt"))
p.r = ggplot(data=df.sp.melt) + 
             geom_bar(aes(y=factor(trt, levels=trt.names), 
                          x=value, 
                          fill=factor(variable, 
                                      levels=rev(stack.vars.new))),
                          stat="identity",
                          position="stack") +
             labs(fill="",y="",x="Posterior mean richness",title="") + 
             scale_fill_manual(values=r.palette) +
             theme(text = element_text(size=14),
                   legend.key.size = unit(0.7,'cm'),
                   legend.background = element_rect(fill="transparent", color=NA),
                   plot.margin=unit(c(0,0,1,0),"lines"),
                   plot.background=element_rect(color="black",linewidth=1),
                   legend.key=element_rect(fill="darkgrey")) +
             coord_cartesian(xlim=c(0,21),clip="off") +
             geom_label(x=40.5,y=6.6,label="d",
                        color="black",fill=alpha("white",0.9),
                        label.r=unit(0,"pt"),label.size=0,
                        size=10,fontface="bold")
p.r

# combine all plots
p.c.stocks = (p.l+theme(plot.margin=unit(c(0,0,0,0),"pt"))+p.d)/(p.e+theme(plot.margin=unit(c(0,16,0,0),"pt"))+p.r)
p.c.stocks
ggsave("Figures/Figure6_Veg_Ecosystem_Cstocks_Richness.jpeg", 
       plot=p.c.stocks, width=40, height=23, units="cm", dpi=1000)


