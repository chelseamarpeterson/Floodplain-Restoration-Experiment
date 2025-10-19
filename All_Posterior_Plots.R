path_to_repo= "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch2_Floodplain_Experiment/Floodplain-Experiment-Repo"
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
# define treatment and variable names, then load datasets 

# treatment names
trt.letters = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")
n.t = length(trt.letters)

# read posterior distributions for soil properties, carbon stocks, and species richness
soil.df = read.csv("Soil_Analysis/Posteriors/All_Soil_Posterior_Intervals_ChainCount5.csv")
stock.df = read.csv("Tree_Analysis/Posteriors/Carbon_Stocks_Richness_Means_Intervals_5Chains.csv")

# soil variable labels
soil.vars = unique(soil.df$variable)
soil.var.labs = unique(soil.df$variable.label)

# stock and richness variables
stock.vars = unique(stock.df$variable)
stock.var.labs = unique(stock.df$variable.label)

################################################################################
# print results for appendix tables

# C stock and richness means and 90% intervals
for (i in 1:length(stock.vars)) {
  var.rows = which(stock.df$variable == stock.vars[i])
  print(cbind(stock.df[var.rows,c("variable","treatment")],
              signif(stock.df[var.rows,c("posterior.mean","X5","X95")],3)))
}

# soil means and 90% intervals
for (i in 1:length(soil.vars)) {
  print(cbind(soil.df[which(soil.df$variable == soil.vars[i]),c("variable","treatment")],
              signif(soil.df[which(soil.df$variable == soil.vars[i]),c("posterior.mean","X5","X95")],3)))
}

################################################################################
# make combined plot for carbon concentrations and stocks

# IPCC estimates
ipcc.df = read.csv("Carbon_Calculations/IPCC_Carbon_Estimates.csv")

# stacked plot for carbon fractions
c.conc.labs = c("tic.percent","poc.percent","moc.percent")
c.conc.labs.new = c("TIC","POC","MOC")
soil.df$lab = rep(0, nrow(soil.df))
for (i in 1:3) { soil.df$lab[which(soil.df$variable == c.conc.labs[i])] = c.conc.labs.new[i] }
s.palette <- c(brewer.pal(9,"Greys")[5],brewer.pal(11,"BrBG")[c(2,1)])
p.c.concentrations = ggplot(soil.df[which(soil.df$lab %in% c.conc.labs.new),], 
                            aes(y=factor(treatment, levels=trt.names), 
                                x=posterior.mean, 
                                fill=factor(lab, levels=c.conc.labs.new))) + 
                            geom_bar(stat="identity",
                                     position="stack") +
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
p.c.concentrations
ipcc.soil.vars = c("Annual crops","Revegetated cropland","Natural wetland")
p.c.stocks = ggplot(stock.df[which(stock.df$variable %in% c.conc.labs.new),], 
                    aes(y=factor(treatment, levels=trt.names), 
                        x=posterior.mean, 
                        fill=factor(variable, levels=c.conc.labs.new))) +
                    geom_bar(stat="identity",position="stack") + 
                    scale_fill_manual(values=s.palette) +
                    labs(fill="",y="",x="Stock (Mg C/ha)",title="") +
                    geom_vline(data=ipcc.df, 
                               aes(xintercept=soc.value, 
                                   color=factor(soc.type,levels=ipcc.soil.vars),
                                   linetype=factor(soc.type,levels=ipcc.soil.vars)), 
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
# plot non-carbon chemical and physical variables

# plot all variables
#ggplot(soil.df, aes(y=factor(t, levels=trt.names), 
#                   x=mean)) + 
#        geom_errorbarh(aes(xmin=X5, xmax=X95), color="black") +
#        geom_point() + 
#        facet_wrap(.~v, scales="free_x",ncol=5)

# plot select variables 
phys.chem.vars = c("Temperature (C)","Gravitational moisture (%)","Bulk density (g/cm3)",
                   "TN (%)","C:N Ratio","Sand (%)","Silt (%)","Clay (%)",
                   "pH","P (ppm)","NO3-N (ppm)","NH4-N (ppm)",
                   "K (meq/100 g)","Ca (meq/100 g)","Mg (meq/100 g)",
                   ">=4.75 mm","2-4.75 mm","0.250-2 mm","0.053-0.250 mm","<0.053 mm")
phys.chem.labs = c("Temperature (C)","Gravitational moisture (%)","Bulk density (g/cm3)",
                   "TN (%)","C:N ratio","Sand (%)","Silt (%)","Clay (%)",
                   "pH","P (ppm)","NO3-N (ppm)","NH4-N (ppm)",
                   "K (meq/100 g)","Ca (meq/100 g)","Mg (meq/100 g)",
                   "\u2265 4.75 mm (%)","2-4.75 mm (%)","0.250-2 mm (%)",
                   "0.053-0.250 mm (%)","< 0.053 mm (%)")
df.plot = soil.df[which(soil.df$variable.label %in% phys.chem.vars),]
df.plot$label.new = rep(0, nrow(df.plot))
for (i in 1:length(phys.chem.vars)) { df.plot$label.new[which(df.plot$variable.label == phys.chem.vars[i])] = phys.chem.labs[i] }
p.chem.phys = ggplot(df.plot, 
                     aes(y=factor(treatment, levels=trt.names), 
                         x=posterior.mean)) + 
                     geom_errorbar(aes(xmin=`X5`, xmax=`X95`), 
                                   width=0.25, color="black", 
                                   position=position_dodge(width=0.5)) +
                     geom_point(size=1, position=position_dodge(width=0.5)) +
                     facet_wrap(.~factor(label.new, 
                                         levels=phys.chem.labs), 
                                ncol=5, scales="free_x") +
                     theme(panel.spacing.x = unit(0.4, "cm"),
                           text = element_text(size=12)) +
                     labs(y="",x="Posterior estimate") + 
                     scale_x_continuous(breaks = breaks_extended(n = 4)) 
p.chem.phys
ggsave("Figures/Figure4_All_Except_CEC_MWD.jpeg", 
       plot=p.chem.phys,width=27,height=18,units="cm",dpi=600)

################################################################################
# make combined plot for C stocks in biomass, debris, and soil + vegetation,
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
ipcc.vars.l = c("Restored temperate forest","Natural temperate forest")
p.l = ggplot(stock.df[which(stock.df$variable.label %in% stack.vars.l),], 
             aes(y=factor(treatment, levels=trt.names), 
                 x=posterior.mean, 
                 fill=factor(variable.label, levels=stack.vars.l))) +
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
abg.df = stock.df[which(stock.df$variable.label == "Aboveground woody biomass"),]
bg.df = stock.df[which(stock.df$variable.label == "Belowground woody biomass"),]
round(bg.df$posterior.mean/(bg.df$posterior.mean+abg.df$posterior.mean)*100,1)

# debris plot
stack.vars.d.old = c("Herbaceous litter","Fine woody debris (< 2.5 cm)",
                     "Fine woody debris (2.5-7.5 cm)","Coarse woody debris (>= 7.6 cm)",
                     "Dead stems (< 2.5 cm)","Standing dead trees (>= 2.5 cm)")
stack.vars.d.new = c("Herbaceous litter","Fine woody debris (< 2.5 cm)",
                     "Fine woody debris (2.5-7.5 cm)","Coarse woody debris (\u2265 7.6 cm)",
                     "Dead stems (< 2.5 cm)","Standing dead trees (\u2265 2.5 cm)")
ipcc.vars.d = c("Restored temperate forest","Natural temperate forest")
stock.df.debris = stock.df[which(stock.df$variable.label %in% stack.vars.d.old),]
stock.df.debris$label.new = rep(0, nrow(stock.df.debris))
for (i in 1:6) { stock.df.debris$label.new[which(stock.df.debris$variable.label == stack.vars.d.old[i])] = stack.vars.d.new [i] }
p.d = ggplot(stock.df.debris, 
             aes(y=factor(treatment, levels=trt.names),
                 x=posterior.mean, 
                 fill=factor(label.new, levels=stack.vars.d.new))) +
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
stack.vars.e = c("Soil inorganic carbon","Soil organic carbon","Total dead vegetation","Total live vegetation")
ipcc.vars.e = c("Restored forested wetland","Natural forested wetland")
p.e = ggplot(stock.df[which(stock.df$variable.label %in% stack.vars.e),], 
             aes(y=factor(treatment, levels=trt.names),
                 x=posterior.mean,  
                 fill=factor(variable.label, 
                             levels=stack.vars.e))) +
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
                              labels=c("TIC",
                                       "SOC",
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
stack.vars.r = c("Total richness","Herbaceous layer richness","Tree layer richness")
stack.vars.new = c("Herbaceous layer only","Tree layer only","Tree and herbaceous layer")
df.sp = stock.df[which(stock.df$variable.label %in% stack.vars.r),]
df.sp.wide = pivot_wider(df.sp, 
                         id_cols=treatment, 
                         names_from=variable.label, 
                         values_from=posterior.mean)
df.sp.wide$`Tree and herbaceous layer` = df.sp.wide$`Herbaceous layer richness` + df.sp.wide$`Tree layer richness` - df.sp.wide$`Total richness`
df.sp.new = df.sp.wide[,-which(colnames(df.sp.wide) == "Total richness")]
colnames(df.sp.new) = c("treatment",stack.vars.new)
df.sp.melt = melt(df.sp.new, id.vars=c("treatment"))
p.r = ggplot(data=df.sp.melt) + 
             geom_bar(aes(y=factor(treatment, levels=trt.names), 
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

################################################################################
# plot C stocks by species for snags, coarse woody debris, live trees,
# and herbaceous layer species

# read in posterior intervals
post.df.spp = read.csv("Tree_Analysis/Posteriors/Vegetation_Carbon_Stocks_Species_Means_Intervals_5Chains.csv")

# plot stacked means for live woody, snag/CWD area, and biomass stocks by species
type.rows = c("live.trees","dead.trees","coarse.woody.debris")
type.labs = c("Live trees (\u2265 2.5 cm)",
              "Standing dead trees (\u2265 2.5 cm)",
              "Coarse woody debris (\u2265 7.6 cm)")
post.df.woodC = post.df.spp[which(post.df.spp$type %in% type.rows),]
post.df.woodC$type.lab = rep(0, nrow(post.df.woodC))
for (i in 1:3) { post.df.woodC$type.lab[which(post.df.woodC$type == type.rows[i])] = type.labs[i] }
post.df.woodC$species = trimws(post.df.woodC$species)
all.sp = sort(unique(post.df.woodC$species))
p.woodyC = ggplot(post.df.woodC, 
                  aes(x=posterior.mean, 
                      y=factor(treatment, levels=trt.names),
                      fill=factor(species, levels=all.sp))) + 
                  geom_bar(stat="identity") +
                  facet_wrap(.~factor(type.lab, levels=type.labs), 
                             scales="free_x", ncol=1) +
                  labs(x="Posterior mean C stock (Mg/ha)", y="", fill="Species") +
                  guides(fill=guide_legend(ncol=1)) + 
                  theme(text = element_text(size=14))
p.woodyC
ggsave("Figures/FigureB1_Woody_C_Stocks_By_Species.jpeg", 
       plot=p.woodyC, width=17, height=19, units="cm")

# plot stacked means for herbaceous species grounds
post.df.herbC = post.df.spp[which(post.df.spp$type == "herbaceous.biomass"),]
spp.herbC = levels(factor(post.df.herbC$species))
p.herbC = ggplot(post.df.herbC, 
                 aes(x=posterior.mean,
                     y=factor(treatment, levels=trt.names), 
                     fill=factor(species, levels=spp.herbC))) + 
                 geom_bar(stat = "identity") + 
                 labs(y="", x="Posterior mean C stock (Mg/ha)", 
                      title = "",
                      fill="Species group") + 
                 theme(text = element_text(size=14))
p.herbC
ggsave("Figures/FigureB2_Herbaceous_C_Stocks_By_Species.jpeg", 
       plot=p.herbC, width=18, height=10, units="cm")

