path_to_soil_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Public-Repo/Soil_Data"
setwd(path_to_soil_folder)

library(ggplot2)
library(patchwork)

### plot Bayesian model results for all soil variables

# treatment names
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

# read in HPDIs
df.int = read.csv("Posteriors/All_Soil_Posterior_Intervals_ChainCount5.csv")

# print results for tables
cbind(df.int[which(df.int$v == var.labels[32]),c("v","t")],round(df.int[which(df.int$v == var.labels[32]),c("mean","5","95")],2))

# stacked plot for aggregate sizes
ag.labs = c(">4.75 mm","2 - 4.75 mm","250 \u03bcm - 2 mm","53 - 250 \u03bcm","<53 \u03bcm")
df.stack.ag = df.int[which(df.int$v %in% ag.labs),]
df.stack.ag$position = rep(0, nrow(df.stack.ag))
for (i in 1:6) {
  df.trt = df.stack.ag[which(df.stack.ag$t == trt.names[i]),]
  trt.id = which(df.stack.ag$t == trt.names[i])
  df.trt$position = cumsum(df.trt[seq(5,1,-1),"mean"])[seq(5,1,-1)] - 1.5*df.trt$mean
  df.stack.ag[trt.id,"position"] = df.trt$position
}
p.ag = ggplot(df.stack.ag, aes(y=factor(t, levels=trt.names), 
                               x=mean, 
                               fill=factor(v, levels=ag.labs))) + 
              geom_bar(stat="identity",position="stack") +
              geom_text(aes(label = round(mean,1)), 
                        nudge_x = df.stack.ag$position,
                        color="white",
                        size=3) +
              labs(x="Posterior mean of relative mass (%)",
                   y="",fill="Aggregate size class") +
              theme(text = element_text(size=14))

# stacked plot for texture
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
for (i in 1:3) { df.stack.text[which(df.stack.text$v == text.labs[i]),"v"]  = text.labs.new[i] }

p.text = ggplot(df.stack.text,
                aes(y=factor(t, levels=trt.names), 
                    x=mean, 
                    fill=factor(v, levels=text.labs.new))) + 
                geom_bar(stat = "identity", position = "stack") +
                geom_text(aes(label = round(mean,1)), 
                          nudge_x = df.stack.text$position,
                          color = "white",
                          size=3) +
                labs(x="Posterior mean of relative mass (%)",y="", 
                     fill="Particle size class") +
                scale_fill_manual(values=c("burlywood3","seashell4","coral3"))  +
                theme(text = element_text(size=14))
p.structure = p.text + p.ag
p.structure
ggsave("Figures/Particle_And_Aggregate_Size.jpeg", 
       plot = p.structure, width = 38, height = 10, units="cm")

# compare particulate organic, total organic, and total carbon
c.labs = c("SOC (%)","POC (%)","TC (%)")
df.stack.c = df.int[which(df.int$v %in% c.labs),]
omit.cols = c("v","t","sd")
omit.inds = which(colnames(df.stack.c) %in% omit.cols)
stat.vars = c("mean","X5","X95","X15","X85","X25","X75")
for (i in 1:6) {
  # get row for treatment
  df.trt = df.stack.c[which(df.stack.c$t == trt.names[i]),]
  
  # make row for inorganic c
  df.tic = data.frame(matrix(nrow=1, ncol=length(df.trt)))
  colnames(df.tic) = colnames(df.trt); df.tic$v = "TIC (%)"; df.tic$t = trt.names[i]
  df.tc = df.trt[which(df.trt$v == "TC (%)"),]
  df.soc = df.trt[which(df.trt$v == "SOC (%)"),]
  df.tic[1,stat.vars] = pmax(df.tc[1,stat.vars]-df.soc[1,stat.vars],0)
  df.stack.c = rbind(df.stack.c, df.tic)
  
  # make row for mineral c
  df.moc = data.frame(matrix(nrow=1, ncol=length(df.trt)))
  colnames(df.moc) = colnames(df.trt); df.moc$v = "MAOC (%)"; df.moc$t = trt.names[i]
  df.poc = df.trt[which(df.trt$v == "POC (%)"),]
  df.moc[1,stat.vars] = pmax(df.soc[1,stat.vars]-df.poc[1,stat.vars],0)
  df.stack.c = rbind(df.stack.c, df.moc)
}

df.stack.c[df.stack.c$v == "SOM (%)","v"] = "Organic (%)"
c.plot = c("Inorganic (%)","Organic (%)")
p.c = ggplot(df.stack.c[which(df.stack.c$v %in% c.plot),], 
             aes(y=factor(t, levels=trt.names), x=mean, 
                 fill=factor(v, levels=c.plot))) + 
  geom_bar(stat="identity",position="stack") +
  labs(x="Posterior mean (%)",y="",
       fill="", title="Carbon concentrations") + 
  scale_fill_manual(values=c("darkgray","darkorange4")) +
  theme(text = element_text(size=14))
p.c

c.plot2 = c("Organic (%)","Inorganic (%)","TC (%)")
ggplot(df.stack.c[which(df.stack.c$v == c.plot2[2]),], 
       aes(x=mean, y=factor(t, levels=trt.names),
           color=factor(v, levels=c.plot2[2]))) +
  geom_point(position=position_dodge(0.5)) +
  geom_errorbarh(aes(xmin=`5`,xmax=`95`,
                     y=factor(t, levels=trt.names),
                     color=factor(v, levels=c.plot2[2])),
                 position=position_dodge(0.5))

# stacked plot for 
meq.labs = c("K (meq)","Mg (meq)","Ca (meq)")
df.stack.meq = df.int[which(df.int$v %in% meq.labs),]
df.stack.meq$position = rep(0, nrow(df.stack.meq))
for (i in 1:6) {
  df.trt = df.stack.meq[which(df.stack.meq$t == trt.names[i]),]
  trt.id = which(df.stack.meq$t == trt.names[i])
  df.trt$position = cumsum(df.trt[c(2,3,1),"mean"])[c(3,1,2)] - 1.5*df.trt$mean
  df.stack.meq[trt.id,"position"] = df.trt$position
}
p.meq = ggplot(df.stack.meq, aes(y=factor(t, levels=trt.names), x=mean, fill=factor(v, levels=meq.labs))) + 
  geom_bar(stat="identity",position="stack") +
  geom_text(aes(label = round(mean,1)), 
            nudge_x = df.stack.meq$position,
            color="white") +
  labs(x="Posterior mean of relative mass (%)",y="",
       fill="Aggregate size class") +
  theme(axis.text.y = element_blank())
p.meq


var.set1 = c("Temperature (C)","Gravitational moisture (%)","Bulk density (g/cm3)",
             "pH","TN (%)","C:N Ratio","NH4 (ppm)",
             "NO3 (ppm)","P (ppm","CEC (meq/100 g)")
ggplot(df.int[which(df.int$v %in% var.set1),], 
       aes(y=factor(t, levels=trt.names), x=mean)) + 
  geom_errorbar(aes(xmin=`5`, xmax=`95`), width=0.25, color="black", position=position_dodge(width=0.5)) +
  geom_errorbar(aes(xmin=`25`, xmax=`75`), width=0.25, color="blue", position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width=0.5)) +
  facet_wrap(.~factor(v, levels=var.set1), scales="free_x", ncol=3) +
  theme(panel.spacing.x = unit(0.4, "cm")) +
  labs(y="",x="") + theme(text = element_text(size=12))

leave.out = c("Mean-weight diameter (mm)")
chem.vars = var.labels[-which(var.labels %in% c(text.labs, ag.labs[1:6],leave.out))][c(seq(1:8),16,seq(9,15))]
df.chem = df.int[which(df.int$v %in% chem.vars),]
ggplot(df.chem, aes(y=factor(t, levels=trt.names), x=mean)) + 
  geom_errorbar(aes(xmin=`5`, xmax=`95`), width=0.25, color="black", position=position_dodge(width=0.5)) +
  geom_errorbar(aes(xmin=`25`, xmax=`75`), width=0.25, color="blue", position=position_dodge(width=0.5)) +
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