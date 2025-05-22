path_to_tree_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Public-Repo/Tree_Analysis"
path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Public-Repo"
setwd(path_to_tree_folder)

library(rethinking)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(reshape2)
library(RColorBrewer)

# read in posterior intervals
int.df = read.csv("Posteriors/CarbonRichness_Means_Intervals_5Chains.csv")

# variables
var.order = c("n.total","N.herb","N.tree",
              "MAOM.C","POM.C","SOC","TIC","TC",
              "int.fwd.carbon","cwd.carbon",
              "snag.carbonmin","dead.stem.carbonmin",
              "live.woody.c.stock","live.stem.carbon",
              "FWD.C.Mg.ha","HL.C.Mg.ha","total.C.Mg.ha",
              "total.live.carbon","total.dead.carbon",
              "total.veg.carbon","total.ecosystem.carbon")
var.labels = c("Total richness","Herbaceous species","Tree species",
               "MAOM-C (< 53 \U00B5m)","POM-C (\u2265 53 \U00B5m)",
               "TOC","Carbonates","TC",
               "Fine woody debris (2.5-7.5 cm)","Coarse woody debris (\u2265 7.6 cm)",
               "Standing dead trees (\u2265 2.5 cm)","Dead stems (< 2.5 cm)",
               "Live trees (\u2265 2.5 cm)","Live stems (< 2.5 cm)",
               "Fine woody debris (< 2.5 cm)","Herbaceous litter","Herbaceous biomass",
               "Total live vegetation","Total dead vegetation",
               "Total vegetation","Total ecosystem")

# plot intervals for all variables
ggplot(int.df, aes(y=factor(trt, levels=trt.names), x=mean)) +
       geom_errorbar(aes(xmin=X5, xmax=X95), width=0.25, color="black") +
       geom_errorbar(aes(xmin=X25, xmax=X75), width=0.25, color="blue") +
       geom_point() + labs(y="",x="Posterior mean") +
       facet_wrap(.~factor(var, levels=var.labels), 
                  ncol=3, scales="free_x")

# print intervals
i = 21 # left off on 13
print(cbind(int.df[which(int.df$var == var.labels[i]),c("var","trt")],
            round(int.df[which(int.df$var == var.labels[i]),c("mean","X5","X95")],1)))
print(cbind(int.df[which(int.df$var == var.labels[i]),c("var","trt")],
            round(int.df[which(int.df$var == var.labels[i]),c("mean","X5","X95")],2)))
print(cbind(int.df[which(int.df$var == var.labels[i]),c("var","trt")],
            round(int.df[which(int.df$var == var.labels[i]),c("mean","X5","X95")],3)))

## plot stacked means of biomass variables
SOC.est.df = data.frame(matrix(nrow=3,ncol=2))
SOC.est.df$`IPCC SOC` = c("Natural wetland","Revegetated cropland","Annual crops")
SOC.est.df$value = c(88.0, 77.2, 55.9)

ABC.est.df = data.frame(matrix(nrow=2,ncol=2))
ABC.est.df$`IPCC woody biomass` = c("Restored temperate forest","Mature temperate forest")
ABC.est.df$value = c(46.1, 62.4)

dead.est.df = data.frame(matrix(nrow=2,ncol=2))
dead.est.df$`IPCC dead wood and litter` = c("Restored temperate forest","Mature temperate forest")
dead.est.df$value = c(14.2, 29.7)

TC.est.df = data.frame(matrix(nrow=2,ncol=2))
TC.est.df$`IPCC total organic C` = c("Restored forested wetland","Mature forested wetland")
TC.est.df$value = c(132.5, 180.1)

# total ecosystem
stack.vars1 = c("TIC","Total dead vegetation","TOC","Total live vegetation")
v.palette <- brewer.pal(11,"RdYlGn")[c(7,9,11)]
s.palette <- c(brewer.pal(9,"Greys")[5],
               brewer.pal(11,"BrBG")[c(2,1)])
d.palette <- brewer.pal(9,"YlOrBr")[seq(3,8)]
p1 = ggplot(int.df[which(int.df$var %in% stack.vars1),], 
            aes(y=factor(trt, levels=trt.names), x=mean, 
                fill=factor(var, levels=stack.vars1))) +
            geom_bar(stat="identity",position="stack") + 
            labs(fill="",y="",
                 x="Posterior mean (Mg/ha)",title="") +
            geom_vline(data=TC.est.df, aes(xintercept=value, 
                                           color=`IPCC total organic C`,
                                           linetype=`IPCC total organic C`), 
                       linewidth=1) +
            scale_color_manual(values=c("blue","red")) +
            scale_linetype_manual(values=c("dashed","solid")) +
            scale_fill_manual(values=c(brewer.pal(9,"Greys")[5],d.palette[4],s.palette[3],v.palette[3]),
                              labels=c("Carbonates","Litter and woody debris","Soil organic matter","Living biomass")) +
            theme(text=element_text(size=14), legend.key.size=unit(0.7,'cm'))
p1

#  live vegetation plot
stack.vars2 = c("Herbaceous biomass","Live stems (< 2.5 cm)","Live trees (\u2265 2.5 cm)")
ipcc.vars2 = c("Restored temperate forest","Mature temperate forest")
p2 = ggplot(int.df[which(int.df$var %in% stack.vars2),], 
            aes(y=factor(trt, levels=trt.names), x=mean, 
                fill=factor(var, levels=stack.vars2))) +
            geom_bar(stat="identity",position="stack") + 
            labs(fill="",y="",
                 x="Posterior mean (Mg/ha)",title="") + 
            geom_vline(data=ABC.est.df,
                       aes(xintercept=value, 
                           color=factor(`IPCC woody biomass`,levels=ipcc.vars2),
                           linetype=factor(`IPCC woody biomass`,levels=ipcc.vars2)), 
                       linewidth=1) +
            scale_color_manual(values=c("red","blue")) +
            scale_linetype_manual(values=c("solid","dashed")) + 
            scale_fill_manual(values=c("olivedrab3","olivedrab4",v.palette[3])) +
            guides(color=guide_legend(title="IPCC woody biomass"),
                   linetype=guide_legend(title="IPCC woody biomass")) +
            theme(text=element_text(size=14), legend.key.size=unit(0.7,'cm'))
p2

# debris plot
stack.vars3 = c("Herbaceous litter","Fine woody debris (< 2.5 cm)",
                "Fine woody debris (2.5-7.5 cm)","Coarse woody debris (\u2265 7.6 cm)",
                "Dead stems (< 2.5 cm)","Standing dead trees (\u2265 2.5 cm)")
ipcc.vars3 = c("Restored temperate forest","Mature temperate forest")
p3 = ggplot(int.df[which(int.df$var %in% stack.vars3),], 
            aes(y=factor(trt, levels=trt.names), x=mean, 
                fill=factor(var, levels=stack.vars3))) +
            geom_bar(stat="identity",position="stack") + 
            labs(fill="",y="",
                 x="Posterior mean (Mg/ha)",title="") + 
            geom_vline(data=dead.est.df, 
                       aes(xintercept=value,
                           color=factor(`IPCC dead wood and litter`,levels=ipcc.vars3),
                           linetype=factor(`IPCC dead wood and litter`,levels=ipcc.vars3)), 
                       linewidth=1) +
            scale_color_manual(values=c("red","blue")) +
            scale_linetype_manual(values=c("solid","dashed")) + 
            scale_fill_manual(values=d.palette) + 
            guides(color=guide_legend(title="IPCC dead wood and litter"),
                   linetype=guide_legend(title="IPCC dead wood and litter")) +
            theme(text=element_text(size=14), legend.key.size=unit(0.7,'cm'))
p3

# soil plot
stack.vars4 = c("Carbonates","POM-C (\u2265 53 \U00B5m)","MAOM-C (< 53 \U00B5m)")
ipcc.vars4 = c("Natural wetland","Revegetated cropland","Annual crops")
p4 = ggplot(int.df[which(int.df$var %in% stack.vars4),], 
            aes(y=factor(trt, levels=trt.names), 
                x=mean, fill=factor(var, levels=stack.vars4))) +
            geom_bar(stat="identity",position="stack") + 
            scale_fill_manual(values=s.palette) +
            labs(fill="",y="",
                 x="Posterior mean (Mg/ha)",title="") +
            geom_vline(data=SOC.est.df, 
                       aes(xintercept=value, 
                           color=factor(`IPCC SOC`,levels=ipcc.vars4),
                           linetype=factor(`IPCC SOC`,levels=ipcc.vars4)), 
                       linewidth=1) +
            scale_color_manual(values=c("blue","purple","red")) + 
            scale_linetype_manual(values=c("dashed","dotted","solid")) +
            guides(color=guide_legend(title="IPCC soil organic carbon"),
                   linetype=guide_legend(title="IPCC soil organic carbon")) +
            theme(text=element_text(size=14), legend.key.size=unit(0.7,'cm'))
p4
p.c.stocks = (p4+theme(plot.margin=unit(c(0,32,0,0),"pt"))+p2)/(p3+theme(plot.margin=unit(c(0,2,0,0),"pt"))+p1)
p.c.stocks
setwd(path_to_repo)
ggsave("Figures/Soil_Veg_Ecosystem_Cstocks.jpeg", 
       plot = p.c.stocks, width = 38, height = 20, units="cm")

# richness plot with error bars
stack.vars5 = c("Total richness","Herbaceous species","Tree species")
p.r = ggplot(int.df[which(int.df$var %in% stack.vars5),], 
       aes(y=factor(trt, levels=trt.names), 
           x=mean, 
           fill=factor(var, levels=stack.vars5))) +
       geom_bar(stat="identity", position=position_dodge()) +
       geom_errorbarh(aes(xmin=`X5`, xmax=`X95`,
                          y=factor(trt,levels=trt.names)), 
                       height=0.2, position=position_dodge(0.9)) + 
       labs(fill="",y="",x="Posterior mean count",title="") + 
       scale_fill_manual(values=c("purple","sienna","darkolivegreen4")) +
       theme(text = element_text(size=14))
setwd(path_to_repo)
ggsave("Figures/Richness.jpeg", 
       plot = p.r, width = 20, height = 12, units="cm")



################################################################################
### run PCA on biomass data
library(ggfortify)
library(ggrepel)  
library(ggalt)
library(plyr)

# call prcomp()
bm.pc = bm.df
omit.vars = c("num",
              "H. jap C stock","P. arun C stock","Mixed biomass C stock",
              "FWD C:N ratio","Litter C:N ratio","Herbaceous C:N ratio")
bm.pc = bm.pc[-which(bm.pc$trt=="R"), -which(colnames(bm.pc) %in% omit.vars)]
pca.biomass = prcomp(bm.pc[,-which(colnames(bm.pc) %in% c("trt","trt.full"))], 
                     center=T, scale=T)
#calculate total variance explained by each principal component
var_explained = pca.biomass$sdev^2 / sum(pca.biomass$sdev^2)

#create scree plot
plot(seq(1:length(var_explained)), var_explained, type="b")

# get loadings matrix
CAloadings <- data.frame(Variables = rownames(pca.biomass$rotation), pca.biomass$rotation)

# add PCA scores to the dataset
bm.pc[, c('PC1', 'PC2')] = pca.biomass$x[, 1:2]

# save variable loadings in a separate dataset
rot = as.data.frame(pca.biomass$rotation[, 1:2])
rot$var = rownames(pca.biomass$rotation)

# rescale the loadings to fit nicely within the scatterplot of our data
mult = max(abs(bm.pc[, c('PC1', 'PC2')])) / max(abs(rot[, 1:2])) / 2
rot[, 1:2] = rot[, 1:2] * mult

# if there are many variables to plot, you can play with ggrepel 
var1 = paste("PC1", paste("(", 100*round(var_explained[1], 4), "%)", sep=""), sep=" ")
var2 = paste("PC2", paste("(", 100*round(var_explained[2], 4), "%)", sep=""), sep=" ")
ggplot(data = rot, aes(x=0, y=0, xend=PC1, yend=PC2, label=var)) +
  geom_point(data = bm.pc, aes(x=PC1, y=PC2, color=factor(trt.full, levels=trt.names[1:5])), inherit.aes=FALSE, linewidth=4) +
  geom_segment(color = 'red', arrow = arrow(length = unit(0.03, "npc"))) +
  geom_text_repel(aes(PC1 * 1, PC2 * 1),color = 'red') +
  labs(x = var1, y = var2, color="Treatment") + 
  theme(legend.title = element_text(linewidth=14),
        legend.text = element_text(linewidth=12), 
        axis.title = element_text(linewidth=14), 
        axis.text = element_text(linewidth=10))
barplot((pca.biomass$sdev^2) / sum(pca.biomass$sdev^2))

################################################################################
# linear analysis code

# fit model with reference for each y var
#lm.mods = list()
#for (i in 1:(n.v-2)) {
#  v = c.vars[i]
#  m.v = ulam(alist(y ~ normal(mu, sigma),
#                   log(mu) <- a + b*N,
#                   a ~ dnorm(0, 1),
#                   b ~ dnorm(0, 1),
#                   sigma ~ dexp(1)),
#             data=c.lists[[v]], chains=10, log_lik=T)
#  lm.mods[[v]] = m.v
#}

# make plot of effects
c.var.labs2 = c("Total ecosystem","Bulk soil","Soil organic carbon","Soil inorganic carbon",
                "Total biomass","Woody biomass","Herbaceous biomass",
                "Herbaceous litter","Fine woody debris")
n.v2 = n.v-2
all.b = data.frame(matrix(nrow=0, ncol=6))
all.b.nr = data.frame(matrix(nrow=0, ncol=6))
colnames(all.b) = c("var","mean","25%","75%","5%","95%")
colnames(all.b.nr) = c("var","mean","25%","75%","5%","95%")
for (i in 1:n.v2) {
  v = c.vars[i]
  b.df = data.frame(matrix(nrow=0, ncol=6))
  b.df.nr = data.frame(matrix(nrow=0, ncol=6))
  colnames(b.df) = c("var","mean","25%","75%","5%","95%")
  colnames(b.df.nr) = c("var","mean","25%","75%","5%","95%")
  b.df[1,c("mean","25%","75%")] = precis(lm.mods[[v]], depth=2, prob=0.50)[2, c("mean","25%","75%")]
  b.df.nr[1,c("mean","25%","75%")] = precis(lm.mods.nr[[v]], depth=2, prob=0.50)[2, c("mean","25%","75%")]
  b.df[1,c("5%","95%")] = precis(lm.mods[[v]], depth=2, prob=0.90)[2, c("5%","95%")]
  b.df.nr[1,c("5%","95%")] = precis(lm.mods.nr[[v]], depth=2, prob=0.90)[2, c("5%","95%")]
  b.df[1,"var"] = c.var.labs2[i]
  b.df.nr[1,"var"] = c.var.labs2[i]
  all.b = rbind(all.b, b.df)
  all.b.nr = rbind(all.b.nr, b.df.nr)
}
all.b$Reference = "Included"
all.b.nr$Reference = "Excluded"
b.r.nr = rbind(all.b, all.b.nr)
b.r.nr$Reference = factor(b.r.nr$Reference)
ggplot(b.r.nr, aes(y=factor(var, levels=c.var.labs2), x=mean, shape=factor(Reference))) + 
  geom_vline(xintercept = 0, linetype=1, color = "darkgray", linewidth=0.5) + 
  geom_errorbar(aes(xmin=`5%`,xmax=`95%`), width=.1, color="black", position=position_dodge(width=0.5)) +
  geom_errorbar(aes(xmin=`25%`,xmax=`75%`), width=.1, color="blue", position=position_dodge(width=0.5)) +
  geom_point(position=position_dodge(width = 0.5), linewidth=1.2) +
  labs(x="Effect linewidth", y="") +
  theme(legend.title = element_blank()) + 
  xlim(-2.7,1.6) 

# try plotting distributions
post.b = data.frame(matrix(nrow=0, ncol=2))
post.b.nr = data.frame(matrix(nrow=0, ncol=2))
colnames(post.b) = c("var","value")
colnames(post.b.nr) = c("var","value")
for (i in 1:n.v2) {
  v = c.vars[i]
  b.samples = extract.samples(lm.mods[[v]])$b
  b.samples.nr = extract.samples(lm.mods.nr[[v]])$b
  b.df = data.frame(matrix(nrow=length(b.samples), ncol=2))
  b.df.nr = data.frame(matrix(nrow=length(b.samples.nr), ncol=2))
  colnames(b.df) = c("var","value")
  colnames(b.df.nr) = c("var","value")
  b.df$var = v
  b.df.nr$var = v
  b.df$value = as.numeric(b.samples)
  b.df.nr$value = as.numeric(b.samples.nr)
  post.b = rbind(post.b, b.df)
  post.b.nr = rbind(post.b.nr, b.df.nr)
}
post.b$Reference = "Included"
post.b.nr$Reference = "Excluded"
post.b.all = rbind(post.b, post.b.nr)

ggplot(post.b.all, aes(x=value, fill=factor(Reference))) + 
  geom_density() + facet_wrap(.~var, ncol=3, scales="free_y") +
  geom_vline(xintercept = 0, linetype=1, color = "darkgray", linewidth=0.5)

# plot posterior intervals for C stocks v. species richness
n = 500
df.all = data.frame(matrix(nrow=0, ncol=7))
df.all.nr = data.frame(matrix(nrow=0, ncol=7))
colnames(df.all) = c("var","N","mu","lb.5","ub.95","lb.25","ub.75")
colnames(df.all.nr) = c("var","N","mu","lb.5","ub.95","lb.25","ub.75")
for (i in 1:n.v2) {
  v = c.vars[i]
  N.seq = seq(0.34, 1.65, length.out=n)
  
  df.v = data.frame(matrix(nrow=n, ncol=7))
  df.v.nr = data.frame(matrix(nrow=n, ncol=7))
  colnames(df.v) = c("var","N","mu","lb.5","ub.95","lb.25","ub.75")
  colnames(df.v.nr) = c("var","N","mu","lb.5","ub.95","lb.25","ub.75")
  
  df.v$N = N.seq*var.means[["N.herb"]]
  df.v.nr$N = N.seq*var.means.nr[["N.herb"]]
  df.v$var = c.var.labs2[i]
  df.v.nr$var = c.var.labs2[i]
  
  l.v <- link( lm.mods[[v]] , data=list( N=N.seq ) )
  l.v.nr <- link( lm.mods.nr[[v]] , data=list( N=N.seq ) )
  
  mu.v = exp(apply( l.v , 2 , mean ))*var.means[[v]]
  mu.v.nr = exp(apply( l.v.nr , 2 , mean ))*var.means.nr[[v]]
  ci.50 <- apply( l.v , 2 , HPDI, prob=0.50)
  ci.50.nr <- apply( l.v.nr , 2 , HPDI, prob=0.50)
  ci.90 <- apply( l.v , 2 , HPDI, prob=0.90)
  ci.90.nr <- apply( l.v.nr , 2 , HPDI, prob=0.90)
  
  df.v$mu = mu.v
  df.v$lb.5 = exp(ci.90[1,])*var.means[[v]]
  df.v$ub.95 = exp(ci.90[2,])*var.means[[v]]
  df.v$lb.25 = exp(ci.50[1,])*var.means[[v]]
  df.v$ub.75 = exp(ci.50[2,])*var.means[[v]]
  
  df.v.nr$mu = mu.v.nr
  df.v.nr$lb.5 = exp(ci.90.nr[1,])*var.means.nr[[v]]
  df.v.nr$ub.95 = exp(ci.90.nr[2,])*var.means.nr[[v]]
  df.v.nr$lb.25 = exp(ci.50.nr[1,])*var.means.nr[[v]]
  df.v.nr$ub.75 = exp(ci.50.nr[2,])*var.means.nr[[v]]
  
  df.all = rbind(df.all, df.v)
  df.all.nr = rbind(df.all.nr, df.v.nr)
}
df.all$Reference = "Included"
df.all.nr$Reference = "Excluded"
df.r.nr = rbind(df.all, df.all.nr)
df.r.nr$Reference = factor(df.r.nr$Reference)

ggplot(df.r.nr, aes(x=N, y=mu, color=Reference)) + 
       geom_line() +
       geom_ribbon(aes(ymin=lb.5,ymax=ub.95),alpha=0.1,linewidth=0) +
       #geom_ribbon(aes(ymin=lb.25,ymax=ub.75),alpha=0.1,linewidth=0,fill="blue") +
       facet_wrap(.~factor(var, levels=c.var.labs2), scales="free_y") + xlim(1.8,8.2) + 
       labs(y="C stock (Mg/ha)",x="Understory Species Richness")

par(mfrow=c(1,1))
df.all$var = factor(df.all$var, levels=c.var.labs)
bm.plot = bm.df[,c("trt","N",c.vars)]
colnames(bm.plot)[3:10] = c.var.labs
bm.plot$N = factor(bm.plot$N)
bm.melt = melt(bm.plot, id.variables=c("trt","N"))
ggplot(df.all, aes(x=N, y=mu)) + geom_line() +
       geom_ribbon(aes(ymin=lb.5,ymax=ub.95),alpha=0.1,linewidth=0) +
       geom_ribbon(aes(ymin=lb.25,ymax=ub.75),alpha=0.1,linewidth=0,fill="blue") +
       facet_wrap(.~factor(var, levels=c.var.labs), scales="free_y") + 
       labs(y="C stock (Mg/ha)",x="Understory Species Richness")
