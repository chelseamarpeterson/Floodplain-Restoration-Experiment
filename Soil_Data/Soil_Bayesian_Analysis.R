setwd("C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/SoilData")

library(rethinking)
library(dplyr)
library(ggplot2)
library(reshape2)

# treatment names
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

# read in soil data
soil.dat = read.csv("All_Soil_Data_2023.csv", header=T)

# labels
ag.cols = c("fpom","g475mm","g2mm","g250um","g53um","l53um","mwd")
meq.cols = c("k.meq","ca.meq","mg.meq","h.meq")

# update column names
colnames(soil.dat) = tolower(colnames(soil.dat))
colnames(soil.dat) = c("plot","trt","trt.full","num","quad",
                       "temp","gmoist","vmoist","bd",
                       "n","c","cn","sand","silt","clay",
                       "ph","p","k","ca","mg","som",
                       "no3","nh4","cec",meq.cols,"pom",ag.cols) 
soil.dat$gmoist = soil.dat$gmoist * 100

# make index for treatment, quadrat, and plot number
soil.dat$tid = as.numeric(factor(soil.dat$trt))
soil.dat$qid = as.numeric(factor(soil.dat$quad))
soil.dat$pid = rep(0, nrow(soil.dat))
k = 1
for (t in 1:6) {
  for (p in 1:3) {
    tp.ind = which(soil.dat$tid == t & soil.dat$num == p)
    soil.dat$pid[tp.ind] = k
    k = k + 1
  }
}

# isolate fpom data
fpom.dat = soil.dat[-which(is.na(soil.dat$fpom)),c("plot","trt","trt.full","num","quad","fpom")]

# make index for treatment, quadrat, and plot number
fpom.dat$tid = as.numeric(factor(fpom.dat$trt))
fpom.dat$qid = as.numeric(factor(fpom.dat$quad))
fpom.dat$pid = rep(0, nrow(fpom.dat))
k = 1
for (t in 1:6) {
  for (p in 1:3) {
    tp.ind = which(fpom.dat$tid == t & fpom.dat$num == p)
    fpom.dat$pid[tp.ind] = k
    k = k + 1
  }
}

# make simple list for Bayes analysis
vars = c("temp","gmoist","bd","n","c","cn","sand","silt","clay",
         "ph","p","k","ca","mg","som","pom","no3","nh4",
         meq.cols[1:3],"cec",ag.cols)
n.v = length(vars)
all.lists = list()
all.means = list()
for (i in 1:n.v) {
  v = vars[i]
  if (v == "fpom") {
    all.means[[v]] = mean(fpom.dat[,v])
    d.v =  list(tid=as.integer(fpom.dat$tid), 
                pid=as.integer(fpom.dat$num), 
                y=fpom.dat[,v])
  } else {
    all.means[[v]] = mean(soil.dat[,v])
    d.v =  list(tid=as.integer(soil.dat$tid), 
                pid=as.integer(soil.dat$num), 
                y=soil.dat[,v])
  }
  d.v$y = log(d.v$y/mean(d.v$y))
  all.lists[[v]] = d.v    
}

# plot histograms of transformed variables
#par(mfrow=c(5,3), mar=c(2.5,3,1.5,1))
#for (i in 1:n.v) {
#  v = vars[i]
  #hist(log(all.lists1[[v]]$y/mean(all.lists1[[v]]$y)), main=v, xlab="")
#  hist(all.lists1[[v]]$y, main=v, xlab="")
#}

# model code w/o random effect
#set.seed(14)
#soil.models1 = list()
#soil.models1[[vars[1]]] = model1
#for (i in 2:n.v) {
#  v = vars[i]
#  m.v1 = ulam(alist(y ~ normal(mu, sigma),
#                    mu <- a[tid],
#                    a[tid] ~ dnorm(0, 1),
#                    sigma ~ dexp(1)), 
#              data=all.lists1[[v]], chains=5, log_lik=T)
#  soil.models1[[v]] = m.v1
#}

# model code w/ random effect
soil.models2 = list()
for (i in 1:n.v) {
  v = vars[i]
  m.v = ulam(alist(y ~ normal(mu, sigma),
                    mu <- a[tid] + b[tid, pid],
                    a[tid] ~ dnorm(0, 1),
                    matrix[tid, pid]: b ~ dnorm(0, sigma_b),
                    sigma_b ~ dexp(1),
                    sigma ~ dexp(1)),
              data=all.lists[[v]], chains=10, log_lik=T)
  soil.models2[[v]] = m.v
}

text.list = list(tid=as.integer(soil.dat$tid), 
                 pid=as.integer(soil.dat$num), 
                 y_clay=soil.dat[,"clay"],
                 y_silt=soil.dat[,"silt"],
                 y_sand=soil.dat[,"sand"])
m.text = ulam(alist(y_clay ~ normal(mu_clay, sigma_clay),
                    y_silt ~ normal(mu_silt, sigma_silt),
                    y_sand ~ normal(mu_sand, sigma_sand),
                    mu_silt <- a_silt[tid] + b_silt[tid, pid],
                    mu_sand <- a_sand[tid] + b_sand[tid, pid],
                    mu_clay <- a_clay[tid] + b_clay[tid, pid],
                    mu_clay <- 100 - mu_sand - mu_silt,
                    a_clay[tid] ~ dnorm(0, 1),
                    a_silt[tid] ~ dnorm(0, 1),
                    a_sand[tid] ~ dnorm(0, 1),
                    matrix[tid, pid]: b_silt ~ dnorm(0, sigma_b_silt),
                    matrix[tid, pid]: b_sand ~ dnorm(0, sigma_b_sand),
                    matrix[tid, pid]: b_clay ~ dnorm(0, sigma_b_clay),
                    sigma_b_clay ~ dexp(1),
                    sigma_b_silt ~ dexp(1),
                    sigma_b_sand ~ dexp(1),
                    sigma_clay ~ dexp(1),
                    sigma_silt ~ dexp(1),
                    sigma_sand ~ dexp(1)),
           data=text.list, chains=1, log_lik=T)
precis(m.text, depth=2)
precis(soil.models2[["silt"]], depth=2)

# variable labels
var.new = c("som","pom","c","n","cn","temp",
            "gmoist","bd","sand","silt","clay",
            "no3","nh4","p","k","ca","mg",
            meq.cols[1:3],"cec","ph",ag.cols)
var.labels = c("SOM (%)","POM (%)","TC (%)","TN (%)","C:N Ratio","Temperature (C)",
               "Gravitational moisture (%)","Bulk density (g/cm3)",
               "Sand (%)","Silt (%)","Clay (%)",
               "NO3 (ppm)","NH4 (ppm)","P (ppm)",
               "K (ppm)","Ca (ppm)","Mg (ppm)",
               "K (meq)","Ca (meq)","Mg (meq)",
               "CEC (meq/100 g)","pH","Floating POM",
               ">4.75 mm","2 - 4.75 mm","250 \u03bcm - 2 mm","53 - 250 \u03bcm",
               "<53 \u03bcm","Mean-weight diameter (mm)")

#### compare WAIC and PSIS for all soil variables
#df.waic = data.frame(matrix(nrow=0, ncol=8))
#df.psis = data.frame(matrix(nrow=0, ncol=8))
#colnames(df.waic)[1:7] = c(colnames(compare(soil.models2[[vars[1]]], soil.models1[[vars[1]]], func="WAIC")),"var")
#colnames(df.psis)[1:7] = c(colnames(compare(soil.models2[[vars[1]]], soil.models1[[vars[1]]], func="PSIS")),"var")
#colnames(df.waic)[8] = "model"
#colnames(df.psis)[8] = "model"
#for (i in 1:n.v) {
#  v = var.new[i]
#  m1 = soil.models1[[v]]
#  m2 = soil.models2[[v]]
#  v.waic = compare(m1, m2, func="WAIC")
#  v.psis = compare(m1, m2, func="PSIS")
#  v.waic$var = var.labels[i]
#  v.psis$var = var.labels[i]
#  v.waic$model[which(row.names(v.waic) == "m1")] = "No random effects"
#  v.waic$model[which(row.names(v.waic) == "m2")] = "Random effects"
#  v.psis$model[which(row.names(v.psis) == "m1")] = "No random effects"
#  v.psis$model[which(row.names(v.psis) == "m2")] = "Random effects"
#  df.waic = rbind(df.waic, v.waic)
#  df.psis = rbind(df.psis, v.psis)
#}
#df.waic$metric = "Widely Applicable Information Criterion (WAIC)"
#df.psis$metric = "Pareto-smoothed Importance Sampling Cross-validation (PSIS))"
#df.all.met = rbind(df.waic, df.psis)

#mod.labs = c("No random effects","Random effects")
#df.waic$model = factor(df.waic$model, levels=mod.labs)
#p1 = ggplot(df.waic[which(df.waic$model=="Random effects"),], aes(y=factor(var, levels=var.labels), x=1-weight)) + 
#      geom_point() + xlim(c(0,1)) +
#      theme(axis.title.y=element_blank(), legend.position="none") + 
#      labs(x = "1 - (Random effects weight)") + scale_y_discrete(limits = rev)
#p2 = ggplot(df.waic, aes(y=factor(var, levels=var.labels), x=WAIC, color=model)) + 
#            geom_point(position=position_dodge(width=0.5)) + 
#            geom_errorbar(aes(y=var, xmin=WAIC-SE, xmax=WAIC+SE), width=.1, position=position_dodge(width=0.5)) +
#            theme(axis.title.y=element_blank(), legend.title=element_blank(), axis.text.y=element_blank()) + 
#            labs(x = "WAIC") + scale_y_discrete(limits = rev)
#p1 + p2
#p1 = ggplot(df.psis, aes(y=factor(var, levels=var.labels), x=PSIS, color=factor(model, levels=mod.labs))) + 
#            geom_point(position=position_dodge(width=0.5)) + 
#            geom_errorbar(aes(y=var, xmin=PSIS-SE, xmax=PSIS+SE), width=.1, position=position_dodge(width=0.5)) +
#            theme(axis.title.y=element_blank(), legend.title=element_blank(), axis.text.y=element_blank()) + 
#            labs(x = "PSIS") + scale_y_discrete(limits = rev)
#p2 = ggplot(df.waic[which(df.waic$model=="Random effects"),], aes(y=factor(var, levels=var.labels), x=1-weight)) + 
#            geom_point() + xlim(c(0,1)) +
#            theme(axis.title.y=element_blank(), legend.position="none") + 
#            labs(x = "1- (Random effects weight)") + scale_y_discrete(limits = rev)
#p2 + p1

### plot posteriors

# extract posterior samples
all.samples = list()
all.a = list()
all.sig = list()
for (i in 1:n.v) {
  v = vars[i]
  all.samples[[v]] = extract.samples(soil.models2[[v]])
  all.a[[v]] = exp(all.samples[[v]]$a) * all.means[[v]]
  all.sig[[v]] = exp(all.samples[[v]]$sigma) * all.means[[v]]
}

# get posterior intervals and write to file
df.int = data.frame(matrix(nrow=0, ncol=10))
colnames(df.int) = c("v","t","mean","sd","5","95","15","85","25","75")
for (i in 1:n.v) {
  df.v = data.frame(matrix(nrow=6, ncol=10))
  colnames(df.v) = c("v","t","mean","sd","5","95","15","85","25","75")
  v = var.new[[i]]
  a.means = apply(all.a[[v]], 2, mean)
  a.sds = apply(all.a[[v]], 2, sd)
  a.hpdi.90 = apply(all.a[[v]], 2, HPDI, prob=0.90)
  a.hpdi.70 = apply(all.a[[v]], 2, HPDI, prob=0.70)
  a.hpdi.50 = apply(all.a[[v]], 2, HPDI, prob=0.50)
  df.v[1:6,"v"] = var.labels[i]
  df.v[1:6,"t"] = trt.names
  df.v[1:6,"mean"] = a.means
  df.v[1:6,"sd"] = a.sds
  df.v[1:6,"5"] = a.hpdi.90[1,]
  df.v[1:6,"95"] = a.hpdi.90[2,]
  df.v[1:6,"15"] = a.hpdi.70[1,]
  df.v[1:6,"85"] = a.hpdi.70[2,]
  df.v[1:6,"25"] = a.hpdi.50[1,]
  df.v[1:6,"75"] = a.hpdi.50[2,]
  df.int = rbind(df.int, df.v)
}
write.csv(df.int, "All_Soil_Posterior_Intervals_ChainCount10.csv", row.names=F)
cbind(df.int[,c("v","t")],round(df.int[,c("mean","5","95")],1))

### plotting

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
p.ag = ggplot(df.stack.ag, aes(y=factor(t, levels=trt.names), x=mean, 
                               fill=factor(v, levels=ag.labs))) + 
              geom_bar(stat="identity",position="stack") +
              geom_text(aes(label = round(mean,1)), 
                        nudge_x = df.stack.ag$position,
                        color="white") +
              labs(x="Posterior mean of relative mass (%)",
                   y="",fill="Aggregate size class") +
              theme(text = element_text(size=14))
              #theme(axis.text.y = element_blank())
p.ag

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
                aes(y=factor(t, levels=trt.names), x=mean, fill=factor(v, levels=text.labs.new))) + 
                geom_bar(stat = "identity", position = "stack") +
                geom_text(aes(label = round(mean,1)), 
                          nudge_x = df.stack.text$position,
                          color = "white") +
                labs(x="Posterior mean of relative mass (%)",y="", 
                     fill="Particle size class") +
                scale_fill_manual(values=c("burlywood3","seashell4","coral3"))  +
                theme(text = element_text(size=14))
                #geom_errorbar(aes(xmin=lower, xmax=upper), width=0.1)
p.text  

library(patchwork)
p.text + p.ag

# compare organic and total carbon
c.labs = c("SOM (%)","TC (%)")
df.stack.c = df.int[which(df.int$v %in% c.labs),]
omit.cols = c("v","t","sd")
omit.inds = which(colnames(df.stack.c) %in% omit.cols)
df.stack.c[which(df.stack.c$v == "SOM (%)"),-omit.inds] = 0.5*df.stack.c[which(df.stack.c$v == "SOM (%)"),-omit.inds]
for (i in 1:6) {
  df.trt = df.stack.c[which(df.stack.c$t == trt.names[i]),]
  df.tic = data.frame(matrix(nrow=1, ncol=length(df.trt)))
  colnames(df.tic) = colnames(df.trt)
  df.tic$v = "Inorganic (%)"
  df.tic$t = trt.names[i]
  df.tc = df.trt[which(df.trt$v == "TC (%)"),]
  df.soc = df.trt[which(df.trt$v == "SOM (%)"),]
  df.tic[1,c("mean","5","95","15","85","25","75")] = df.tc[1,c("mean","5","95","15","85","25","75")] - df.soc[1,c("mean","5","95","15","85","25","75")]
  df.stack.c = rbind(df.stack.c, df.tic)
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

## plot continuous distributions
new.df = all.a[["som"]]
colnames(new.df) = trt.names
all.melt = melt(new.df, id.vars=trt.names)
all.melt$var = "SOM (%)"
all.melt = all.melt[,2:4]
colnames(all.melt) = c("trt","value","var")
for (i in 2:n.v) {
  v = var.new[i]
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

## test simple hierarchical models
model3.code = "
data{
     vector[90] y;
     array[90] int pid;
     array[90] int tid;
}
parameters{
     matrix[6,3] b;
     vector[6] a;
     real<lower=0> sigma_t;
     real<lower=0> sigma;
}
model{
     vector[90] mu;
     matrix[6,3] mu_t;
     sigma_t ~ exponential( 1 );
     sigma ~ exponential( 1 );
     a ~ normal( 0 , 1 );
     for ( i in 1:6 ) {
        for (j in 1:3) {
           mu_t[i, j] = a[i];
        }
     }
     to_vector( b ) ~ normal( to_vector(mu_t) , sigma_t);
     for ( i in 1:90 ) {
         mu[i] = b[tid[i], pid[i]];
     }
     y ~ normal( mu , sigma );
}"

stan_test <- stan( model_code=model3.code , data=all.lists[["som"]] , chains=1)
precis(stan_test, depth=3)
compare( stan_test , model1 )


