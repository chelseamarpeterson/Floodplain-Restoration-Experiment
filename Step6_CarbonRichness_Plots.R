path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Public-Repo"
setwd(path_to_tree_folder)

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
