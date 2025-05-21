library(ggplot2)
library(ggfortify)
library(patchwork)

## read in all soil data
soil.dat = read.table("All_Soil_Data_2023.csv", header=T, sep=",")
colnames(soil.dat) = tolower(colnames(soil.dat))

# treatment names
trts = c("A","B","C","D","E","R")
trt.names = c("Balled-and-burlapped","Bareroot","Seedling","Acorn","Seedbank","Reference")

# labels
ag.cols = c("fpom","g475mm","g2mm","g250um","g53um","l53um","mwd")
meq.cols = c("k.meq","ca.meq","mg.meq","h.meq")

# update column names
colnames(soil.dat) = tolower(colnames(soil.dat))
colnames(soil.dat) = c("plot","trt","trt.full","num","quad",
                       "temp","gmoist","vmoist","bd",
                       "n","c","cn","sand","silt","clay",
                       "ph","p","k","ca","mg","som",
                       "no3","nh4","cec",
                       meq.cols,"pom",ag.cols)

# MWD & no ref
soil.dat$soc = 0.5*soil.dat$som
pc.dat = soil.dat[-which(soil.dat$trt=="R"),-which(colnames(soil.dat) %in% c(ag.cols[1:6],"som","pom","vmoist",meq.cols))]
colnames(pc.dat)[6:24] = c("Temp","Moisture","BD","TN","TC","CN",
                           "Sand","Silt","Clay","pH","P","K","Ca",
                           "Mg","NO3","NH4","CEC","MWD","SOC")
pca.all = prcomp(~ SOC+TC+TN+CN+MWD+CEC+Sand+Silt+Clay+pH+P+K+Ca+Mg+NO3+NH4+BD+Moisture+Temp, 
                 data=pc.dat, scale=T, center=T)
pc.dat$trt.full = factor(pc.dat$trt.full, levels=trt.names)
pc.dat$num = factor(pc.dat$num)
autoplot(pca.all, x=1, y=2, data=pc.dat, 
         loadings=T, loadings.label=T, size=4, loadings.label.colour='black',
         col="trt.full", shape="num", loadings.label.size=4, 
         loadings.label.vjust = -0.15, loadings.label.hjust = 0.1) +
         labs(color="Treatment", shape="Plot")

autoplot(pca.all, x=2, y=3, data=all.dat[-which(all.dat$trt=="R"),], 
         loadings=T, loadings.label=T, size=4, loadings.label.colour='black',
         col="trt.full", shape="num", loadings.label.size=5, 
         loadings.label.vjust = -0.15, loadings.label.hjust = 0.1) +
         labs(color="Treatment", shape="Plot")
p1 + p2
summary(pca.all)
pca.all$rotation[,1:5]

sort(pca.all$rotation[,1]) # Silt/BD - SOM/CEC
sort(pca.all$rotation[,2]) # Sand - Silt - Clay/P
sort(pca.all$rotation[,3]) # Clay/sand - CEC/Silt

# Extract the eigenvalues from the PCA object
eigenvalues <- pca.all$sdev^2

# Percentage of variance explained
plot(eigenvalues/sum(eigenvalues), type = "b",
     xlab = "Principal Component",
     ylab = "Percentage of Variance Explained", ylim=c(0,0.2))
#abline(v = 2, col = "red")


round(eigenvalues/sum(eigenvalues)*100,2)
round(cumsum(eigenvalues/sum(eigenvalues)*100),2)

################################################################################
### simple Bayes models: SOM ~ All other vars ----

# put simple dfs into lists
lm.lists = list()
xvars = c("MWD","g475mm","g2mm","g250um","g53um","l53um",     
          "temp","gmoist","bd","sand","silt","clay","ph",
          "no3","nh4","p","k","ca","mg","cec")
xvar.labs = c("MWD",">4.75mm",">2mm",">250um",">53um","<53um",     
              "Temperature","Gravimetric moisture","Bulk density",
              "Sand","Silt","Clay","pH",
              "NO3","NH4","P","K","Ca","Mg","CEC")
n.xv = length(xvars)
all.dat$tid = as.numeric(factor(all.dat$trt))
all.dat$pid = as.numeric(factor(all.dat$num))
for (i in 1:n.xv) {
  x = xvars[i]
  all.x = as.numeric(unlist(all.dat[,x]))
  d.x = list(tid = as.integer(all.dat$tid),
             pid = as.integer(all.dat$pid),
             som = log(all.dat$som/mean(all.dat$som)),
             x = log(all.x/mean(all.x)))
  hist(log(all.x/mean(all.x)), main=x)
  lm.lists[[x]] = d.x
}

# fit model for each x var
lm.mods1 = list()
for (i in 3:n.xv) {
  x = xvars[i]
  m.x = ulam(alist(som ~ normal(mu, sigma),
                   mu <- a + ga[tid] + (b + gb[tid])*x,
                   a ~ dnorm(0, 1),
                   ga[tid] ~ dnorm(0, 1),
                   b ~ dnorm(0, 1),
                   gb[tid] ~ dnorm(0, 1),
                   sigma ~ dexp(1)),
             data=lm.lists[[x]], chains=1, log_lik=T)
  lm.mods1[[x]] = m.x
}

# evaluate model goodness of fit
som.df.waic = compare(lm.mods1$g475mm,lm.mods1$g2mm,lm.mods1$g250um,lm.mods1$g53um,lm.mods1$l53um,lm.mods1$MWD,
                      lm.mods1$temp,lm.mods1$gmoist,lm.mods1$ph,lm.mods1$bd,lm.mods1$sand,lm.mods1$silt,lm.mods1$clay,
                      lm.mods1$no3,lm.mods1$nh4,lm.mods1$p,lm.mods1$k,lm.mods1$ca,lm.mods1$mg,lm.mods1$cec)
som.df.waic$var = c("k","mg","nh4","p","bd","ph","temp","silt","sand","cec","no3","ca","gmoist",
                    "g250um","g475mm","l53um","g2mm","clay","MWD","g53um")
som.df.waic$var = factor(som.df.waic$var, levels=som.df.waic$var)
ggplot(som.df.waic, aes(y=var, x=WAIC)) + 
       geom_point() + 
       geom_errorbar(aes(y=var, xmin=WAIC-SE, xmax=WAIC+SE), width=.1, position=position_dodge(width=0.5))

# make plot of effects
all.ga = data.frame(matrix(nrow=0, ncol=7))
all.gb = data.frame(matrix(nrow=0, ncol=7))
colnames(all.ga) = c("mean","25%","75%","5%","95%","trt","var")
colnames(all.gb) = c("mean","25%","75%","5%","95%","trt","var")

for (i in 1:n.xv) {
  ga.df.50 = precis(lm.mods1[[xvars[i]]], depth=2, prob=0.50)[1:7,]
  ga.df.50[2:7,"mean"] = ga.df.50[1,"mean"] + ga.df.50[2:7,"mean"]
  ga.df.90 = precis(lm.mods1[[xvars[i]]], depth=2, prob=0.90)[1:7,]
  ga.df.90[2:7,"mean"] = ga.df.90[1,"mean"] + ga.df.90[2:7,"mean"]
  gb.df.50 = precis(lm.mods1[[xvars[i]]], depth=2, prob=0.50)[8:14,]
  gb.df.50[2:7,"mean"] = gb.df.50[1,"mean"] + gb.df.50[2:7,"mean"]
  gb.df.90 = precis(lm.mods1[[xvars[i]]], depth=2, prob=0.90)[8:14,]
  gb.df.90[2:7,"mean"] = gb.df.90[1,"mean"] + gb.df.90[2:7,"mean"]
  ga.df.90$var = xvar.labs[i]
  gb.df.90$var = xvar.labs[i]
  ga.df.90$trt = c("Overall",trt.names)
  gb.df.90$trt = c("Overall",trt.names)
  ga.df = cbind(ga.df.50[,c("mean","25%","75%")], ga.df.90[,c("5%","95%","trt","var")])
  gb.df = cbind(gb.df.50[,c("mean","25%","75%")], gb.df.90[,c("5%","95%","trt","var")])
  all.ga = rbind(all.ga, ga.df)
  all.gb = rbind(all.gb, gb.df)
}

all.gb$trt = factor(all.gb$trt, levels=c("Overall",trt.names))
ggplot(all.gb, aes(y=trt, x=mean)) + 
       geom_point() +
       geom_errorbar(aes(xmin=`5%`,xmax=`95%`), width=.3, color="black") +
       geom_errorbar(aes(xmin=`25%`,xmax=`75%`), width=.3, color="blue") +
       facet_wrap(.~factor(var, levels=xvar.labs), scales="free_x") + 
       xlim(c(-2.18,2.18)) +
       labs(x="Effect size", y="") +
       theme(legend.title = element_blank()) 

# plot results v. som
n = 200
df.all = data.frame(matrix(nrow=0, ncol=6))
colnames(df.all) = c("var","x","trt","mu","lb","ub")
for (j in 1:n.xv) {
  x = xvars[j]
  x.seq = seq(min(lm.lists[[x]]$x), max(lm.lists[[x]]$x), length.out=n)
  for (i in 1:6) {
    df.t = data.frame(matrix(nrow=n, ncol=6))
    colnames(df.t) = c("var","x","trt","mu","lb","ub")
    df.t$var = xvar.labs[j]
    df.t$x = x.seq
    df.t$trt = trt.names[i]
    l.t <- link( lm.mods1[[x]] , data=list( tid=rep(i, n), x=x.seq ) )
    mu.t = apply( l.t , 2 , mean )
    ci.t <- apply( l.t , 2 , HPDI, prob=0.50)
    df.t$mu = mu.t
    df.t$lb = ci.t[1,]
    df.t$ub = ci.t[2,]
    df.all = rbind(df.all, df.t)
  }
  post <- extract.samples(lm.mods1[[x]], 1000)
  df.net = data.frame(matrix(nrow=n, ncol=6))
  colnames(df.net) = c("var","x","trt","mu","lb","ub")
  df.net$var = xvar.labs[j]
  df.net$x = x.seq
  df.net$trt = "Overall"
  for (k in 1:n) {
    net.pred = post$a + post$b*x.seq[k]
    df.net[k,"mu"] = apply( net.pred , 2 , mean )
    ci.net <- apply( net.pred , 2 , HPDI, prob=0.50)
    df.net[k,"lb"] = as.numeric(ci.net[1,])
    df.net[k,"ub"] = as.numeric(ci.net[2,])
  }
  df.all = rbind(df.all, df.net)
}
sort(unique(df.all$var))
par(mfrow=c(1,1))
df.all$var = factor(df.all$var, levels=xvar.labs)
df.all$trt = factor(df.all$trt, levels=c(trt.names,"Overall"))
ggplot(df.all, aes(x=x, y=mu, fill=trt, color=trt)) + geom_line() +
       geom_ribbon(aes(ymin=lb,ymax=ub),alpha=0.1,linewidth=0) +
       facet_wrap(.~var, scales="free")

# SOM v. aggregate
ag.vars = c("MWD","g475mm","g2mm","g250um","g53um","l53um")
ag.df = df.all[which(df.all$var %in% ag.vars),]
ag.df$trt = factor(ag.df$trt, levels=trt.names)
ggplot(ag.df, aes(x=x, y=mu, fill=trt, color=trt)) + geom_line() +
       geom_ribbon(aes(ymin=lb,ymax=ub),alpha=0.1,linewidth=0) +
       facet_wrap(.~var, scales="free", ncol=3) +
       theme(legend.title = element_blank())

# SOM v. all other variables
nonag.df = df.all[-which(df.all$var %in% ag.vars),]
nonag.df$trt = factor(nonag.df$trt, levels=trt.names)
ggplot(nonag.df, aes(x=x, y=mu, fill=trt, color=trt)) + geom_line() +
  geom_ribbon(aes(ymin=lb, ymax=ub),alpha=0.1,linewidth=0) +
  facet_wrap(.~var, scales="free", ncol=4) +
  theme(legend.title = element_blank())












