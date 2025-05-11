library(pwr)

pwr.t.test(n=8, d=1.32, sig.level=0.05, type="two.sample", alternative="greater")
pwr.t.test(n=8, d=0.66, sig.level=0.05, type="two.sample", alternative="greater")

# simulation

# import skeleton
skeleton <- read.csv("C:/Users/dzaya1/Box Sync/INHS/miscInhsStuff/huang_annie/treatmentSkeleton.csv")
n <- 6

dfr <- skeleton
for(i in 2:n){
	dfr <- rbind(dfr, skeleton)

}

# start with humulus
h <- subset(dfr, spp=="hj")
# monoc, full water, full light
temp.mean <- 20
h$y <- rnorm(nrow(h),temp.mean,temp.mean*0.215)
h$y[h$comp=="comp.p"] <- h$y[h$comp=="comp.p"]*0.75
h$y[h$light=="partial"] <- h$y[h$light=="partial"]*0.45
h$y[h$water=="lo"] <- h$y[h$water=="lo"]*0.8

with(h, tapply(y, list(comp, light, water), mean))

m.h <- lm(y ~ comp + light + water, data=h)
drop1(m.h, test="F")



h$y <- rnorm(nrow(h),temp.mean,temp.mean*0.215)
h$y[h$comp=="comp.p"] <- h$y[h$comp=="comp.p"]*0.75
h$y[h$light=="partial"] <- h$y[h$light=="partial"]*0.45
h$y[h$water=="lo"] <- h$y[h$water=="lo"]*0.8

with(h, tapply(y, list(comp, light, water), mean))

m.h <- lm(y ~ comp + light + water, data=h)
drop1(m.h, test="F")

m.h.int <- lm(y ~ comp + light + water + comp:light, data=h)
drop1(m.h.int, test="F")


## Humulus simulations, with sd estimated from Perry & Galatowitsch
n <- 3
dfr <- skeleton
for(i in 2:n){
	dfr <- rbind(dfr, skeleton)

}
# start with humulus
h <- subset(dfr, spp=="hj")

nsims <- 1000
h.results <- data.frame(spp=rep("hj", nsims), n=n, comp=NA, light=NA, water=NA) 
for(j in 1:nsims){
	
	# extract empty data for one species
	h <- subset(dfr, spp=="hj")
	
	# special rules for each species
	temp.mean <- 20
	h$y <- rnorm(nrow(h),temp.mean,temp.mean*0.215)
	h$y[h$comp=="comp.p"] <- h$y[h$comp=="comp.p"]*0.75
	h$y[h$light=="partial"] <- h$y[h$light=="partial"]*0.45
	h$y[h$water=="lo"] <- h$y[h$water=="lo"]*0.8
	
	# anova
	m.h <- lm(y ~ comp + light + water, data=h)
	# extract pvalues
	p.vals <- coef(summary(m.h))[,4]
	# store pvalues
	h.results[j,"comp"] <- p.vals[2]
	h.results[j,"light"] <- p.vals[3]
	h.results[j,"water"] <- p.vals[4]

}

sum(h.results$comp < 0.05)
sum(h.results$light < 0.05)
sum(h.results$water < 0.05)


## Repeat, but with double the std dev
n <- 5
nsims <- 1000
h.results <- data.frame(spp=rep("hj", nsims), n=n, comp=NA, light=NA, water=NA) 
for(j in 1:nsims){
	
	# extract empty data for one species
	h <- subset(dfr, spp=="hj")
	
	# special rules for each species
	temp.mean <- 20
	h$y <- rnorm(nrow(h),temp.mean,temp.mean*0.43)
	h$y[h$comp=="comp.p"] <- h$y[h$comp=="comp.p"]*0.75
	h$y[h$light=="partial"] <- h$y[h$light=="partial"]*0.45
	h$y[h$water=="lo"] <- h$y[h$water=="lo"]*0.8
	
	# anova
	m.h <- lm(y ~ comp + light + water, data=h)
	# extract pvalues
	p.vals <- coef(summary(m.h))[,4]
	# store pvalues
	h.results[j,"comp"] <- p.vals[2]
	h.results[j,"light"] <- p.vals[3]
	h.results[j,"water"] <- p.vals[4]

}

sum(h.results$comp < 0.05)
sum(h.results$light < 0.05)
sum(h.results$water < 0.05)


