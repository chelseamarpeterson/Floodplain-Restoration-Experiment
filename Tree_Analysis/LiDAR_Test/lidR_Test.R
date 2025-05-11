setwd("C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Maps/LIDAR")

library(lidR)
library(RCSF)

# display map
las <- readLAS("SW_Henry_2020.las")
plot(las)

# Khosravipour et al. pitfree algorithm
thr <- c(0,2,5,10,15)
edg <- c(0, 1.5)
chm <- rasterize_canopy(las, 1, pitfree(thr, edg))

plot(chm)

# individual tree segmentation
las <- segment_trees(las, li2012())
col <- random.colors(200)
plot(las, color = "treeID", colorPalette = col)

# calculate metrics
hmean <- pixel_metrics(las, ~max(Z), 10) # calculate mean at 10 m
plot(hmean, col = height.colors(50))


las <- classify_ground(las, algorithm = csf())
gnd <- filter_ground(las)
plot(gnd, size = 3, bg = "white") 
