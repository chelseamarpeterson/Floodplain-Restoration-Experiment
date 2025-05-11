setwd('C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch1/Maps')

library(tidyr)
library(ggplot2)

## compare planned and actual tree survey areas

# read data
actual.areas = read.csv('Woody_Plots_Actual_Aug2022.csv', header=T)
planned.areas = read.csv('Woody_Plots_Planned_Aug2022.csv', header=T)

# sort data
actual.areas = actual.areas %>% unite("Plot", Treatment:Plot_Num, sep= "", remove = FALSE)
planned.areas = planned.areas %>% unite("Plot", Treatment:Plot_Num, sep= "", remove = FALSE)

actual.sort.id = sort(actual.areas$Plot, index.return=T)
actual.areas.sort = actual.areas[actual.sort.id$ix,]

planned.sort.id = sort(planned.areas$Plot, index.return=T)
planned.areas.sort = planned.areas[planned.sort.id$ix,]

# calculate percent difference
actual.areas.sort$area_diff = actual.areas.sort$Area_m2 - planned.areas.sort$Shape_Area
actual.areas.sort$percent_diff = (actual.areas.sort$Area_m2 - planned.areas.sort$Shape_Area)/planned.areas.sort$Shape_Area*100

# plots
ggplot(actual.areas.sort, aes(x=Plot, y=percent_diff, color=Treatment)) + geom_point() + 
       labs(y="Percent Difference in Area (Actual v. Planned)")

ggplot(actual.areas.sort, aes(x=Treatment, y=percent_diff)) + geom_point() + 
       labs(y="Percent Difference in Area (Actual v. Planned)")

actual.areas.sort$planned_area = planned.areas.sort$Shape_Area
ggplot(actual.areas.sort, aes(x=Plot, y=Area_m2, color=Treatment)) + geom_point() + 
       geom_line(aes(x=Plot, y=planned_area, group=Treatment))
