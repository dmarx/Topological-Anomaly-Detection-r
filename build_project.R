library("devtools")
library(roxygen2)

setwd("C:/Users/davidmarx/Documents/Projects/Toy Projects/Topological-Anomaly-Detection_R")
create("TADr")

setwd('./TADr')
document('.')

setwd('..')
install('TADr')