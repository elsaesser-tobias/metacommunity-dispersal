
library(ggplot2);
library(deSolve);
library(simecol);
library(OCNet);
library(reticulate);
library(readxl);
library(pls);

library(openxlsx);
library(dplyr);

#clear enviroment
rm(list = ls())
#setwd("C:/Users/elsae/Desktop/Server")
start_time <- Sys.time()
#b will later be used to set different random seeds in each loop
b<- as.numeric(read.table("numbers.txt"))

#track time

#set random seed for reproducibility
set.seed(b)

#set landscape distribution
urban_landscape_percent <- 0.33
forest_landscape_percent <- 0.33
agriculture_landscape_percent <- 0.34

if(urban_landscape_percent + forest_landscape_percent + agriculture_landscape_percent != 1){
  stop("Not 100% of landscape assigned")
}

#in step one the OCN will be generated, habitats and migration routes will be calclulated
#the second step will then be the ODE model of the metacommunity

#size of the OCN will be OCN_size^2, cellsz determines the size of a single cell of the OCN matrix in meters
OCN_size <- 100
cellsz <- 100
#percentage of cells of the OCN that will later be counted as stream, the rest is considered drainage area
percentage_of_stream_points <- 0.09
number_of_stream_points <- (OCN_size * OCN_size) * percentage_of_stream_points

#percentage of the stream points that get used as as habitat
#percentage_habitat_points <- 0.02
nhabitats <- 20

#set dispersal parameters
hab_radius <- 100
n_steps <- 5000
n_dragonflies <- 1000

#set ODE parameters
dispersal_percent <- 0.45





#create the optimal channel network -------------------------------------------------------
OCN <- create_OCN(OCN_size, OCN_size, outletPos = OCN_size / 2, cellsize = cellsz)
OCN <- landscape_OCN(OCN, zMin = 50000)
OCN2 <- OCN

#calculate thresholds
threshold <- find_area_threshold_OCN(OCN)

#find threshold_area so that the number of stream points is slightly less than number_of_stream_points
threshold_area <- 0
for(i in 1:length(threshold$nNodesRN)){
  if(threshold$nNodesRN[i] < number_of_stream_points){
    threshold_area <- threshold$thrValues[i]
    break
  }
}

#aggregate OCN to reduce stream cells, use river network (RN) information from aggregated OCN
OCN <- aggregate_OCN(landscape_OCN(OCN, zMin = 50000), thrA = threshold_area)



#problem: no stream order for each RN node
#get AG-stream order to RN where available
RNSO <- c()

for(i in 1:length(OCN$AG$toRN)){
  RNSO[OCN$AG$toRN[i]] <- OCN$AG$streamOrder[i]
}
#RNSO[which(OCN$RN$downNode == 0)] <- max(RNSO)
#fill the rest through downstream node
#not well optimized but works

#the code takes all the places where the above node is available, but not the node itself
#the streamorder of the node itself is then filled with the one above
#this works because at AG level all changes in stream order are accounted for
for(j in 1:length(OCN$RN$downNode)){
  for(i in 1:length(OCN$RN$downNode)){
    if(is.na(RNSO[i]) == FALSE && is.na(RNSO[OCN$RN$downNode[i]]) && OCN$RN$downNode[i] != 0){
      RNSO[OCN$RN$downNode[i]] <- RNSO[i]
    }
  }
}
#Problem: outlet ist na, deshalb stoppt es nicht?



RNX <- OCN$RN$X
RNY <- OCN$RN$Y
RNDN <- OCN$RN$downNode


#add habitats ----------------------------------------------------------
#nhabitats <- round(length(RNX) * percentage_habitat_points)
set.seed(b)
habitat_points <- sample(1:length(RNX), nhabitats)

habitatsx <- RNX[habitat_points]
habitatsy <- RNY[habitat_points]



#export relevant data to files, so that they can be read by the python script -------------------------------
write.table(RNX, "x_points.txt", row.names = FALSE, col.names = FALSE)
write.table(RNY, "y_points.txt", row.names = FALSE, col.names = FALSE)
write.table(RNSO, "SO.txt", row.names = FALSE, col.names = FALSE)
write.table(RNDN, "DN.txt", row.names = FALSE, col.names = FALSE)
write.table(b, "py_seed.txt", row.names = FALSE, col.names = FALSE)
write.table(cellsz, "py_cellsz.txt", row.names = FALSE, col.names = FALSE)
write.table(habitatsx, "habitatsx.txt", row.names = FALSE, col.names = FALSE)
write.table(habitatsy, "habitatsy.txt", row.names = FALSE, col.names = FALSE)

write.table(urban_landscape_percent, "urban_landscape_percent.txt", row.names = FALSE, col.names = FALSE)
write.table(forest_landscape_percent, "forest_landscape_percent.txt", row.names = FALSE, col.names = FALSE)
write.table(agriculture_landscape_percent, "agriculture_landscape_percent.txt", row.names = FALSE, col.names = FALSE)

save.image("workspace.RData")