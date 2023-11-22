
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
setwd("C:/Users/elsae/Documents/R_Model")

#b will later be used to set different random seeds in each loop
b <- 14
#track time
start_time <- Sys.time()
#set random seed for reproducibility
set.seed(b)

#set landscape distribution
urban_landscape_percent <- 0.8
forest_landscape_percent <- 0.1
agriculture_landscape_percent <- 0.1

if(urban_landscape_percent + forest_landscape_percent + agriculture_landscape_percent != 1){
  stop("Not 100% of landscape assigned")
}

#in step one the OCN will be generated, habitats and migration routes will be calclulated
#the second step will then be the ODE model of the metacommunity

#size of the OCN will be OCN_size^2, cellsz determines the size of a single cell of the OCN matrix in meters
OCN_size <- 100
cellsz <- 100
#percentage of cells of the OCN that will later be counted as stream, the rest is considered drainage area
percentage_of_stream_points <- 0.1
number_of_stream_points <- (OCN_size * OCN_size) * percentage_of_stream_points

#percentage of the stream points that get used as as habitat
percentage_habitat_points <- 0.02


#create the optimal channel network -------------------------------------------------------
OCN <- create_OCN(OCN_size, OCN_size, outletPos = OCN_size / 2, cellsize = cellsz)
OCN <- landscape_OCN(OCN, zMin = 50000)

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
nhabitats <- round(length(RNX) * percentage_habitat_points)
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




#In the python script the landscape matrix gets create, the stream is drawn on top of it----------------------------
#and the best dispersal pathways are determined and the cost is calculated.
#For details check the python script "create_matrix.py".
py_run_file("create_matrix.py")



image(py$nlmElement, useRaster=TRUE, axes=FALSE, col =gray(0:255/255))


alt <- matrix(OCN$FD$Z, nrow = 100, ncol = 100)
alt <- alt / 250
expanded_alt <- kronecker(alt, matrix(1, nrow = 100, ncol = 100))
dim(expanded_alt)



# image(py$nlmElement, useRaster=TRUE, axes=FALSE, col =gray(0:255/255))
# image(expanded_alt, useRaster=TRUE, axes=FALSE, col =gray(0:255/255))
# 
# write.csv(expanded_alt, file = "alt.csv")
# write.csv(py$nlmElement, file = "landscape.csv")


#Zonenhistogramm




#x coordinate, y coordinate, radius
count_landscape <- function(mat, x, y, r, landscape_type){
  x <- x + 3000
  y <- y + 3000
  
  
  # Create a matrix of distances from the center
  distances <- sqrt((row(mat) - x)^2 + (col(mat) - y)^2)
  
  # Select elements within the circular region
  circle_elements <- mat[distances <= r]
  
  # Print the elements within the circular region
  #print(circle_elements)
  if(landscape_type == "stream"){return(sum(circle_elements == 0.1))}
  
  if(landscape_type == "agriculture"){return(sum(circle_elements == 0.5))}
  
  if(landscape_type == "forrest"){return(sum(circle_elements == 0.75))}
  
  if(landscape_type == "urban"){return(sum(circle_elements == 1))}
  
}

start_time <- Sys.time()

urban_radius <- 3000
agriculture_radius <- 750
forrest_radius <- 1000

urban_sums <- c()
agriculture_sums <- c()
forrest_sums <- c()
altitudes <- c()

for(i in 1:length(habitatsx)){
  urban_sums[i] <- count_landscape(py$nlmElement_orig, habitatsx[i],habitatsy[i], urban_radius, "urban")
  agriculture_sums[i] <- count_landscape(py$nlmElement_orig, habitatsx[i],habitatsy[i], agriculture_radius, "agriculture")
  forrest_sums[i] <- count_landscape(py$nlmElement_orig, habitatsx[i],habitatsy[i], forrest_radius, "forrest")
  altitudes[i] <- expanded_alt[habitatsx[i],habitatsy[i]]
}


end_time <- Sys.time()
print(end_time - start_time)


habitat_df <- data.frame(urban3000m.m2 = urban_sums, agriculture750m.m2 = agriculture_sums, forest1000m.m2 = forrest_sums, 
                         alt = altitudes)

#                     Rhine Valley Pesticides 2022                          #
#                   Ken Mauser, PhD Project - SystemLink                    #
#                                                                           #
#                 Build prediction maps Tobias                              #
#                                                                           #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#


# preparation -----------------------------------------------------------------


# import data ------------------------------------------------------------------


histo_cat <- read_excel("histo_cat.xlsx")
veg <- read_excel("veg.xls")

# change format ----------------------------------------------------------------


histo_cat$agriculture750m.m2 <- histo_cat$agriculture750m*100
histo_cat$forest1000m.m2 <- histo_cat$forest1000m*100
histo_cat$urban3000m.m2 <- histo_cat$forest3000m*100

## join histogram data with measurement points ---------------------------------

icol <- c("agriculture750m.m2", "forest1000m.m2", "urban3000m.m2", "sample_number")

veg <- left_join(x = veg, y = histo_cat[,icol], by = "sample_number")

# prediction Tobias ------------------------------------------------------------


set.seed(1)

pcr1 <- pcr(num ~ alt + 
              forest1000m.m2 +
              urban3000m.m2 +
              agriculture750m.m2,
            data = veg[veg$borderdistance > 3000,], scale = T, validation = "CV")


# pcr1 <- lm(num ~ alt + 
#               forest1000m.m2 +
#               urban3000m.m2 +
#               agriculture750m.m2,
#             data = veg[veg$borderdistance > 3000,])


summary(pcr1)
plot(pcr1)
validationplot(object = pcr1, val.type = "R2")

max(R2(pcr1)[[1]])

pred1 <- predict(pcr1, newdata = habitat_df, ncomp = 1)


summary(pred1)

predhabitat <- habitat_df

predhabitat$prediction <- as.vector(pred1)

habitat_quality <- 100 / predhabitat$prediction
#----------------------------------------------------------------------------------



# draw_simple_OCN(OCN)
# mycol <- rgb(0, 0, 255, max = 255, alpha = 125, names = "red50")
# points(RNX[habitat_points], RNY[habitat_points], pch = 21, col = "blue", bg = mycol, cex = 5)
# for(i in 1:length(py$startx)){
#      x_plot <- c(py$startx[i],py$aimx[i])
#      y_plot <- c(py$starty[i],py$aimy[i])
#      lines(x_plot,y_plot)
#    }



#step 2: preparations for the ODE model start here
#set up starting states:
#which / how many habitats should be occupied by what species at the start -------------------------------------
percentage_Pred1 <- 0.5
percentage_Pred2 <- 0.5
percentage_Prey1 <- 0.5
percentage_Prey2 <- 0.5
res_spots <- rep(1, times = length(habitatsx))

set.seed(b)
Pred1_spots <- runif(length(habitatsx))

set.seed(b + 100000000)
Pred2_spots <- runif(length(habitatsx))

set.seed(b + 200000000)
Prey1_spots <- runif(length(habitatsx))

set.seed(b + 300000000)
Prey2_spots <- runif(length(habitatsx))

for(i in 1:length(Pred1_spots)){
  if(Pred1_spots[i] <= percentage_Pred1){
    Pred1_spots[i] <- 20
  }
  else{
    Pred1_spots[i] <- 0
  }
}

for(i in 1:length(Pred2_spots)){
  if(Pred2_spots[i] <= percentage_Pred2){
    Pred2_spots[i] <- 20
  }
  else{
    Pred2_spots[i] <- 0
  }
}

for(i in 1:length(Prey1_spots)){
  if(Prey1_spots[i] <= percentage_Prey1){
    Prey1_spots[i] <- 20
  }
  else{
    Prey1_spots[i] <- 0
  }
}

for(i in 1:length(Prey2_spots)){
  if(Prey2_spots[i] <= percentage_Prey2){
    Prey2_spots[i] <- 20
  }
  else{
    Prey2_spots[i] <- 0
  }
}



#create names for all the habitat spots ---------------------------------------------------
habitat_names <- c()
for(i in 1:length(habitatsx)){
  habitat_names[i] <- paste("X", habitatsx[i], "Y", habitatsy[i], sep = "")
}

start_names <- c()
aim_names <- c()
for(i in 1:length(py$startx)){
  start_names[i] <- paste("X", py$startx[i], "Y", py$starty[i], sep = "")
  aim_names[i] <- paste("X", py$aimx[i], "Y", py$aimy[i], sep = "")
}


#set the parameters for the model -----------------------------------------------------
parameters <- c(aPred1 = 0.00444,
                aPred2 = 0.00444,
                hPred1 = 0.0192,
                hPred2 = 0.0192,
                mPred1 = 0.11,
                mPred2 = 0.1,
                modPred1 = 1,
                modPred2 = 1,
                prey1_growth_speed = 0.005,
                prey2_growth_speed = 0.005,
                capacity_prey1 = 50,
                capacity_prey2 = 50,
                n_loc = length(habitatsx))





#create dispersal matrix D----------------------------------------------------------------------------

migration_percent <- 0.4
n_habitat_connections <- c()
D <- matrix(0, nrow = length(habitatsx),ncol = length(habitatsx))
for(i in 1:length(aim_names)){
    
    start_pos <- which(grepl(start_names[i], habitat_names))
    aim_pos <- which(grepl(aim_names[i], habitat_names))
    
    D[start_pos,aim_pos] <- (migration_percent / sum(grepl(start_names[i], start_names))) *
        ((1250 - py$route_weights[[i]]) / 1250)
    n_habitat_connections[start_pos] <- sum(grepl(start_names[i], start_names))
  
}
diag(D) <- - rowSums(D)
n_habitat_connections[is.na(n_habitat_connections)] <- 0

#set up the ODEs as vector system-----------------------------------------------------------------------

ODE_functions<-function(t, state, parameters) {
  #state[state<0] = 0
  with(as.list(parameters),{
  
  Pred1 <- state[1:n_loc]
  Pred2 <- state[(n_loc+1):(2*n_loc)]
  Prey1 <- state[(2*n_loc+1):(3*n_loc)]
  Prey2 <- state[(3*n_loc +1):(4*n_loc)]

  
  
  dPred1 <- (((aPred1 * Pred1 * Prey1) / (1 + aPred1 * hPred1 * Prey1)) +
              ((aPred1 * Pred1 * Prey2) / (1 + aPred1 * hPred1 * Prey2))) * modPred1 -
              Pred1 * mPred1 + D %*% Pred1
  
  dPred2 <- (((aPred2 * Pred2 * Prey1) / (1 + aPred2 * hPred2 * Prey1)) +
              ((aPred2 * Pred2 * Prey2) / (1 + aPred2 * hPred2 * Prey2))) * modPred1 -
              Pred2 * mPred2  + D %*% Pred2
  
  
  dPrey1 <- prey1_growth_speed * (habitat_quality * capacity_prey1 - Prey1) -
    ((aPred1 * Pred1 * Prey1) / (1 + aPred1 * hPred1 * Prey1)) - 
    ((aPred2 * Pred2 * Prey1) / (1 + aPred2 * hPred2 * Prey1))
  
  dPrey2 <- prey2_growth_speed * (habitat_quality * capacity_prey1 - Prey2) -
    ((aPred1 * Pred1 * Prey2) / (1 + aPred1 * hPred1 * Prey2)) - 
    ((aPred2 * Pred2 * Prey2) / (1 + aPred2 * hPred2 * Prey2))

  return(list(c(dPred1,dPred2,dPrey1,dPrey2)))
      })
  
}
  

initial_state <- c(Pred1_spots, Pred2_spots, Prey1_spots, Prey2_spots)

times <- seq(0, 100, by = 0.1)

out <- lsoda(y = initial_state, times = times, func = ODE_functions, parms = parameters, 
             rtol = 1e-12, atol = 1e-12, maxsteps = 10000)
#rtol = 1e-12, atol = 1e-12, maxsteps = 1000000
#head(out)
#plot(out)


#check if there are negative numbers
if(any(out < 0)){
  print("There are negative numbers.")
}



#in case there are negative numbers this can be used to analyse where they are
# negative_indices <- which(out < 0, arr.ind = TRUE)
# 
# if (length(negative_indices) > 0) {
#   negative_values <- out[negative_indices]
#   print(paste("Negative elements found at indices:", toString(negative_indices)))
#   print(paste("Negative values:", toString(negative_values)))
# } else {
#   print("No negative elements found in the matrix.")
# }


#make a nice visualization plot
result <- list()
#sort results by habitat
for(i in 1:length(habitatsx)){
  result[[paste("habitat", i, sep = "")]][["time"]] <- out[,1]
  for(j in 0:3){
    if(j == 0)
    result[[paste("habitat", i, sep = "")]][["Pred1"]] <- out[,j*20 + i + 1]
    if(j == 1)
      result[[paste("habitat", i, sep = "")]][["Pred2"]] <- out[,j*20 + i + 1]
    if(j == 2)
      result[[paste("habitat", i, sep = "")]][["Prey1"]] <- out[,j*20 + i + 1]
    if(j == 3)
      result[[paste("habitat", i, sep = "")]][["Prey2"]] <- out[,j*20 + i + 1]
  }
}




#create a function to easily plot each habitat----------------------------------------
plot_habitat <- function(x){
habitat_data <- as.data.frame(result[[paste("habitat", x, sep = "")]])

pred_color = "red"
pray_color = "blue"


ggplot(habitat_data, aes(x = time)) +
  geom_line(aes(y = Pred1, linetype = "Pred1", color = "Pred")) +
  geom_line(aes(y = Pred2, linetype = "Pred2", color = "Pred")) +
  geom_line(aes(y = Prey1, linetype = "Prey1", color = "Prey")) +
  geom_line(aes(y = Prey2, linetype = "Prey2", color = "Prey")) +
  labs(title = paste("Habitat ", x, " Data", ", Number of connections to other habitats: ", 
                     n_habitat_connections[x], sep = ""),
       x = "Time",
       y = "Species Density",
       linetype = "Variable") +
  scale_linetype_manual(values = c("Pred1" = "solid", "Pred2" = "dashed", "Prey1" = "solid",
                                   "Prey2" = "dashed")) +
  scale_color_manual(values = c(Pred = pred_color, Pray = pray_color)) +
  theme_minimal()
}

plot_habitat(4)
#use plot_habitat(x) to plot habitat x



#check how long one iteration of the script takes
end_time <- Sys.time()
print(end_time - start_time)

#checken ob Prozentuale Verteilung richtig ist
<<<<<<< HEAD
#sum(py$nlmElement_orig == 0.1)
#dim(py$nlmElement_orig)

#image(py$nlmElement, useRaster=TRUE, axes=FALSE, col =gray(0:255/255))

sum_pred1 <- 0
sum_pred2 <- 0
sum_prey1 <- 0
sum_prey2 <- 0

l <- 1
for(i in result){
  sum_pred1 <- i$Pred1[length(times)] + sum_pred1
  sum_pred2 <- i$Pred2[length(times)] + sum_pred2
  sum_prey1 <- i$Prey1[length(times)] + sum_prey1
  sum_prey2 <- i$Prey2[length(times)] + sum_prey2
  l <- l +1
}

barplot(c(sum_pred1,sum_pred2, sum_prey1, sum_prey2))

=======
# comment
>>>>>>> 8d585cee4a26c116f5d38490a1d09dfd62685689
