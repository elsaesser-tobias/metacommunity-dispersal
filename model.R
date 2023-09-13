
library(ggplot2);
library(deSolve);
library(simecol);
library(OCNet);
library(reticulate);


#clear enviroment
rm(list = ls())
setwd("C:/Users/elsae/Documents/R_Model")

#b will later be used to set different random seeds in each loop
b <- 13
#track time
start_time <- Sys.time()
#set random seed for reproducibility
set.seed(b)

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


#create the optimal channel network
OCN <- create_OCN(OCN_size, OCN_size, outletPos = 1, cellsize = cellsz)
OCN <- landscape_OCN(OCN)

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
OCN <- aggregate_OCN(landscape_OCN(OCN), thrA = threshold_area)



#problem: no stream order for each RN node
#get AG-stream order to RN where available
RNSO <- 0
for(i in 1:length(OCN$AG$toRN)){
  RNSO[OCN$AG$toRN[i]] <- OCN$AG$streamOrder[i]
}

#fill the rest through downstream node
#not well optimized but works

#the code takes all the places where the above node is available, but not the node itself
#the streamorder of the node itself is then filled with the one above
#this works because at AG level all changes in stream order are accounted for
for(j in 1:length(OCN$RN$downNode)){
  for(i in 2:length(OCN$RN$downNode)){
    if(is.na(RNSO[i]) == FALSE & is.na(RNSO[OCN$RN$downNode[i]])){
      RNSO[OCN$RN$downNode[i]] <- RNSO[i]
    }
  }
}


RNX <- OCN$RN$X
RNY <- OCN$RN$Y
RNDN <- OCN$RN$downNode


#add habitats
nhabitats <- round(length(RNX) * percentage_habitat_points)
set.seed(b)
habitat_points <- sample(1:length(RNX), nhabitats)

habitatsx <- RNX[habitat_points]
habitatsy <- RNY[habitat_points]



#export relevant data to files, so that they can be read by the python script
write.table(RNX, "x_points.txt", row.names = FALSE, col.names = FALSE)
write.table(RNY, "y_points.txt", row.names = FALSE, col.names = FALSE)
write.table(RNSO, "SO.txt", row.names = FALSE, col.names = FALSE)
write.table(RNDN, "DN.txt", row.names = FALSE, col.names = FALSE)
write.table(b, "py_seed.txt", row.names = FALSE, col.names = FALSE)
write.table(cellsz, "py_cellsz.txt", row.names = FALSE, col.names = FALSE)
write.table(habitatsx, "habitatsx.txt", row.names = FALSE, col.names = FALSE)
write.table(habitatsy, "habitatsy.txt", row.names = FALSE, col.names = FALSE)





#In the python script the landscape matrix gets create, the stream is drawn on top of it
#and the best dispersal pathways are determined and the cost is calculated.
#For details check the python script "create_matrix.py".
py_run_file("create_matrix.py")
image(py$nlmElement, useRaster=TRUE, axes=FALSE, col =gray(0:255/255))


draw_simple_OCN(OCN)
mycol <- rgb(0, 0, 255, max = 255, alpha = 125, names = "red50")
points(RNX[habitat_points], RNY[habitat_points], pch = 21, col = "blue", bg = mycol, cex = 5)
for(i in 1:length(py$startx)){
     x_plot <- c(py$startx[i],py$aimx[i])
     y_plot <- c(py$starty[i],py$aimy[i])
     lines(x_plot,y_plot)
   }





#step 2: preparations for the ODE model start here


#set up starting states:
#which / how many habitats should be occupied by what species at the start
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
    Pred1_spots[i] <- 1
  }
  else{
    Pred1_spots[i] <- 0.001
  }
}

for(i in 1:length(Pred2_spots)){
  if(Pred2_spots[i] <= percentage_Pred2){
    Pred2_spots[i] <- 1
  }
  else{
    Pred2_spots[i] <- 0.001
  }
}

for(i in 1:length(Prey1_spots)){
  if(Prey1_spots[i] <= percentage_Prey1){
    Prey1_spots[i] <- 1
  }
  else{
    Prey1_spots[i] <- 0.001
  }
}

for(i in 1:length(Prey2_spots)){
  if(Prey2_spots[i] <= percentage_Prey2){
    Prey2_spots[i] <- 1
  }
  else{
    Prey2_spots[i] <- 0.001
  }
}



#create names for all the habitat spots
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


#set the parameters for the model
parameters <- c(aPrey1 = 10.5,
                aPrey2 = 1.5,
                aPred1 = 5.75,
                aPred2 = 0.75,
                capacity_resource = 10,
                resource_growth_speed = 100,
                hPrey1 = 0.0001,
                hPrey2 = 0.01,
                hPred1 = 0.01,
                hPred2 = 0.1,
                mPrey1 = 0.52,
                mPrey2 = 0.12,
                mPred1 = 0.4,
                mPred2 = 0.2,
                modPred1 = 0.05,
                modPred2 = 0.05,
                modPrey1 = 0.05,
                modPrey2 = 0.05,
                n_loc = length(habitatsx))



migration_percent <- 0.4
n_habitat_connections <- c()
#create dispersal matrix D
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

#set up the ODEs as vector system
ODE_functions<-function(t, state, parameters) {
  #state[state<0] = 0
  with(as.list(parameters),{
  
  Pred1 <- state[1:n_loc]
  Pred2 <- state[(n_loc+1):(2*n_loc)]
  Prey1 <- state[(2*n_loc+1):(3*n_loc)]
  Prey2 <- state[(3*n_loc +1):(4*n_loc)]
  Res <- state[(4*n_loc+1):(5*n_loc)]
  
  
  dPred1 <- (((aPred1 * Pred1 * Prey1) / (1 + aPred1 * hPred1 * Prey1)) +
              ((aPred1 * Pred1 * Prey2) / (1 + aPred1 * hPred1 * Prey2))) * modPred1 -
              Pred1 * mPred1 + D %*% Pred1
  
  dPred2 <- (((aPred2 * Pred2 * Prey1) / (1 + aPred2 * hPred2 * Prey1)) +
              ((aPred2 * Pred2 * Prey2) / (1 + aPred2 * hPred2 * Prey2))) * modPred1 -
              Pred2 * mPred2  + D %*% Pred2
  
  dPrey1 <- ((aPrey1 * Prey1 * Res) / (1 + aPrey1 * hPrey1 * Res)) * modPrey1 -
    ((aPred1 * Pred1 * Prey1) / (1 + aPred1 * hPred1 * Prey1)) - 
    ((aPred2 * Pred2 * Prey1) / (1 + aPred2 * hPred2 * Prey1)) - 
    Prey1 * mPrey1 + D %*% Prey1
  
  dPrey2 <- ((aPrey2 * Prey2 * Res) / (1 + aPrey2 * hPrey2 * Res)) * modPrey2 -
    ((aPred1 * Pred1 * Prey2) / (1 + aPred1 * hPred1 * Prey2)) - 
    ((aPred2 * Pred2 * Prey2) / (1 + aPred2 * hPred2 * Prey2)) - 
    Prey2 * mPrey2 + D %*% Prey2
  
  dRes <- resource_growth_speed * (capacity_resource - Res) -
    ((aPrey1 * Prey1 * Res) / (1 + aPrey1 * hPrey1 * Res)) - 
    ((aPrey2 * Prey2 * Res) / (1 + aPrey2 * hPrey2 * Res))

  return(list(c(dPred1,dPred2,dPrey1,dPrey2,dRes)))
      })
  
}
  

initial_state <- c(Pred1_spots, Pred2_spots, Prey1_spots, Prey2_spots, res_spots)

times <- seq(0, 20, by = 0.1)

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
  for(j in 0:4){
    if(j == 0)
    result[[paste("habitat", i, sep = "")]][["Pred1"]] <- out[,j*20 + i + 1]
    if(j == 1)
      result[[paste("habitat", i, sep = "")]][["Pred2"]] <- out[,j*20 + i + 1]
    if(j == 2)
      result[[paste("habitat", i, sep = "")]][["Pray1"]] <- out[,j*20 + i + 1]
    if(j == 3)
      result[[paste("habitat", i, sep = "")]][["Pray2"]] <- out[,j*20 + i + 1]
    if(j == 4)
      result[[paste("habitat", i, sep = "")]][["Res"]] <- out[,j*20 + i + 1]
  }
}




#create a function to easily plot each habitat
plot_habitat <- function(x){
habitat_data <- as.data.frame(result[[paste("habitat", x, sep = "")]])

pred_color = "red"
pray_color = "blue"
res_color = "green"

ggplot(habitat_data, aes(x = time)) +
  geom_line(aes(y = Pred1, linetype = "Pred1", color = "Pred")) +
  geom_line(aes(y = Pred2, linetype = "Pred2", color = "Pred")) +
  geom_line(aes(y = Pray1, linetype = "Pray1", color = "Pray")) +
  geom_line(aes(y = Pray2, linetype = "Pray2", color = "Pray")) +
  geom_line(aes(y = Res, linetype = "Res", color = "Res")) +
  labs(title = paste("Habitat ", x, " Data", ", Number of connections to other habitats: ", 
                     n_habitat_connections[x], sep = ""),
       x = "Time",
       y = "Species Density",
       linetype = "Variable") +
  scale_linetype_manual(values = c("Pred1" = "solid", "Pred2" = "dashed", "Pray1" = "solid",
                                   "Pray2" = "dashed", "Res" = "solid")) +
  scale_color_manual(values = c(Pred = pred_color, Pray = pray_color, Res = res_color)) +
  theme_minimal()
}

plot_habitat(2)
#use plot_habitat(x) to plot habitat x



#check how long one iteration of the script takes
end_time <- Sys.time()
print(end_time - start_time)



