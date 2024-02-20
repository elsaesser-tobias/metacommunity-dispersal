
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

#read data creeated in previous R and python files
load("workspace.RData")

nlmElement <- read.table("nlmElement.csv", header = FALSE, sep =",")
nlmElement <- as.numeric(unlist(nlmElement[sapply(nlmElement, is.numeric)]))
nlmElement <- matrix(nlmElement, nrow = 10000, ncol = 10000)


#---------------------------------------------------------------------------------------------------------------#
#create dispersal paths and dispersal matrix D


create_paths <- function(start_point_x, start_point_y, n_paths, distance, hab_x, hab_y, rad){
  n_match <- rep(0, length(habitatsx))
  for(i in 1:n_paths){
    
    
    route_x <- c(start_point_x)
    route_y <- c(start_point_y)
    route_paste <- c()
    route_paste[1] <- paste("x",start_point_x,"y",start_point_y, sep = "")
    route_dir <- c()
    route_which <- 2
    
    current_x <- start_point_x
    current_y <- start_point_y
    dir <- 0
    for(j in 1:distance){
      #check if close to edge, if yes, abandon path
      if(current_x == 1 | current_x == nrow(nlmElement) | 
         current_y == 1 | current_y == nrow(nlmElement)){
        break
      }
      
      
      #get surrounding
      around <- c()
    
      around [1] <- nlmElement[current_x - 1,current_y - 1] 
      around [2] <- nlmElement[current_x,current_y - 1]
      around [3] <- nlmElement[current_x + 1,current_y - 1]
      around [4] <- nlmElement[current_x + 1,current_y]
      around [5] <- nlmElement[current_x + 1,current_y + 1]
      around [6] <- nlmElement[current_x,current_y + 1]
      around [7] <- nlmElement[current_x - 1,current_y + 1]
      around [8] <- nlmElement[current_x - 1,current_y]
      possible_dir <- c()
      if(dir == 1){
        possible_dir <- c(7,8,1,2,3)
      } else if(dir == 2){
        possible_dir <- c(8,1,2,3,4)
      } else if(dir == 3){
        possible_dir <- c(1,2,3,4,5)
      } else if(dir == 4){
        possible_dir <- c(2,3,4,5,6)
      } else if(dir == 5){
        possible_dir <- c(3,4,5,6,7)
      } else if(dir == 5){
        possible_dir <- c(3,4,5,6,7)
      } else if(dir == 6){
        possible_dir <- c(4,5,6,7,8)
      } else if(dir == 7){
        possible_dir <- c(5,6,7,8,1)
      } else if(dir == 8){
        possible_dir <- c(6,7,8,1,2)
      } else {
        possible_dir <- c(1,2,3,4,5,6,7,8)
      }
      
      sum(around[possible_dir])
      #set multiplier for direction
      multiplier <- c()
      multiplier <- rep(0.01 , length(possible_dir))
      
      multiplier[which(possible_dir == dir)] <- 0.9
      multiplier[which(possible_dir == dir) + 1] <- 0.04
      multiplier[which(possible_dir == dir) - 1] <- 0.04
      
      around[possible_dir] <- around[possible_dir] * multiplier
      
      
      set.seed((j + i + b))
      dir <- sample(possible_dir, size = 1, prob = around[possible_dir])
      
      if(dir == 1){
        current_x <- current_x - 1
        current_y <- current_y -1
      } else if(dir == 2){
        current_x <- current_x
        current_y <- current_y -1
      } else if(dir == 3){
        current_x <- current_x + 1
        current_y <- current_y -1
      } else if(dir == 4){
        current_x <- current_x + 1
        current_y <- current_y
      } else if(dir == 5){
        current_x <- current_x + 1
        current_y <- current_y + 1
      } else if(dir == 6){
        current_x <- current_x
        current_y <- current_y + 1
      } else if(dir == 7){
        current_x <- current_x - 1
        current_y <- current_y + 1
      } else if(dir == 8){
        current_x <- current_x - 1
        current_y <- current_y
      }
      
      route_x[route_which] <- current_x
      route_y[route_which] <- current_y
      route_dir[route_which - 1] <- dir
      route_which <- route_which + 1  
      
    }
    matrix_x <- matrix(route_x, nrow = length(route_x), ncol = length(hab_x), byrow = FALSE)
    matrix_y <- matrix(route_y, nrow = length(route_y), ncol = length(hab_y), byrow = FALSE)
    matrix_habx <- matrix(hab_x, nrow = length(route_x), ncol = length(hab_x), byrow = TRUE)
    matrix_haby <- matrix(hab_y, nrow = length(route_y), ncol = length(hab_y), byrow = TRUE)
    
    arrived <- rep(NA,length(hab_x))
    dist <- sqrt(((matrix_x - matrix_habx)^2) + ((matrix_y - matrix_haby)^2)) <rad
    for(q in 1:length(hab_x)){
      if(any(dist[,q])){
        arrived[q] <- which.max(dist[,q])
      }
    }
    
    n_match[which.min(arrived)] <- n_match[which.min(arrived)] + 1

  }
  return(n_match)
}

D <- matrix(NA, nrow = length(habitatsx), ncol = length(habitatsx))


for(i in 1: length(habitatsx)){
  matches <- c()
  new_habitatsx <- habitatsx
  new_habitatsx[i] <- 100000
  new_habitatsy <- habitatsy
  new_habitatsy[i] <- 100000
  matches <- create_paths(habitatsx[i],habitatsy[i], n_dragonflies, n_steps,new_habitatsx, new_habitatsy, hab_radius)
  D[i,] <- matches
  

}
D <- D *0.001 * dispersal_percent
diag(D) <- - rowSums(D)

n_habitat_connections <- c()
for(i in 1:length(D[1,])){
  n_habitat_connections[i] <- sum(D[i,] > 0)
}
#--------------------------------------------------------------------------------------------------------------#
#create catchment areas from matrix
threshold_area2 <- 0
OCN2 <- aggregate_OCN(landscape_OCN(OCN2, zMin = 50000), thrA = threshold_area2)

# draw_simple_OCN(OCN = OCN2, thrADraw = threshold_area2)
# draw_simple_OCN(OCN = OCN, thrADraw = threshold_area)


FD_habitat_points <- OCN$RN$toFD[habitat_points]
catchment_areas <- OCN2$RN$upstream[FD_habitat_points]

OCN$RN$X[habitat_points]
OCN2$RN$X[FD_habitat_points]

stream_amount <- numeric(length(habitat_points))
agri_amount <- numeric(length(habitat_points))
total_amount <- numeric(length(habitat_points))
a = 1
cell <- cellsz /2
for(i in catchment_areas){
  for(j in i){
    x <- OCN2$RN$X[j]
    y <- OCN2$RN$X[j]
    stream_amount[a] <- sum(nlmElement[(x-cell+1):(x+cell), (y-cell+1):(y+cell)] == 0.1) + stream_amount[a]
    agri_amount[a] <- sum(nlmElement[(x-cell+1):(x+cell), (y-cell+1):(y+cell)] == 0.5) + agri_amount[a]
    total_amount[a] <- sum(nlmElement[(x-cell+1):(x+cell), (y-cell+1):(y+cell)] > 0) + total_amount[a]
  }
  a = a + 1
}


habitat_quality <- 1 - (agri_amount / (total_amount - stream_amount)) * 0.75





# alt <- matrix(OCN$FD$Z, nrow = 100, ncol = 100)
# alt <- alt / 250
# expanded_alt <- kronecker(alt, matrix(1, nrow = 100, ncol = 100))
# dim(expanded_alt)
# 
# 
# 
# # image(py$nlmElement, useRaster=TRUE, axes=FALSE, col =gray(0:255/255))
# # image(expanded_alt, useRaster=TRUE, axes=FALSE, col =gray(0:255/255))
# # 
# # write.csv(expanded_alt, file = "alt.csv")
# # write.csv(py$nlmElement, file = "landscape.csv")
# 
# 
# #Zonenhistogramm
# 
# 
# 
# 
# #x coordinate, y coordinate, radius
# count_landscape <- function(mat, x, y, r, landscape_type){
#   x <- x + 3000
#   y <- y + 3000
#   
#   
#   # Create a matrix of distances from the center
#   distances <- sqrt((row(mat) - x)^2 + (col(mat) - y)^2)
#   
#   # Select elements within the circular region
#   circle_elements <- mat[distances <= r]
#   
#   # Print the elements within the circular region
#   #print(circle_elements)
#   if(landscape_type == "stream"){return(sum(circle_elements == 0.1))}
#   
#   if(landscape_type == "agriculture"){return(sum(circle_elements == 0.5))}
#   
#   if(landscape_type == "forrest"){return(sum(circle_elements == 0.75))}
#   
#   if(landscape_type == "urban"){return(sum(circle_elements == 1))}
#   
# }
# 
# 
# 
# 
# 
# urban_radius <- 3000
# agriculture_radius <- 750
# forrest_radius <- 1000
# 
# urban_sums <- c()
# agriculture_sums <- c()
# forrest_sums <- c()
# altitudes <- c()
# 
# for(i in 1:length(habitatsx)){
#   urban_sums[i] <- count_landscape(nlmElement_orig, habitatsx[i],habitatsy[i], urban_radius, "urban")
#   agriculture_sums[i] <- count_landscape(nlmElement_orig, habitatsx[i],habitatsy[i], agriculture_radius, "agriculture")
#   forrest_sums[i] <- count_landscape(nlmElement_orig, habitatsx[i],habitatsy[i], forrest_radius, "forrest")
#   altitudes[i] <- expanded_alt[habitatsx[i],habitatsy[i]]
# }
# 
# 
# 
# 
# habitat_df <- data.frame(urban3000m.m2 = urban_sums, agriculture750m.m2 = agriculture_sums, forest1000m.m2 = forrest_sums, 
#                          alt = altitudes)
# 
# #                     Rhine Valley Pesticides 2022                          #
# #                   Ken Mauser, PhD Project - SystemLink                    #
# #                                                                           #
# #                 Build prediction maps for Tobias                          #
# #                                                                           #
# #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# 
# 
# # preparation -----------------------------------------------------------------
# 
# 
# # import data ------------------------------------------------------------------
# 
# 
# histo_cat <- read_excel("histo_cat.xlsx")
# veg <- read_excel("veg.xls")
# 
# # change format ----------------------------------------------------------------
# 
# 
# histo_cat$agriculture750m.m2 <- histo_cat$agriculture750m*100
# histo_cat$forest1000m.m2 <- histo_cat$forest1000m*100
# histo_cat$urban3000m.m2 <- histo_cat$forest3000m*100
# 
# ## join histogram data with measurement points ---------------------------------
# 
# icol <- c("agriculture750m.m2", "forest1000m.m2", "urban3000m.m2", "sample_number")
# 
# veg <- left_join(x = veg, y = histo_cat[,icol], by = "sample_number")
# 
# # prediction Tobias ------------------------------------------------------------
# 
# 
# set.seed(b)
# 
# pcr1 <- pcr(num ~ alt + 
#               forest1000m.m2 +
#               urban3000m.m2 +
#               agriculture750m.m2,
#             data = veg[veg$borderdistance > 3000,], scale = T, validation = "CV")
# 
# 
# # pcr1 <- lm(num ~ alt + 
# #               forest1000m.m2 +
# #               urban3000m.m2 +
# #               agriculture750m.m2,
# #             data = veg[veg$borderdistance > 3000,])
# 
# 
# summary(pcr1)
# plot(pcr1)
# validationplot(object = pcr1, val.type = "R2")
# 
# max(R2(pcr1)[[1]])
# 
# pred1 <- predict(pcr1, newdata = habitat_df, ncomp = 1)
# 
# 
# summary(pred1)
# 
# predhabitat <- habitat_df
# 
# predhabitat$prediction <- as.vector(pred1)
# 
# habitat_quality <- 100 / predhabitat$prediction
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
percentage_Pred1 <- 1
percentage_Pred2 <- 1
percentage_Prey1 <- 1
percentage_Prey2 <- 1
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
# habitat_names <- c()
# for(i in 1:length(habitatsx)){
#   habitat_names[i] <- paste("X", habitatsx[i], "Y", habitatsy[i], sep = "")
# }
# 
# start_names <- c()
# aim_names <- c()
# for(i in 1:length(startx)){
#   start_names[i] <- paste("X", startx[i], "Y", starty[i], sep = "")
#   aim_names[i] <- paste("X", aimx[i], "Y", aimy[i], sep = "")
# }


#set the parameters for the model -----------------------------------------------------
parameters <- c(aPred1 = 0.00444,
                aPred2 = 0.00444,
                hPred1 = 0.0192,
                hPred2 = 0.0192,
                mPred1 = 0.11,
                mPred2 = 0.11,
                modPred1 = 0.05,
                modPred2 = 0.05,
                prey1_growth_speed = 0.05,
                prey2_growth_speed = 0.05,
                capacity_prey1 = 4000,
                capacity_prey2 = 4000,
                n_loc = length(habitatsx))





#create dispersal matrix D----------------------------------------------------------------------------
# 
# migration_percent <- 0.5
# n_habitat_connections <- c()
# D <- matrix(0, nrow = length(habitatsx),ncol = length(habitatsx))
# for(i in 1:length(aim_names)){
# 
#     start_pos <- which(grepl(start_names[i], habitat_names))
#     aim_pos <- which(grepl(aim_names[i], habitat_names))
# 
#     D[start_pos,aim_pos] <- (migration_percent / sum(grepl(start_names[i], start_names))) *
#         ((1250 - route_weights[[i]]) / 1250)
#     n_habitat_connections[start_pos] <- sum(grepl(start_names[i], start_names))
# 
# }
# diag(D) <- - rowSums(D)
# n_habitat_connections[is.na(n_habitat_connections)] <- 0

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
              ((aPred2 * Pred2 * Prey2) / (1 + aPred2 * hPred2 * Prey2))) * modPred2 -
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

plot_habitat(3)
#use plot_habitat(x) to plot habitat x



#check how long one iteration of the script takes


#checken ob Prozentuale Verteilung richtig ist
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




#points(5000, 5000, pch = 21, col = "blue", cex = 515)


png(file=paste("barplot",b,".png", sep =""), width=1000, height=1000, bg = "white")


barplot(c(sum_pred1,sum_pred2, sum_prey1, sum_prey2), names.arg = c("Pred1", "Pred2", "Prey1", "Prey2"))


dev.off() 



end_time <- Sys.time()
print(end_time - start_time)

save.image(paste("workspace_end",b, ".RData", sep = ""))































# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# create_paths2 <- function(start_point_x, start_point_y, n_paths, distance, hab_x, hab_y){
#   n_match <- rep(0, length(habitatsx))
#   for(i in 1:n_paths){
# 
# 
#     route_x <- c(start_point_x)
#     route_y <- c(start_point_y)
#     route_paste <- c()
#     route_paste[1] <- paste("x",start_point_x,"y",start_point_y, sep = "")
#     route_dir <- c()
#     route_which <- 2
# 
#     current_x <- start_point_x
#     current_y <- start_point_y
#     dir <- 0
#     for(j in 1:distance){
#       #check if close to edge, if yes, abandon path
#       if(current_x == 1 | current_x == nrow(nlmElement) |
#          current_y == 1 | current_y == nrow(nlmElement)){
#         break
#       }
# 
# 
#       #get surrounding
#       around <- c()
# 
#       around [1] <- nlmElement[current_x - 1,current_y - 1]
#       around [2] <- nlmElement[current_x,current_y - 1]
#       around [3] <- nlmElement[current_x + 1,current_y - 1]
#       around [4] <- nlmElement[current_x + 1,current_y]
#       around [5] <- nlmElement[current_x + 1,current_y + 1]
#       around [6] <- nlmElement[current_x,current_y + 1]
#       around [7] <- nlmElement[current_x - 1,current_y + 1]
#       around [8] <- nlmElement[current_x - 1,current_y]
#       possible_dir <- c()
#       if(dir == 1){
#         possible_dir <- c(7,8,1,2,3)
#       } else if(dir == 2){
#         possible_dir <- c(8,1,2,3,4)
#       } else if(dir == 3){
#         possible_dir <- c(1,2,3,4,5)
#       } else if(dir == 4){
#         possible_dir <- c(2,3,4,5,6)
#       } else if(dir == 5){
#         possible_dir <- c(3,4,5,6,7)
#       } else if(dir == 5){
#         possible_dir <- c(3,4,5,6,7)
#       } else if(dir == 6){
#         possible_dir <- c(4,5,6,7,8)
#       } else if(dir == 7){
#         possible_dir <- c(5,6,7,8,1)
#       } else if(dir == 8){
#         possible_dir <- c(6,7,8,1,2)
#       } else {
#         possible_dir <- c(1,2,3,4,5,6,7,8)
#       }
# 
#       sum(around[possible_dir])
#       #set multiplier for direction
#       multiplier <- c()
#       multiplier <- rep(0.01 , length(possible_dir))
# 
#       multiplier[which(possible_dir == dir)] <- 0.9
#       multiplier[which(possible_dir == dir) + 1] <- 0.04
#       multiplier[which(possible_dir == dir) - 1] <- 0.04
# 
#       around[possible_dir] <- around[possible_dir] * multiplier
# 
# 
#       set.seed((j + i +b))
#       dir <- sample(possible_dir, size = 1, prob = around[possible_dir])
# 
#       if(dir == 1){
#         current_x <- current_x - 1
#         current_y <- current_y -1
#       } else if(dir == 2){
#         current_x <- current_x
#         current_y <- current_y -1
#       } else if(dir == 3){
#         current_x <- current_x + 1
#         current_y <- current_y -1
#       } else if(dir == 4){
#         current_x <- current_x + 1
#         current_y <- current_y
#       } else if(dir == 5){
#         current_x <- current_x + 1
#         current_y <- current_y + 1
#       } else if(dir == 6){
#         current_x <- current_x
#         current_y <- current_y + 1
#       } else if(dir == 7){
#         current_x <- current_x - 1
#         current_y <- current_y + 1
#       } else if(dir == 8){
#         current_x <- current_x - 1
#         current_y <- current_y
#       }
# 
#       route_x[route_which] <- current_x
#       route_y[route_which] <- current_y
#       route_dir[route_which - 1] <- dir
#       route_which <- route_which + 1
# 
#     }
#     matrix_x <- matrix(route_x, nrow = length(route_x), ncol = length(hab_x), byrow = FALSE)
#     matrix_y <- matrix(route_y, nrow = length(route_y), ncol = length(hab_y), byrow = FALSE)
#     matrix_habx <- matrix(hab_x, nrow = length(route_x), ncol = length(hab_x), byrow = TRUE)
#     matrix_haby <- matrix(hab_y, nrow = length(route_y), ncol = length(hab_y), byrow = TRUE)
# 
#     arrived <- rep(NA,length(hab_x))
#     dist <- sqrt(((matrix_x - matrix_habx)^2) + ((matrix_y - matrix_haby)^2)) <500
#     for(q in 1:length(hab_x)){
#       if(any(dist[,q])){
#         arrived[q] <- which.max(dist[,q])
#       }
#     }
# 
#     n_match[which.min(arrived)] <- n_match[which.min(arrived)] + 1
# 
#   }
#   return(data.frame(route_x,route_y))
# }
# 
# new_habitatsx <- habitatsx
# new_habitatsx[15] <- 100000
# new_habitatsy <- habitatsy
# new_habitatsy[15] <- 100000
# #---------------------------------------------------------------------#
# 
# A <- create_paths2(habitatsx[15],habitatsy[15], 1, 10000,new_habitatsx, new_habitatsy)
# p <- sum(nlmElement[cbind(A[[1]],A[[2]])] == 0.1) / length(A[[1]])
# 
# min_x <- min(A[[1]])
# max_x <- max(A[[1]])
# 
# max_y <- max(A[[2]])
# min_y <- min(A[[2]])
# 
# colors <- c("red", "green", "yellow", "blue")
# 
# # Define breaks for color scaling
# breaks <- c(0, 0.00001, 0.001, 0.01, 1)
# image(min_x:max_x, min_y:max_y, nlmElement[min_x:max_x,min_y:max_y], breaks = breaks, col = colors)
# lines(A)
# A[1,2]
# points(A[1,1],A[1,2], lwd = 5, col = "purple")
# 
