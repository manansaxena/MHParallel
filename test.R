library(MHParallel)

# pred_eta <- as.array(fit.mc_d[[1]][["Eta"]])
# #true_eta <- as.array(fit.mc_d[[1]][["mdataset"]][["Lambda_true"]])
# 
# B <- as.array(fit.mc_d[[1]][["mdataset"]][["Theta"]] %*% fit.mc_d[[1]][["mdataset"]][["X"]])
# K <- as.array(fit.mc_d[[1]][["mdataset"]][["Xi"]])
# A <- as.array(diag(fit.mc_d[[1]][["N"]]) + t(fit.mc_d[[1]][["mdataset"]][["X"]]) %*% fit.mc_d[[1]][["mdataset"]][["Gamma"]] %*% fit.mc_d[[1]][["mdataset"]][["X"]])
# Y <- as.array(fit.mc_d[[1]][["mdataset"]][["Y"]])
# mean_eta <- matrix(0,fit.mc_d[[1]][["D"]]-1,fit.mc_d[[1]][["N"]])
# 
# for(i in c(1:fit.mc_d[[1]][["N"]])){
#   mean_eta[,i] <- as.matrix(rowMeans(fit.mc_d[[1]][["Eta"]][,i,]))
# }
# #stepped_eta_1 <- array(0, dim=c(1000,dim(pred_eta)[1],dim(pred_eta)[2]))
# 
# #stepped_eta_1[1,,] <- pred_eta[,,1]
# 
# #for(i in c(2:1000)){
# #  stepped_eta_1[i,,] <- MHParallel::mh_step_interface(stepped_eta_1[i-1,,],B,K,A,"normal",i)
# #}
# 
# result <- MHParallel::mh_interface(pred_eta[,,1],B,K,A,Y,mean_eta,2000,1)

# one chain long time
num_chains <- 1
initial_state <- as.array(c(1:num_chains))
result <- MHParallel::metro_interface(initial_state,num_chains,10000,0,4,0,as.array(c(0,1)))

(result[["Timer"]][["Overalln Stop"]]-result[["Timer"]][["Overall Start"]])/10^6

data <- result[["particles"]][[1]]

hist(data, main = "Histogram with Overlayed Density Curve", prob = TRUE, col = "lightblue")
 
# Fit a density curve
density_data <- density(data)

# Overlay the density curve
lines(density_data, col = "darkgreen", lwd = 2)

#multiple chains decent time
num_chains <- 10
initial_state <- as.array(c(1:num_chains))
result <- MHParallel::metro_interface(initial_state,num_chains,2000,0,4,0,as.array(c(0,1)))


data <- result[["particles"]][[2]]

hist(data, main = "Histogram with Overlayed Density Curve", prob = TRUE, col = "lightblue")

# Fit a density curve
density_data <- density(data)

# Overlay the density curve
lines(density_data, col = "darkgreen", lwd = 2)

(result[["Timer"]][["Overalln Stop"]]-result[["Timer"]][["Overall Start"]])/10^6




