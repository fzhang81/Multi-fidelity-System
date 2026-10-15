#setwd("~/Downloads")
library(laGP)
library(plgp)
library(lhs)
library(Matrix)
library(mvtnorm) 

# Iterate over both "M" and "F" types
for (type in c("M", "F")) {
  for (i in 1:50) {
    
    # Construct filename dynamically
    filename <- paste0("2025618", type, "_p5_noise02_", i)
    readpath <- paste0("Calidata/", filename, ".csv")
    
    ball2 <- read.csv(readpath)
    ball2_0 <- ball2
    
    # Standardize first two columns
    ball2[, 1:2] <- scale(ball2[, 1:2]) 
    
    yM <- as.matrix(ball2[, 3])
    XU <- as.matrix(ball2[, 1:2])
    
    # Calculate original training data parameters
    mu_x1 <- mean(ball2_0$x1)
    sd_x1 <- sd(ball2_0$x1)
    mu_x2 <- mean(ball2_0$x2)
    sd_x2 <- sd(ball2_0$x2)
    
    # Fit Gaussian Process
    yMhat <- newGPsep(XU, yM, d = 0.1, g = 1e-7, dK = TRUE)
    mle <- mleGPsep(yMhat, tmin = 0.0001, tmax = 10)
    
    obj_min_indx <- vector(mode = "numeric", length = 10000)
    lhd_save <- list()
    results <- list()
    obj_min_x <- list()
    
    for (rep in 1:10000) {
      
      # Generate Latin Hypercube Sample
      lhd_unit <- randomLHS(n = 500, k = 2)
      
      # Map to original scale
      lhd_unit[, 1] <- lhd_unit[, 1] * (max(ball2_0$x1) - min(ball2_0$x1)) + min(ball2_0$x1)
      lhd_unit[, 2] <- lhd_unit[, 2] * (max(ball2_0$x2) - min(ball2_0$x2)) + min(ball2_0$x2)
      
      # Standardize using training data parameters
      lhd_unit[, 1] <- (lhd_unit[, 1] - mu_x1) / sd_x1
      lhd_unit[, 2] <- (lhd_unit[, 2] - mu_x2) / sd_x2
      
      xnew_dm <- as.matrix(lhd_unit)
      lhd_save <- append(lhd_save, list(xnew_dm))
      
      pred <- predGPsep(yMhat, lhd_unit, lite = FALSE)
      
      # Ensure positive definite covariance matrix
      Sigma_pd <- as.matrix(nearPD(pred$Sigma)$mat)
      obj <- rmvnorm(1, mean = pred$mean, sigma = Sigma_pd)^2
      
      # Get index of minimum objective
      obj_min_indx[rep] <- which.min(obj)
      obj_min_x <- append(obj_min_x, list(xnew_dm[obj_min_indx[rep], ]))
    }
    
    # Extract posterior x coordinates
    post_x <- t(matrix(unlist(obj_min_x), ncol = 10000, nrow = 2))
    x1 <- post_x[, 1]
    x2 <- post_x[, 2]
    
    # Recover original scale
    X_ori <- ball2_0[, 1:2]
    x1_recov <- x1 * sd(X_ori[, 1]) + mean(X_ori[, 1])
    x2_recov <- x2 * sd(X_ori[, 2]) + mean(X_ori[, 2])
    
    sd_x1 <- sd(x1_recov)
    sd_x2 <- sd(x2_recov)
    med_x1 <- median(x1_recov)
    med_x2 <- median(x2_recov)
    
    # Free C-side memory to prevent RAM crashes (CRITICAL FIX)
    deleteGPsep(yMhat)
    
    # Save workspace
    save_path <- paste0('GP_', filename, '_u100rep100lhd500.RData')
    save.image(file = save_path)
  }
}
