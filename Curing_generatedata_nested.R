#setwd("~/Downloads")
library(MASS)
library(knitr)
library(laGP)
library(plgp)
library(lhs)

##################################################
# nlhd.R: Generate nested Latin hypercube designs
##################################################

nlhd <- function(S, k) {
  k <- as.integer(k) 
  u <- length(S)    # number of layers 
  n <- prod(S)      # total number of runs
  
  # Initializations
  pmat <- matrix(0, nrow = n, ncol = k)
  xmat <- matrix(0, nrow = n, ncol = k)

  for (p in 1:k) {
    C <- 0  
    D <- rep(0, n)
    oldm <- 0
    
    for (i in 1:u) {
      newm <- prod(S[1:i])
      t <- n / newm
      if (i > 1) { 
        oldm <- prod(S[1:(i-1)])
        C <- ceiling(D[1:oldm] / t)
      }
      
      PI <- setdiff(1:newm, C)
      PI <- sample(PI, newm - oldm)
      
      if (t == 1) {
        D[(oldm+1):newm] <- sample(PI, newm - oldm)
      } else {
        for (j in 1:(newm - oldm)) {
          D[j + oldm] <- sample(((PI[j]-1)*t + 1):(PI[j]*t), 1)
        }
      }
    }
    U <- runif(n)
    pmat[, p] <- D
    xmat[, p] <- (D - U) / n
  }
  output <- list(xmat = xmat, pmat = pmat)
  return(output)
}

set.seed(10092026)

nlhd_result <- nlhd(S = c(50, 4), k = 4)
nested_xmat <- nlhd_result$xmat
lhd_50_nested <- nested_xmat[1:50, ]  
lhd_200_total <- nested_xmat

# Read and filter datasets
GenerateM <- read.csv("L_Q44_complete.csv")
GenerateM <- GenerateM[GenerateM$DoC > 0.96, c("T1", "T2", "t1", "t2", "DoC")]

GenerateF <- read.csv("L_Q33_complete.csv")
GenerateF <- GenerateF[GenerateF$DoC > 0.96, c("T1", "T2", "t1", "t2", "u")] 

# Combine independent variables (first 4 columns)
combined_X <- rbind(GenerateM[, 1:4], GenerateF[, 1:4])

# Calculate global physical boundaries
global_bounds <- data.frame(
  min = apply(combined_X, 2, min),
  max = apply(combined_X, 2, max)
)

##################################################
# Low Fidelity Data Generation
##################################################

scaled_GenerateM <- scale(GenerateM[, c("T1", "T2", "t1", "t2")])
ball4_scaled <- data.frame(scaled_GenerateM, y = GenerateM$DoC)

poly_model <- lm(y ~ T1 + T2 + t1 + t2 +                    
                   I(T1^2) + I(T2^2) + I(t1^2) + I(t2^2) +                   
                   T1:T2 + T1:t1 + T1:t2 + T2:t1 + T2:t2 + t1:t2,                  
                 data = ball4_scaled)
summary(poly_model)$r.squared

stepwise_model <- stepAIC(poly_model, direction = "both", trace = TRUE)
summary(stepwise_model)$r.squared

# Save scaling parameters for training data
train_means_M <- attr(scaled_GenerateM, "scaled:center")
train_sds_M <- attr(scaled_GenerateM, "scaled:scale")

lhd_original_scale <- lhd_200_total
colnames(lhd_original_scale) <- c("T1", "T2", "t1", "t2")

# Transform to appropriate physical scale
for (i in 1:4) {
  lhd_original_scale[, i] <- lhd_original_scale[, i] * 
    (global_bounds$max[i] - global_bounds$min[i]) + global_bounds$min[i]
}

# Standardize
lhd_scaled <- scale(lhd_original_scale, center = train_means_M, scale = train_sds_M)
lhd_scaled <- as.data.frame(lhd_scaled)

# Generate predictions with noise
predictions <- predict(stepwise_model, newdata = lhd_scaled)
noise_sd <- 0  # Adjust noise standard deviation as needed
predictions_noisy <- predictions + rnorm(length(predictions), mean = 0, sd = noise_sd)

ball4 <- data.frame(lhd_original_scale, y = predictions_noisy)
write.csv(ball4, file = "Curing_44DoCS_noise0_lhd200.csv", row.names = FALSE)

##################################################
# High Fidelity Data Generation
##################################################

scaled_GenerateF <- scale(GenerateF[, c("T1", "T2", "t1", "t2")])
ball2_scaled <- data.frame(scaled_GenerateF, y = GenerateF$u)

poly_model <- lm(y ~ T1 + T2 + t1 + t2 +                   
                   I(T1^2) + I(T2^2) + I(t1^2) + I(t2^2) +                  
                   T1:T2 + T1:t1 + T1:t2 + T2:t1 + T2:t2 + t1:t2,                
                 data = ball2_scaled)

stepwise_model <- stepAIC(poly_model, direction = "both", trace = FALSE)

# Save scaling parameters for training data
train_means_F <- attr(scaled_GenerateF, "scaled:center")
train_sds_F   <- attr(scaled_GenerateF, "scaled:scale")

lhd_original_scale <- lhd_50_nested
colnames(lhd_original_scale) <- c("T1", "T2", "t1", "t2")

# Transform to appropriate physical scale
for (i in 1:4) {
  lhd_original_scale[, i] <- lhd_original_scale[, i] * 
    (global_bounds$max[i] - global_bounds$min[i]) + global_bounds$min[i]
}

# Standardize
lhd_scaled <- scale(lhd_original_scale, center = train_means_F, scale = train_sds_F)
lhd_scaled <- as.data.frame(lhd_scaled)

# Generate predictions with noise
predictions <- predict(stepwise_model, newdata = lhd_scaled)
noise_sd <- 0.1  # Adjust noise standard deviation as needed
predictions_noisy <- predictions + rnorm(length(predictions), mean = 0, sd = noise_sd)

ball2 <- data.frame(lhd_original_scale, y = predictions_noisy)
write.csv(ball2, file = "Curing_33uF_noise01_lhd50_nested.csv", row.names = FALSE)
