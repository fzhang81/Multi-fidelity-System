setwd("~/Downloads")
library(lhs)
set.seed(2025619)

# Define parameter vectors
as <- c(-0.6, 0.2)
ap <- c(-0.8, 0.4)

# Generate 50 datasets
for (i in 1:50) {
  
  # Generate Latin Hypercube Designs
  lhd1 <- randomLHS(200, 2)
  lhd2 <- randomLHS(50, 2)
  
  # Convert to data frames and assign column names
  lhd1 <- as.data.frame(lhd1)
  lhd2 <- as.data.frame(lhd2)
  colnames(lhd1) <- c("x1", "x2")
  colnames(lhd2) <- c("x1", "x2")
  
  # Transform variables to [-1, 1] range
  lhd1_transformed <- lhd1 * 2 - 1
  lhd2_transformed <- lhd2 * 2 - 1
  
  # Calculate ys for the first LHD
  ys <- rowSums((lhd1_transformed - rep(as, each = nrow(lhd1_transformed)))^2) + 5
  lhd1_transformed$ys <- ys
  
  # Calculate yp for the second LHD
  sum_term <- rowSums((lhd2_transformed - rep(ap, each = nrow(lhd2_transformed)))^2) + 5
  # product_term <- 0.1 * lhd2_transformed$x1 * lhd2_transformed$x2
  error_term <- rnorm(nrow(lhd2_transformed), mean = 0, sd = 0.2)
  yp <- sum_term + error_term # + product_term
  lhd2_transformed$yp <- yp
  
  # Save generated data to CSV files
  write.csv(lhd2_transformed, paste0("Calidata/2025618F_p5_noise02_", i, ".csv"), row.names = FALSE)
  write.csv(lhd1_transformed, paste0("Calidata/2025618M_p5_noise02_", i, ".csv"), row.names = FALSE)
  
}

