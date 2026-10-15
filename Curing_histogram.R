load("~/Downloads/cali_2025618F50M200_p5_noise02_1_u100rep100lhd500.RData")

# Plot results
par(mfrow = c(1, 2)) 
hist(x1_recov, freq = FALSE, main = "", 
     xlab = "x1", ylab = "density", xlim = c(-1, 0),
     breaks = seq(-1.0, 1, by = 0.02))
lines(density(x1_recov), col = "red", lwd = 1)
abline(v = med_x1, col = "red", lwd = 2, lty = 2)  # Mark x1 median
text(
  x = med_x1 + 0.01,  
  y = max(density(x1_recov)$y) * 0.9,
  labels = round(med_x1, 3),
  col = "red",
  pos = 4  
)

hist(x2_recov, freq = FALSE, main = "",
     xlab = "x2", ylab = "density", xlim = c(0, 1),
     breaks = seq(-1, 1, by = 0.02))
lines(density(x2_recov), col = "red", lwd = 1)
abline(v = med_x2, col = "red", lwd = 2, lty = 2)  # Mark x2 median
text(
  x = med_x2 + 0.01,  
  y = max(density(x2_recov)$y) * 0.95,
  labels = round(med_x2, 3),
  col = "red",
  pos = 4  
)
par(mfrow = c(1, 1))

