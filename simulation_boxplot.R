library(ggplot2)
library(tidyr)
library(gridExtra) 

# ==========================================
# Part 1: RMSE Boxplots
# ==========================================

# Read and process the first dataset (MSE)
result_df <- read.csv("x1_mse.csv", header = TRUE)
colnames(result_df) <- c("Low-fidelity", "High-fidelity", "Multi-fidelity Calibration")
result_df[which(result_df[,2] > 0.3), 2] <- median(result_df[,2])

# Convert MSE to RMSE
result_df <- sqrt(result_df)

# Read and process the second dataset (MSE)
result_df2 <- read.csv("x2_mse.csv", header = TRUE)
colnames(result_df2) <- c("Low-fidelity", "High-fidelity", "Multi-fidelity Calibration")
result_df2[which(result_df2[,2] > 0.3), 2] <- median(result_df2[,2])

# Convert MSE to RMSE
result_df2 <- sqrt(result_df2)

# Calculate global Y-axis limits for RMSE plots
min_y <- min(c(min(result_df, na.rm = TRUE), min(result_df2, na.rm = TRUE)))
max_y <- max(c(max(result_df, na.rm = TRUE), max(result_df2, na.rm = TRUE)))

# Convert to long format
result_long1 <- pivot_longer(result_df, 
                             cols = everything(), 
                             names_to = "Method", 
                             values_to = "RMSE")
result_long1$Dataset <- "X1"

result_long2 <- pivot_longer(result_df2, 
                             cols = everything(), 
                             names_to = "Method", 
                             values_to = "RMSE")
result_long2$Dataset <- "X2"

# Fix the order of the methods on the X-axis
fixed_order <- c("Multi-fidelity Calibration", "High-fidelity", "Low-fidelity")
result_long1$Method <- factor(result_long1$Method, levels = fixed_order)
result_long2$Method <- factor(result_long2$Method, levels = fixed_order)

# Define consistent color mapping
fill_colors <- c("Multi-fidelity Calibration" = "lightgreen", 
                 "High-fidelity" = "orange", 
                 "Low-fidelity" = "lightblue")

# Plot RMSE boxplots
p1 <- ggplot(result_long1, aes(x = Method, y = RMSE, fill = Method)) +
  geom_boxplot() +
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(limits = c(min_y, max_y)) +  
  labs(title = "", x = "x1", y = "RMSE") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(result_long2, aes(x = Method, y = RMSE, fill = Method)) +
  geom_boxplot() +
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(limits = c(min_y, max_y)) +  
  labs(title = "", x = "x2", y = "") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Arrange RMSE plots side-by-side
grid.arrange(p1, p2, ncol = 2)

# ==========================================
# Part 2: Bias Boxplots
# ==========================================

# Read median data
result_df <- read.csv("x1_med.csv", header = TRUE)
colnames(result_df) <- c("Low-fidelity", "High-fidelity", "Multi-fidelity Calibration")

result_df2 <- read.csv("x2_med.csv", header = TRUE)
colnames(result_df2) <- c("Low-fidelity", "High-fidelity", "Multi-fidelity Calibration")

# Calculate Bias = Estimated - True Value (x1 true = -0.8, x2 true = 0.4)
bias_df1 <- result_df - (-0.8)
bias_df2 <- result_df2 - 0.4

# Convert Bias data to long format
result_long1 <- pivot_longer(bias_df1, 
                             cols = everything(), 
                             names_to = "Method", 
                             values_to = "Bias")
result_long1$Dataset <- "X1"

result_long2 <- pivot_longer(bias_df2, 
                             cols = everything(), 
                             names_to = "Method", 
                             values_to = "Bias")
result_long2$Dataset <- "X2"

# Apply the fixed X-axis order to the Bias datasets
result_long1$Method <- factor(result_long1$Method, levels = fixed_order)
result_long2$Method <- factor(result_long2$Method, levels = fixed_order)

# Plot Bias boxplots with standard limits (-0.8 to 0.8)
p1 <- ggplot(result_long1, aes(x = Method, y = Bias, fill = Method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7, linewidth = 0.8) +
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(limits = c(-0.8, 0.8)) +  
  labs(title = "", x = "x1", y = "Bias") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(result_long2, aes(x = Method, y = Bias, fill = Method)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7, linewidth = 0.8) +
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(limits = c(-0.8, 0.8)) +  
  labs(title = "", x = "x2", y = "") +     
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Arrange Bias plots side-by-side
grid.arrange(p1, p2, ncol = 2)
