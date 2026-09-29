setwd("/Users/fan/Library/Mobile Documents/com~apple~CloudDocs/Documents/backup/calibration")
library(lhs)
set.seed(2025619)

# 如果没有这个文件夹，先自动创建，防止报错
if (!dir.exists("Calidata")) dir.create("Calidata")

# 修复3：添加50次外层循环，生成论文5.4节所需的50个独立数据集
for (i in 1:50) {
  
  # 修复1：实现嵌套设计 (Nested Design)，回应审稿人
  lhd1_raw <- randomLHS(200, 2)
  # 从200个低保真点中，随机抽取50个点作为高保真点，确保 X_H 是 X_L 的精确子集
  lhd2_raw <- lhd1_raw[sample(nrow(lhd1_raw), 50), ] 
  
  lhd1 <- as.data.frame(lhd1_raw)
  lhd2 <- as.data.frame(lhd2_raw)
  colnames(lhd1) <- c("x1", "x2")
  colnames(lhd2) <- c("x1", "x2")
  
  lhd1_transformed <- lhd1 * 2 - 1
  lhd2_transformed <- lhd2 * 2 - 1
  
  as_center <- c(-0.6, 0.2)
  ap_center <- c(-0.8, 0.4)
  
  # 计算第一个LHD(低保真)的ys
  ys <- rowSums((lhd1_transformed - rep(as_center, each = nrow(lhd1_transformed)))^2) + 5
  lhd1_transformed$ys = ys
  
  # 计算第二个LHD(高保真)的yp
  sum_term <- rowSums((lhd2_transformed - rep(ap_center, each = nrow(lhd2_transformed)))^2) + 5
  
  # 修复2：将噪声标准差从 sd=0.2 改为与论文一致的 sd=0.02
  error_term <- rnorm(nrow(lhd2_transformed), mean = 0, sd = 0.02)
  yp <- sum_term + error_term 
  lhd2_transformed$yp <- yp
  
  # 写入CSV，循环保存50个文件
  write.csv(lhd2_transformed, paste0("Calidata/NEW2025618F_p5_noise02_", i, ".csv"), row.names = FALSE)
  write.csv(lhd1_transformed, paste0("Calidata/NEW2025618M_p5_noise02_", i, ".csv"), row.names = FALSE)
}
