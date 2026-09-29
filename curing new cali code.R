setwd("~/Downloads/calibration")
library(MASS)

# 1. 读取并拟合高保真替代模型 (L_Q33)
filenameF = "L_Q33_complete"
GenerateF = read.csv(paste0(filenameF, ".csv"))

GenerateF = GenerateF[, c("T1", "T2", "t1", "t2", "u")]
scaled_GenerateF <- scale(GenerateF[, c("T1", "T2", "t1", "t2")])
ball2_scaled <- data.frame(scaled_GenerateF, y = GenerateF$u)

poly_model <- lm(y ~ T1 + T2 + t1 + t2 + 
                   I(T1^2) + I(T2^2) + I(t1^2) + I(t2^2) +
                   T1:T2 + T1:t1 + T1:t2 + T2:t1 + T2:t2 + t1:t2, 
                 data = ball2_scaled)

stepwise_model <- stepAIC(poly_model, direction = "both", trace = FALSE)

# 保存训练数据的标准化参数
train_means_F <- attr(scaled_GenerateF, "scaled:center")
train_sds_F   <- attr(scaled_GenerateF, "scaled:scale")


# ==================== 关键修改区域 ====================
# 2. 直接读取 200 点的低保真数据集，确保嵌套关系 (X_H ⊂ X_L)
low_fidelity_200 <- read.csv("44DoCS_noise0_lhd200_3.csv")

# 从 200 个点中提取 50 个点：
# 方法 A（推荐）：设置随机种子后随机抽取 50 行
set.seed(123) 
nested_indices <- sample(1:nrow(low_fidelity_200), 50)
lhd_original_scale <- low_fidelity_200[nested_indices, c("T1", "T2", "t1", "t2")]

# （或者方法 B：直接按顺序取前 50 行：lhd_original_scale <- low_fidelity_200[1:50, c("T1", "T2", "t1", "t2")]）
# ======================================================


# 3. 标准化选出的 50 个点并进行预测
lhd_scaled <- scale(lhd_original_scale, center = train_means_F, scale = train_sds_F)
lhd_scaled <- as.data.frame(lhd_scaled)

# 预测并添加 0.1 标准差的噪声
predictions <- predict(stepwise_model, newdata = lhd_scaled)
noise_sd <- 0.1  # 噪声标准差 0.1
set.seed(123)    # 保证噪声生成也可复现
predictions_noisy <- predictions + rnorm(length(predictions), mean = 0, sd = noise_sd)

# 4. 保存为新的嵌套高保真数据集
ball2 <- data.frame(lhd_original_scale, y = predictions_noisy)
write.csv(ball2, file = "uF_noise01_lhd50_nested_3.csv", row.names = FALSE)

cat("嵌套高保真数据集 uF_noise01_lhd50_nested_3.csv 已生成完成！\n")



###########################main code##########################
set.seed(20250734)
eps <- .Machine$double.eps

setwd("~/Downloads/calibration")

# 1. 修改为您刚刚生成的嵌套数据文件名！
filenameF <- "uF_noise01_lhd50_3"  # 高保真（50个点，嵌套）
filenameM <- "44DoCS_noise0_lhd200_3"     # 低保真（200个点）

ball4 <- read.csv(paste0(filenameF,".csv"))
ball2 <- read.csv(paste0(filenameM,".csv"))

ball2_0 <- ball2
ball4_0 <- ball4

# 合并数据用于标准化
X_oricombine = rbind(ball2_0[,1:4], ball4_0[,1:4])
mean_all = apply(X_oricombine, 2, mean)
sd_all = apply(X_oricombine, 2, sd)

ball2[,1:4] = scale(ball2[,1:4], center = mean_all, scale = sd_all)
ball4[,1:4] = scale(ball4[,1:4], center = mean_all, scale = sd_all)

XU <- as.matrix(ball2[,1:4])
XF <- as.matrix(ball4[,1:4])

u_bound=100
u <- seq(-u_bound, u_bound, length=1000) #grid

uhats_LOOsamples_summary=list()
uhats_LOOsamples=list()
uhat_summary=list()

# 因为数据里第5列是 y，所以 yi=1 刚好读取第5列
for (yi in 1:1) {
  
  yM <- as.matrix(ball2[,(4+yi)])
  yF <- as.matrix(ball4[,(4+yi)])
  
  yMhat <- newGPsep(XU, yM, d=0.1, g=1e-7, dK=TRUE)
  mle <- mleGPsep(yMhat, tmin=0.0001, tmax=10)
  
  formals(bhat.fit)$da <- darg(d=list(mle=TRUE), X=XF)
  formals(bhat.fit)$ga <- garg(g=list(mle=TRUE), y=yF)
  
  unll <- rep(NA, length(u))
  for(i in 1:length(u)) 
    unll[i] <- calib(u[i], XF, yF, yMhat, bhat.fit)$nll
  
  obj <- function(x, XF, yF, yMhat, fit) calib(x, XF, yF, yMhat, fit)$nll
  soln <- optimize(obj, lower=-u_bound, upper=u_bound, XF=XF, yF=yF, 
                   yMhat=yMhat, fit=bhat.fit, tol = .Machine$double.eps^0.5)
  uhat <- soln$minimum
  uhat_summary=append(uhat_summary,uhat)
  print(paste("y",yi,"uhat",uhat))
  
  uhats <- q1 <- q2 <- m <- s2 <- rep(NA, nrow(XF))
  for(i in 1:nrow(XF)) {
    train <- XF[-i,,drop=FALSE]
    test <- XF[i,,drop=FALSE]
    cp <- calib.pred(test, train, yF[-i], yMhat, bhat.fit)
    m[i] <- cp$mean
    s2[i] <- cp$s2 
    q1[i] <- cp$q1
    q2[i] <- cp$q2 
    uhats[i] <- cp$uhat
  }
  uhats_LOOsamples_summary=append(uhats_LOOsamples_summary,list(c(summary(uhats),sd=sd(uhats),quantile(uhats,0.025),quantile(uhats,0.975))))
  uhats_LOOsamples=append(uhats_LOOsamples,list(uhats))
  print(paste("y",yi,"uhats LOO samples"))
  print(c(summary(uhats),quantile(uhats,0.025),quantile(uhats,0.975)))
  
  deleteGPsep(yMhat)
}

#generate 100u for each y
unew100=list()
for (yi in 1:1) {
  kde <- density(uhats_LOOsamples[[yi]])
  unew100 <- append(unew100,list(sample(kde$x, size=100, prob=kde$y, replace=F)))
}

lhd_save=list()
obj_save=list()
results=list()
obj_min_x=list()

# 这个双重循环可能需要跑一段时间（100x100次采样优化）
for (ui in 1:100) {
  
  uhat=unew100[[1]][ui]
  obj_min_indx=vector(mode="numeric",length=100)
  
  for (rep in 1:100) {
    lhd_unit <- randomLHS(n = 200, k = 4)
    
    lhd_unit[,1] = lhd_unit[,1] * (max(X_oricombine[,1]) - min(X_oricombine[,1])) + min(X_oricombine[,1])
    lhd_unit[,2] = lhd_unit[,2] * (max(X_oricombine[,2]) - min(X_oricombine[,2])) + min(X_oricombine[,2])
    lhd_unit[,3] = lhd_unit[,3] * (max(X_oricombine[,3]) - min(X_oricombine[,3])) + min(X_oricombine[,3])
    lhd_unit[,4] = lhd_unit[,4] * (max(X_oricombine[,4]) - min(X_oricombine[,4])) + min(X_oricombine[,4])
    
    lhd_unit[,1:4] = scale(lhd_unit[,1:4], center = mean_all, scale = sd_all)
    
    xnew_dm=as.matrix(lhd_unit)
    lhd_save=append(lhd_save,list(xnew_dm))
    
    XF_new=rbind(XF,xnew_dm)
    
    yi=1
    yM <- as.matrix(ball2[,(4+yi)])
    yF <- as.matrix(ball4[,(4+yi)])
    
    yMhat <- newGPsep(XU, yM, d=0.1, g=1e-7, dK=TRUE)
    mle <- mleGPsep(yMhat, tmin=0.0001, tmax=10)
    
    formals(bhat.fit)$da <- darg(d=list(mle=TRUE), X=XF)
    formals(bhat.fit)$ga <- garg(g=list(mle=TRUE), y=yF)
    
    bhat <- calib(uhat, XF, yF, yMhat, bhat.fit, clean=FALSE)
    
    predictions <- predict_posterior(
      X_new = xnew_dm, XU = XU, XF = XF, yM = yM, yF = yF,
      uhat = uhat, mle = mle, bhat = bhat, yMhat = yMhat
    )
    
    deleteGPsep(yMhat)
    
    mu_ynew_y1=predictions$mean
    sigma2_ynew_y1=predictions$cov
    
    obj=rmvnorm(1,mean=mu_ynew_y1,sigma=sigma2_ynew_y1)^2
    obj_save=append(obj_save,list(obj))
    
    obj_min_indx[rep]=which(obj==min(obj))
    obj_min_x=append(obj_min_x,list(xnew_dm[obj_min_indx[rep],]))
  }
  
  results=append(results,list(obj_min_indx))
}

post_x=t(matrix(unlist(obj_min_x),ncol = 10000,nrow = 4))
x1=post_x[,1]
x2=post_x[,2]
x3=post_x[,3]
x4=post_x[,4]

x1_recov = x1 * sd_all[1] + mean_all[1]
x2_recov = x2 * sd_all[2] + mean_all[2]
x3_recov = x3 * sd_all[3] + mean_all[3]
x4_recov = x4 * sd_all[4] + mean_all[4]

med_x1=median(x1_recov)
med_x2=median(x2_recov)
med_x3=median(x3_recov)
med_x4=median(x4_recov)

# 保存工作区
save_path=paste0('Results_NestedData_Calibration.RData')
save.image(file=save_path)

# 画图：修改了标签为真实的物理含义 (T1, T2, t1, t2)
par(mfrow=c(2,2)) 
hist(x1_recov, freq=FALSE, main="", 
     xlab=expression(T[1]~"("*degree*C*")"), ylab="density", xlim=c(125, 180),
     breaks=seq(125, 180, by=2))
lines(density(x1_recov), col="red", lwd=1)
abline(v=med_x1, col="red", lwd=2, lty=2) 
text(x = med_x1 + 0.01, y = max(density(x1_recov)$y)*0.9, labels = round(med_x1, 2), col = "red", pos = 4)

hist(x2_recov, freq=FALSE, main="",
     xlab=expression(T[2]~"("*degree*C*")"), ylab="density", xlim=c(150, 180),
     breaks=seq(150, 180, by=1))
lines(density(x2_recov), col="red", lwd=1)
abline(v=med_x2, col="red", lwd=2, lty=2) 
text(x = med_x2 + 0.01, y = max(density(x2_recov)$y) * 0.9, labels = round(med_x2, 2), col = "red", pos = 4)

hist(x3_recov, freq=FALSE, main="",
     xlab=expression(t[1]~"(min)"), ylab="density", xlim=c(600, 6600),
     breaks=seq(600, 6600, by=200))
lines(density(x3_recov), col="red", lwd=1)
abline(v=med_x3, col="red", lwd=2, lty=2) 
text(x = med_x3 - 1500, y = max(density(x3_recov)$y) * 0.9, labels = round(med_x3, 2), col = "red", pos = 4)

hist(x4_recov, freq=FALSE, main="",
     xlab=expression(t[2]~"(min)"), ylab="density", xlim=c(7240, 12000),
     breaks=seq(7200, 12000, by=100))
lines(density(x4_recov), col="red", lwd=1)
abline(v=med_x4, col="red", lwd=2, lty=2) 
text(x = med_x4 + 0.01, y = max(density(x4_recov)$y) * 0.9, labels = round(med_x4, 2), col = "red", pos = 4)
par(mfrow=c(1,1))



# 画图：修改了标签，减少了间距，并去除了右侧多余的 y 轴标签
# mar 参数控制单个图的边界空白：c(bottom, left, top, right)
# 默认值是 c(5, 4, 4, 2) + 0.1，这里减小了 top 和 right 的空白使图更紧凑
par(mfrow=c(2,2), mar=c(4.5, 4.5, 1, 1)) 

# 图1：左上（保留 ylab）
hist(x1_recov, freq=FALSE, main="", 
     xlab=expression(T[1]~"("*degree*C*")"), ylab="density", xlim=c(125, 180),
     breaks=seq(125, 180, by=2))
lines(density(x1_recov), col="red", lwd=1)
abline(v=med_x1, col="red", lwd=2, lty=2) 
text(x = med_x1 + 0.01, y = max(density(x1_recov)$y)*0.9, labels = round(med_x1, 2), col = "red", pos = 4)

# 图2：右上（去除 ylab）
hist(x2_recov, freq=FALSE, main="",
     xlab=expression(T[2]~"("*degree*C*")"), ylab="", xlim=c(150, 180),
     breaks=seq(150, 180, by=1))
lines(density(x2_recov), col="red", lwd=1)
abline(v=med_x2, col="red", lwd=2, lty=2) 
text(x = med_x2 + 0.01, y = max(density(x2_recov)$y) * 0.9, labels = round(med_x2, 2), col = "red", pos = 4)

# 图3：左下（保留 ylab）
hist(x3_recov, freq=FALSE, main="",
     xlab=expression(t[1]~"(min)"), ylab="density", xlim=c(600, 6600),
     breaks=seq(600, 6600, by=200))
lines(density(x3_recov), col="red", lwd=1)
abline(v=med_x3, col="red", lwd=2, lty=2) 
text(x = med_x3 - 1500, y = max(density(x3_recov)$y) * 0.9, labels = round(med_x3, 2), col = "red", pos = 4)

# 图4：右下（去除 ylab）
hist(x4_recov, freq=FALSE, main="",
     xlab=expression(t[2]~"(min)"), ylab="", xlim=c(7240, 12000),
     breaks=seq(7200, 12000, by=100))
lines(density(x4_recov), col="red", lwd=1)
abline(v=med_x4, col="red", lwd=2, lty=2) 
text(x = med_x4 + 0.01, y = max(density(x4_recov)$y) * 0.9, labels = round(med_x4, 2), col = "red", pos = 4)

# 恢复默认的画布设置，以免影响后续其他绘图
par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)

