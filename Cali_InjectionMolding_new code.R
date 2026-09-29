library(knitr)
library(laGP)
library(plgp)
library(lhs)

###########################functions########################
bhat.fit <- function(X, yF, Ym, da, ga, clean=TRUE)
{
  bhat <- newGPsep(X, yF - Ym, d=da$start, g=ga$start, dK=TRUE)
  if(ga$mle) cmle <- mleGPsep(bhat, param="both", tmin=c(da$min, ga$min),
                              tmax=c(da$max, ga$max), ab=c(da$ab, ga$ab))
  else cmle <- mleGPsep(bhat, tmin=da$min, tmax=da$max, ab=da$ab)
  cmle$nll <- - llikGPsep(bhat, dab=da$ab, gab=ga$ab)
  if(clean) deleteGPsep(bhat)
  else { cmle$gp <- bhat; cmle$gptype <- "sep" }        
  return(cmle)
}
calib <- function(u, XF, yF, yMhat, fit, clean=TRUE)
{
  #XFu <- cbind(XF, matrix(rep(u, nrow(XF)), ncol=length(u), byrow=TRUE))
  XFu <- XF
  Ym <- u*predGPsep(yMhat, XFu, lite=TRUE)$mean
  #Ym <- u*sqrt(2*XF)
  cmle <- fit(XF, yF, Ym, clean=clean)
  return(cmle)
}
calib.pred <- function(XX, XF, yF, yMhat, fit)
{
  soln <- optimize(obj, lower=-u_bound, upper=u_bound, XF=XF, yF=yF, 
                   yMhat=yMhat, fit=fit)
  bhat <- calib(soln$minimum, XF, yF, yMhat, fit, clean=FALSE)
  if(bhat$gptype == "sep") pb <- predGPsep(bhat$gp, XX, lite=TRUE)
  else pb <- predGP(bhat$gp, XX, lite=TRUE)
  pm <- predGPsep(yMhat, XX, lite=TRUE)
  m <- soln$minimum*pm$mean + pb$mean
  s2 <- soln$minimum^2*pm$s2 + pb$s2 
  q1 <- qnorm(0.95, m, sqrt(s2))
  q2 <- qnorm(0.05, m, sqrt(s2))
  if(bhat$gptype == "sep") deleteGPsep(bhat$gp)
  else deleteGP(bhat$gp)
  return(list(mean=m, s2=s2, q1=q1, q2=q2, uhat=soln$minimum))
}
ViVldet <- function(u, XU, XF, SigmaM, tau2M, mle, SigmaB)
{
  ## build blocks
  #XFu <- cbind(XF, u)
  SMXFu <- tau2M * covar.sep(XU, XF, d=mle$d, g=0)
  SMu <- tau2M * covar.sep(XF, d=mle$d, g=eps)
  
  ## build V from blocks
  V <- cbind(SigmaM, u*SMXFu)
  V <- rbind(V, cbind(u*t(SMXFu), u^2*SMu + SigmaB))
  
  ## return inverse and determinant
  ## (improvements possible with partitioned inverse equations)
  return(list(V=V,
              inv=solve(V), 
              ldet=as.numeric(determinant(V, log=TRUE)$modulus)))
}
predict_posterior <- function(X_new, XU, XF, yM, yF, uhat, mle, bhat, yMhat) {
  # 确保eps已定义
  eps <- sqrt(.Machine$double.eps)
  
  # 计算计算机模型的协方差矩阵
  KM <- covar.sep(XU, d=mle$d, g=eps)
  tau2M <- drop(t(yM) %*% solve(KM) %*% yM / length(yM))
  SigmaM <- tau2M * KM
  
  # 计算偏差函数的协方差矩阵
  KB <- covar.sep(XF, d=bhat$theta[1:ncol(XF)], g=bhat$theta[ncol(XF)+1])
  Ym <- uhat * predGPsep(yMhat, XF, lite=TRUE)$mean
  YmYm <- yF - Ym
  tau2B <- drop(t(YmYm) %*% solve(KB) %*% YmYm / length(YmYm))
  SigmaB <- tau2B * KB
  
  # 获取V矩阵的逆
  VV <- ViVldet(uhat, XU, XF, SigmaM, tau2M, mle, SigmaB)
  Vinv <- VV$inv
  
  # 计算交叉协方差
  k_new_M <- tau2M * covar.sep(XU, X_new, d=mle$d, g=0)  # 维度: nrow(XU) x nrow(X_new)
  k_new_F <- tau2M * covar.sep(XF, X_new, d=mle$d, g=0)  # 维度: nrow(XF) x nrow(X_new)
  
  # 构建联合协方差向量 - 注意uhat的应用
  k_star_block1_corrected = uhat * k_new_M
  k_delta_XF_Xnew = tau2B * covar.sep(XF, X_new, d=bhat$theta[1:ncol(XF)], g=bhat$theta[ncol(XF)+1])
  k_star_block2_corrected = uhat^2 * k_new_F + k_delta_XF_Xnew
  k_star = rbind(k_star_block1_corrected, k_star_block2_corrected)
  #k_star <- rbind(k_new_M, uhat * k_new_F)  # 维度: (nrow(XU) + nrow(XF)) x nrow(X_new)
  
  # 计算X_new处的自协方差矩阵 (包括点之间的协方差)
  #k_new_new <- tau2M * covar.sep(X_new, d=mle$d, g=eps)  # 维度: nrow(X_new) x nrow(X_new)
  k_new_new <- uhat^2 * tau2M * covar.sep(X_new, d=mle$d, g=eps) + tau2B * covar.sep(X_new, d=bhat$theta[1:ncol(XF)], g=bhat$theta[ncol(XF)+1])
  # 合并观测数据
  y_obs <- c(yM, yF)
  
  # 计算后验均值向量
  mu_post <- t(k_star) %*% Vinv %*% y_obs  # 维度: nrow(X_new) x 1
  
  # 计算完整的后验协方差矩阵
  Sigma_post <- k_new_new - t(k_star) %*% Vinv %*% k_star  # 维度: nrow(X_new) x nrow(X_new)
  
  # 确保协方差矩阵是对称正定的
  Sigma_post <- (Sigma_post + t(Sigma_post))/2
  diag(Sigma_post) <- pmax(diag(Sigma_post), eps)
  
  return(list(
    mean = mu_post,
    cov = Sigma_post,
    variance = diag(Sigma_post),
    lower_95 = mu_post - 1.96 * sqrt(diag(Sigma_post)),
    upper_95 = mu_post + 1.96 * sqrt(diag(Sigma_post))
  ))
}

###########################main code##########################
set.seed(20250815)
eps <- .Machine$double.eps

setwd("~/Downloads/calibration")
ball2 <- read.csv("sim.csv")
ball4 <- read.csv("exp.csv")
#ball2 <- read.csv("sim n55.csv")
#ball4 <- read.csv("sim n27.csv")
ball2$Hor=abs(ball2$Horizontal_Left)+abs(ball2$Horizontal_Right)
ball2$Ver=abs(ball2$Vertical_Up)+abs(ball2$Vertical_Down)
ball4$Hor=abs(ball4$Horizontal_Left)+abs(ball4$Horizontal_Right)
ball4$Ver=abs(ball4$Vertical_Up)+abs(ball4$Vertical_Down)

ball2_0 <- ball2
ball4_0 <- ball4

X_oricombine = rbind(ball2_0[,1:4], ball4_0[,1:4])
mean_all = apply(X_oricombine, 2, mean)
sd_all = apply(X_oricombine, 2, sd)

ball2[,1:4] = scale(ball2[,1:4], center = mean_all, scale = sd_all)
ball4[,1:4] = scale(ball4[,1:4], center = mean_all, scale = sd_all)

XU <- as.matrix(ball2[,1:4])
XF <- as.matrix(ball4[,1:4])

u_bound=10
u <- seq(-u_bound, u_bound, length=1000) #grid

uhats_LOOsamples_summary=list()
uhats_LOOsamples=list()
uhat_summary=list()

for (yi in 1:4) {
  
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
for (yi in 1:4) {
  kde <- density(uhats_LOOsamples[[yi]])
  unew100 <- append(unew100,list(sample(kde$x, size=100, prob=kde$y, replace=F)))
}

lhd_save=list()
obj_save=list()
results=list()
obj_min_x=list()

for (ui in 1:100) {
  obj_min_indx=vector(mode="numeric",length=100)
  
  for (rep in 1:100) {
    #generate new test X
    lhd_unit <- randomLHS(n = 200, k = 4)
    
    # 1. 先将 lhd_unit 映射到原始区间
    lhd_unit[,1] = lhd_unit[,1] * (max(X_oricombine[,1]) - min(X_oricombine[,1])) + min(X_oricombine[,1])
    lhd_unit[,2] = lhd_unit[,2] * (max(X_oricombine[,2]) - min(X_oricombine[,2])) + min(X_oricombine[,2])
    lhd_unit[,3] = lhd_unit[,3] * (max(X_oricombine[,3]) - min(X_oricombine[,3])) + min(X_oricombine[,3])
    lhd_unit[,4] = lhd_unit[,4] * (max(X_oricombine[,4]) - min(X_oricombine[,4])) + min(X_oricombine[,4])
    
    # 2. 用合并数据的均值和标准差做标准化
    lhd_unit[,1:4] = scale(lhd_unit[,1:4], center = mean_all, scale = sd_all)

    # generate new x
    xnew_dm=as.matrix(lhd_unit)
    lhd_save=append(lhd_save,list(xnew_dm))
    
    obj=vector(mode="numeric",length=dim(xnew_dm)[1])
    
    XF_new=rbind(XF,xnew_dm)
    
    yi=1
    uhat=unew100[[1]][ui]
    yM <- as.matrix(ball2[,(4+yi)])
    yF <- as.matrix(ball4[,(4+yi)])
    
    yMhat <- newGPsep(XU, yM, d=0.1, g=1e-7, dK=TRUE)
    mle <- mleGPsep(yMhat, tmin=0.0001, tmax=10)
    
    formals(bhat.fit)$da <- darg(d=list(mle=TRUE), X=XF)
    formals(bhat.fit)$ga <- garg(g=list(mle=TRUE), y=yF)
    
    bhat <- calib(uhat, XF, yF, yMhat, bhat.fit, clean=FALSE)
    
    predictions <- predict_posterior(
      X_new = xnew_dm,
      XU = XU,
      XF = XF,
      yM = yM,
      yF = yF,
      uhat = uhat,
      mle = mle,
      bhat = bhat,
      yMhat = yMhat
    )
    
    deleteGPsep(yMhat)
    
    #mean
    mu_ynew_y1=predictions$mean
    #var
    sigma2_ynew_y1=sigma=predictions$cov
    
 
    yi=2
    uhat=unew100[[2]][ui]
    yM <- as.matrix(ball2[,(4+yi)])
    yF <- as.matrix(ball4[,(4+yi)])
    
    yMhat <- newGPsep(XU, yM, d=0.1, g=1e-7, dK=TRUE)
    mle <- mleGPsep(yMhat, tmin=0.0001, tmax=10)
    
    formals(bhat.fit)$da <- darg(d=list(mle=TRUE), X=XF)
    formals(bhat.fit)$ga <- garg(g=list(mle=TRUE), y=yF)
    
    bhat <- calib(uhat, XF, yF, yMhat, bhat.fit, clean=FALSE)
    
    predictions <- predict_posterior(
      X_new = xnew_dm,
      XU = XU,
      XF = XF,
      yM = yM,
      yF = yF,
      uhat = uhat,
      mle = mle,
      bhat = bhat,
      yMhat = yMhat
    )
    
    deleteGPsep(yMhat)
    
    #mean
    mu_ynew_y2=predictions$mean
    #var
    sigma2_ynew_y2=sigma=predictions$cov
    
 
    yi=3
    uhat=unew100[[3]][ui]
    yM <- as.matrix(ball2[,(4+yi)])
    yF <- as.matrix(ball4[,(4+yi)])
    
    yMhat <- newGPsep(XU, yM, d=0.1, g=1e-7, dK=TRUE)
    mle <- mleGPsep(yMhat, tmin=0.0001, tmax=10)
    
    formals(bhat.fit)$da <- darg(d=list(mle=TRUE), X=XF)
    formals(bhat.fit)$ga <- garg(g=list(mle=TRUE), y=yF)
    
    bhat <- calib(uhat, XF, yF, yMhat, bhat.fit, clean=FALSE)
    
    predictions <- predict_posterior(
      X_new = xnew_dm,
      XU = XU,
      XF = XF,
      yM = yM,
      yF = yF,
      uhat = uhat,
      mle = mle,
      bhat = bhat,
      yMhat = yMhat
    )
    
    deleteGPsep(yMhat)
    
    #mean
    mu_ynew_y3=predictions$mean
    #var
    sigma2_ynew_y3=sigma=predictions$cov
    
 
    yi=4
    uhat=unew100[[4]][ui]
    yM <- as.matrix(ball2[,(4+yi)])
    yF <- as.matrix(ball4[,(4+yi)])
    
    yMhat <- newGPsep(XU, yM, d=0.1, g=1e-7, dK=TRUE)
    mle <- mleGPsep(yMhat, tmin=0.0001, tmax=10)
    
    formals(bhat.fit)$da <- darg(d=list(mle=TRUE), X=XF)
    formals(bhat.fit)$ga <- garg(g=list(mle=TRUE), y=yF)
    
    bhat <- calib(uhat, XF, yF, yMhat, bhat.fit, clean=FALSE)
    
    predictions <- predict_posterior(
      X_new = xnew_dm,
      XU = XU,
      XF = XF,
      yM = yM,
      yF = yF,
      uhat = uhat,
      mle = mle,
      bhat = bhat,
      yMhat = yMhat
    )
    
    deleteGPsep(yMhat)
    
    #mean
    mu_ynew_y4=predictions$mean
    #var
    sigma2_ynew_y4=sigma=predictions$cov
    
    
    obj=rmvnorm(1,mean=mu_ynew_y1,sigma=sigma2_ynew_y1)^2+
      rmvnorm(1,mean=mu_ynew_y2,sigma=sigma2_ynew_y2)^2+
      rmvnorm(1,mean=mu_ynew_y3,sigma=sigma2_ynew_y3)^2+
      rmvnorm(1,mean=mu_ynew_y4,sigma=sigma2_ynew_y4)^2
    #obj=mu_ynew_y1^2
    obj_save=append(obj_save,list(obj))

    obj_min_indx[rep]=which(obj==min(obj))
    obj_min_x=append(obj_min_x,list(xnew_dm[obj_min_indx[rep],]))
  }
  
  results=append(results,list(obj_min_indx))
}
save_path=paste0('NEW_InjectionMolding_cali_y1234.RData')
save.image(file=save_path)
load("~/Downloads/calibration/NEW_InjectionMolding_cali_y1234.RData")
post_x=t(matrix(unlist(obj_min_x),ncol = 10000,nrow = 4))
x1=post_x[,1]
x2=post_x[,2]
x3=post_x[,3]
x4=post_x[,4]

x1_recov = x1 * sd_all[1] + mean_all[1]
x2_recov = x2 * sd_all[2] + mean_all[2]
x3_recov = x3 * sd_all[3] + mean_all[3]
x4_recov = x4 * sd_all[4] + mean_all[4]

sd(x1_recov)
sd(x2_recov)
sd(x3_recov)
sd(x4_recov)

med_x1=median(x1_recov)
med_x2=median(x2_recov)
med_x3=median(x3_recov)
med_x4=median(x4_recov)

# 画图：调整布局使图像更紧凑，并去除了右侧列多余的 y 轴标签
# mar=c(bottom, left, top, right)，把 top 从默认 4.1 缩减到 1 以消除两行间的大空白
par(mfrow=c(2,2), mar=c(4.5, 4.5, 1, 1)) 

# 图1：x1（左上，保留 ylab）
hist(x1_recov, freq=FALSE, main="", 
     xlab="Mold.Temperature", ylab="density", xlim=c(30, 50),
     breaks=seq(30, 50, by=1))
lines(density(x1_recov), col="red", lwd=1)
abline(v=med_x1, col="red", lwd=2, lty=2) 
text(
  x = med_x1 -5, 
  y = max(density(x1_recov)$y)*0.9,
  labels = round(med_x1, 2),
  col = "red",
  pos = 4
)

# 图2：x2（右上，去除 ylab）
hist(x2_recov, freq=FALSE, main="",
     xlab="Injection.Speed", ylab="density", xlim=c(20, 70),
     breaks=seq(20, 70, by=2))
lines(density(x2_recov), col="red", lwd=1)
abline(v=med_x2, col="red", lwd=2, lty=2) 
text(
  x = med_x2 + 0.01, 
  y = max(density(x2_recov)$y) * 0.9,
  labels = round(med_x2, 2),
  col = "red",
  pos = 4
)

# 图3：x3（左下，保留 ylab）
hist(x3_recov, freq=FALSE, main="",
     xlab="Packing.Pressure", ylab="density", xlim=c(400, 600),
     breaks=seq(400, 600, by=10))
lines(density(x3_recov), col="red", lwd=1)
abline(v=med_x3, col="red", lwd=2, lty=2) 
text(
  x = med_x3 +1, 
  y = max(density(x3_recov)$y) * 0.9,
  labels = round(med_x3, 2),
  col = "red",
  pos = 4
)

# 图4：x4（右下，去除 ylab）
hist(x4_recov, freq=FALSE, main="",
     xlab="Packing.Time", ylab="density", xlim=c(1, 5),
     breaks=seq(1, 5, by=0.2))
lines(density(x4_recov), col="red", lwd=1)
abline(v=med_x4, col="red", lwd=2, lty=2) 
text(
  x = med_x4 + 0.01, 
  y = max(density(x4_recov)$y) * 0.9,
  labels = round(med_x4, 2),
  col = "red",
  pos = 4
)

# 恢复默认画布边距，避免影响后续绘图
par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)



# =======================================================
# Bivariate Analysis for Optimal X (Fixed & Simplified)
# =======================================================

# 1. 将恢复后的最优参数合并为一个数据框
X_recovered <- data.frame(
  Mold.Temp = x1_recov,
  Inj.Speed = x2_recov,
  Pack.Press = x3_recov,
  Pack.Time = x4_recov
)

# 2. 定义左下角面板：二维平滑散点图 (2D Density)
panel.smoothScatter <- function(x, y, ...) {
  smoothScatter(x, y, add = TRUE, nrpoints = 0, 
                colramp = colorRampPalette(c("white", "lightblue", "blue", "darkblue")))
}

# 3. 绘制双变量散点图矩阵
par(mfrow=c(1,1), mar=c(4, 4, 4, 2) + 0.1) # 恢复默认画布
pairs(X_recovered,
      lower.panel = panel.smoothScatter, # 左下角展示二维密度图
      upper.panel = NULL,                # 右上角留空，避免信息冗余
      labels = c("Mold Temp\nx1", "Inj Speed\nx2", "Pack Press\nx3", "Pack Time\nx4"),
      main = "Bivariate Joint Distribution of Optimal Parameters",
      gap = 0.5)


###############GPT6##############
## Figure 7 后新增：四个最优工艺输入的成对图
## obj_min_x[[i]] 是同一次优化选出的完整四维输入

stopifnot(
  length(obj_min_x) > 0L,
  all(lengths(obj_min_x) == 4L),
  length(mean_all) == 4L,
  length(sd_all) == 4L,
  all(is.finite(mean_all)),
  all(is.finite(sd_all)),
  all(sd_all > 0)
)

# 每一行对应同一次优化；与你原来的 post_x 行列排列一致
Xopt_std <- do.call(rbind, obj_min_x)
stopifnot(ncol(Xopt_std) == 4L, all(is.finite(Xopt_std)))

# 恢复原始工艺单位；列顺序须与 ball2[, 1:4] 一致
Xopt <- sweep(Xopt_std, 2, as.numeric(sd_all), "*")
Xopt <- sweep(Xopt,     2, as.numeric(mean_all), "+")

colnames(Xopt) <- c(
  "Mold temp.", "Injection speed",
  "Packing pressure", "Packing time"
)

cat("Xopt dimension:", nrow(Xopt), "x", ncol(Xopt), "\n")

# 左下角：成对散点。透明度有助于观察 10,000 个点的聚集区域
panel_points <- function(x, y, ...) {
  points(
    x, y, pch = 16, cex = 0.45,
    col = grDevices::adjustcolor("#1B4F72", alpha.f = 0.40)
  )
}

# 右上角：描述性的 Spearman 相关系数
panel_corr <- function(x, y, ...) {
  r <- suppressWarnings(cor(x, y, method = "spearman"))
  usr <- par("usr")
  on.exit(par(usr = usr))
  par(usr = c(0, 1, 0, 1))
  
  if (is.finite(r)) {
    text(0.5, 0.5, bquote(rho[S] == .(round(r, 2))), cex = 1.05)
  } else {
    text(0.5, 0.5, "rho = NA")
  }
}

png("Figure7_pairwise_xopt.png", width = 2400, height = 2400, res = 300)
par(mar = c(1, 1, 1, 1), oma = c(3, 3, 1, 1), cex.axis = 0.8)

pairs(
  Xopt,
  lower.panel = panel_points,
  upper.panel = panel_corr,
  cex.labels = 0.95,
  gap = 0.5
)

dev.off()

# 另输出相关矩阵，便于撰写正文
round(cor(Xopt, method = "spearman"), 2)
