library(knitr)
library(laGP)
library(plgp)
library(lhs)

########################### Functions ########################
bhat.fit <- function(X, yF, Ym, da, ga, clean=TRUE) {
  bhat <- newGPsep(X, yF - Ym, d=da$start, g=ga$start, dK=TRUE)
  if(ga$mle) {
    cmle <- mleGPsep(bhat, param="both", tmin=c(da$min, ga$min),
                     tmax=c(da$max, ga$max), ab=c(da$ab, ga$ab))
  } else {
    cmle <- mleGPsep(bhat, tmin=da$min, tmax=da$max, ab=da$ab)
  }
  cmle$nll <- - llikGPsep(bhat, dab=da$ab, gab=ga$ab)
  if(clean) deleteGPsep(bhat)
  else { cmle$gp <- bhat; cmle$gptype <- "sep" }        
  return(cmle)
}

calib <- function(u, XF, yF, yMhat, fit, clean=TRUE) {
  XFu <- XF
  Ym <- u * predGPsep(yMhat, XFu, lite=TRUE)$mean
  cmle <- fit(XF, yF, Ym, clean=clean)
  return(cmle)
}

calib.pred <- function(XX, XF, yF, yMhat, fit) {
  soln <- optimize(obj, lower=-u_bound, upper=u_bound, XF=XF, yF=yF, yMhat=yMhat, fit=fit)
  bhat <- calib(soln$minimum, XF, yF, yMhat, fit, clean=FALSE)
  if(bhat$gptype == "sep") pb <- predGPsep(bhat$gp, XX, lite=TRUE)
  else pb <- predGP(bhat$gp, XX, lite=TRUE)
  pm <- predGPsep(yMhat, XX, lite=TRUE)
  m <- soln$minimum * pm$mean + pb$mean
  s2 <- soln$minimum^2 * pm$s2 + pb$s2 
  q1 <- qnorm(0.95, m, sqrt(s2))
  q2 <- qnorm(0.05, m, sqrt(s2))
  if(bhat$gptype == "sep") deleteGPsep(bhat$gp)
  else deleteGP(bhat$gp)
  return(list(mean=m, s2=s2, q1=q1, q2=q2, uhat=soln$minimum))
}

ViVldet <- function(u, XU, XF, SigmaM, tau2M, mle, SigmaB) {
  eps <- sqrt(.Machine$double.eps)
  SMXFu <- tau2M * covar.sep(XU, XF, d=mle$d, g=0)
  SMu <- tau2M * covar.sep(XF, d=mle$d, g=eps)
  
  V <- cbind(SigmaM, u*SMXFu)
  V <- rbind(V, cbind(u*t(SMXFu), u^2*SMu + SigmaB))
  
  return(list(V=V, inv=solve(V), ldet=as.numeric(determinant(V, log=TRUE)$modulus)))
}

predict_posterior <- function(X_new, XU, XF, yM, yF, uhat, mle, bhat, yMhat) {
  eps <- sqrt(.Machine$double.eps)
  KM <- covar.sep(XU, d=mle$d, g=eps)
  tau2M <- drop(t(yM) %*% solve(KM) %*% yM / length(yM))
  SigmaM <- tau2M * KM
  
  KB <- covar.sep(XF, d=bhat$theta[1:ncol(XF)], g=bhat$theta[ncol(XF)+1])
  Ym <- uhat * predGPsep(yMhat, XF, lite=TRUE)$mean
  YmYm <- yF - Ym
  tau2B <- drop(t(YmYm) %*% solve(KB) %*% YmYm / length(YmYm))
  SigmaB <- tau2B * KB
  
  VV <- ViVldet(uhat, XU, XF, SigmaM, tau2M, mle, SigmaB)
  Vinv <- VV$inv
  
  k_new_M <- tau2M * covar.sep(XU, X_new, d=mle$d, g=0)
  k_new_F <- tau2M * covar.sep(XF, X_new, d=mle$d, g=0)
  
  k_star_block1_corrected = uhat * k_new_M
  k_delta_XF_Xnew = tau2B * covar.sep(XF, X_new, d=bhat$theta[1:ncol(XF)], g=bhat$theta[ncol(XF)+1])
  k_star_block2_corrected = uhat^2 * k_new_F + k_delta_XF_Xnew
  k_star = rbind(k_star_block1_corrected, k_star_block2_corrected)
  
  k_new_new <- uhat^2 * tau2M * covar.sep(X_new, d=mle$d, g=eps) + tau2B * covar.sep(X_new, d=bhat$theta[1:ncol(XF)], g=bhat$theta[ncol(XF)+1])
  y_obs <- c(yM, yF)
  
  mu_post <- t(k_star) %*% Vinv %*% y_obs
  Sigma_post <- k_new_new - t(k_star) %*% Vinv %*% k_star
  Sigma_post <- (Sigma_post + t(Sigma_post))/2
  diag(Sigma_post) <- pmax(diag(Sigma_post), eps)
  
  return(list(mean = mu_post, cov = Sigma_post, variance = diag(Sigma_post)))
}

########################### Main Loop (1 to 50) ##########################
set.seed(2025631)
setwd("~/Downloads")

for (dataset_i in 1:50) {
  cat("Running Calibration for Dataset:", dataset_i, "\n")
  
  # 1. 读取对应批次生成的数据
  ball2 <- read.csv(paste0("Calidata/2025618M_p5_noise02_", dataset_i, ".csv"))
  ball4 <- read.csv(paste0("Calidata/2025618F_p5_noise02_", dataset_i, ".csv"))
  
  ball2_0 <- ball2
  ball4_0 <- ball4
  X_oricombine <- rbind(ball2_0[,1:2], ball4_0[,1:2])
  
  ball2[, 1:2] <- scale(ball2[, 1:2])
  ball4[, 1:2] <- scale(ball4[, 1:2])
  
  XU <- as.matrix(ball2[, 1:2])
  XF <- as.matrix(ball4[, 1:2])
  
  u_bound <- 10
  u <- seq(-u_bound, u_bound, length=1000)
  
  uhats_LOOsamples_summary <- list()
  uhats_LOOsamples <- list()
  uhat_summary <- list()
  
  # 2. 估计 uhat 及 LOO 交叉验证
  yM <- as.matrix(ball2$ys)
  yF <- as.matrix(ball4$yp)
  
  yMhat <- newGPsep(XU, yM, d=0.1, g=1e-7, dK=TRUE)
  mle <- mleGPsep(yMhat, tmin=0.0001, tmax=10)
  
  formals(bhat.fit)$da <- darg(d=list(mle=TRUE), X=XF)
  formals(bhat.fit)$ga <- garg(g=list(mle=TRUE), y=yF)
  
  obj <- function(x, XF, yF, yMhat, fit) calib(x, XF, yF, yMhat, fit)$nll
  soln <- optimize(obj, lower=-u_bound, upper=u_bound, XF=XF, yF=yF, yMhat=yMhat, tol=.Machine$double.eps^0.5)
  uhat <- soln$minimum
  uhat_summary <- append(uhat_summary, uhat)
  
  uhats <- rep(NA, nrow(XF))
  for(i in 1:nrow(XF)) {
    train <- XF[-i,,drop=FALSE]
    test <- XF[i,,drop=FALSE]
    cp <- calib.pred(test, train, yF[-i], yMhat, bhat.fit)
    uhats[i] <- cp$uhat
  }
  uhats_LOOsamples_summary <- append(uhats_LOOsamples_summary, list(c(summary(uhats), sd=sd(uhats), quantile(uhats,0.025), quantile(uhats,0.975))))
  uhats_LOOsamples <- append(uhats_LOOsamples, list(uhats))
  deleteGPsep(yMhat)
  
  # 3. 抽样 100 个 u 的后验样本
  kde <- density(uhats_LOOsamples[[1]])
  unew100 <- sample(kde$x, size=100, prob=kde$y, replace=TRUE)
  
  obj_min_x <- list()
  
  # 4. 蒙特卡洛优化循环 (100 u * 100 reps = 10,000 次)
  for (ui in 1:100) {
    uhat_curr <- unew100[ui]
    
    for (rep in 1:100) {
      lhd_unit <- randomLHS(n = 500, k = 2)
      lhd_unit[,1] <- lhd_unit[,1] * (max(X_oricombine[,1]) - min(X_oricombine[,1])) + min(X_oricombine[,1])
      lhd_unit[,2] <- lhd_unit[,2] * (max(X_oricombine[,2]) - min(X_oricombine[,2])) + min(X_oricombine[,2])
      
      xnew_dm <- scale(lhd_unit)
      
      yMhat <- newGPsep(XU, yM, d=0.1, g=1e-7, dK=TRUE)
      mle <- mleGPsep(yMhat, tmin=0.0001, tmax=10)
      formals(bhat.fit)$da <- darg(d=list(mle=TRUE), X=XF)
      formals(bhat.fit)$ga <- garg(g=list(mle=TRUE), y=yF)
      
      bhat <- calib(uhat_curr, XF, yF, yMhat, bhat.fit, clean=FALSE)
      predictions <- predict_posterior(xnew_dm, XU, XF, yM, yF, uhat_curr, mle, bhat, yMhat)
      deleteGPsep(yMhat)
      
      obj <- predictions$mean^2
      min_idx <- which.min(obj)
      obj_min_x <- append(obj_min_x, list(xnew_dm[min_idx, ]))
    }
  }
  
  # 5. 反向还原坐标并存盘
  post_x <- t(matrix(unlist(obj_min_x), ncol = 10000, nrow = 2))
  X_ori <- rbind(ball2_0[,1:2], ball4_0[,1:2])
  x1_recov <- post_x[,1] * sd(X_ori[,1]) + mean(X_ori[,1])
  x2_recov <- post_x[,2] * sd(X_ori[,2]) + mean(X_ori[,2])
  
  med_x1 <- median(x1_recov)
  med_x2 <- median(x2_recov)
  sd_x1 <- sd(x1_recov)
  sd_x2 <- sd(x2_recov)
  
  # 保持文件名与 plot 脚本严格匹配
  save_path <- paste0('cali_p5_noise02_', dataset_i, '_u100rep100lhd500.RData')
  save.image(file=save_path)
}