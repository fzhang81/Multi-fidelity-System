library(knitr)
library(laGP)
library(plgp)
library(lhs)

########################### Functions ##########################

bhat.fit <- function(X, yF, Ym, da, ga, clean = TRUE) {
  bhat <- newGPsep(X, yF - Ym, d = da$start, g = ga$start, dK = TRUE)
  
  if (ga$mle) {
    cmle <- mleGPsep(bhat, param = "both", tmin = c(da$min, ga$min),
                     tmax = c(da$max, ga$max), ab = c(da$ab, ga$ab))
  } else {
    cmle <- mleGPsep(bhat, tmin = da$min, tmax = da$max, ab = da$ab)
  }
  
  cmle$nll <- -llikGPsep(bhat, dab = da$ab, gab = ga$ab)
  
  if (clean) {
    deleteGPsep(bhat)
  } else { 
    cmle$gp <- bhat
    cmle$gptype <- "sep" 
  }        
  return(cmle)
}

calib <- function(u, XF, yF, yMhat, fit, clean = TRUE) {
  XFu <- XF
  Ym <- u * predGPsep(yMhat, XFu, lite = TRUE)$mean
  cmle <- fit(XF, yF, Ym, clean = clean)
  return(cmle)
}

calib.pred <- function(XX, XF, yF, yMhat, fit) {
  soln <- optimize(obj, lower = -u_bound, upper = u_bound, XF = XF, yF = yF, 
                   yMhat = yMhat, fit = fit)
  bhat <- calib(soln$minimum, XF, yF, yMhat, fit, clean = FALSE)
  
  if (bhat$gptype == "sep") {
    pb <- predGPsep(bhat$gp, XX, lite = TRUE)
  } else {
    pb <- predGP(bhat$gp, XX, lite = TRUE)
  }
  
  pm <- predGPsep(yMhat, XX, lite = TRUE)
  m <- soln$minimum * pm$mean + pb$mean
  s2 <- soln$minimum^2 * pm$s2 + pb$s2 
  q1 <- qnorm(0.95, m, sqrt(s2))
  q2 <- qnorm(0.05, m, sqrt(s2))
  
  if (bhat$gptype == "sep") {
    deleteGPsep(bhat$gp)
  } else {
    deleteGP(bhat$gp)
  }
  return(list(mean = m, s2 = s2, q1 = q1, q2 = q2, uhat = soln$minimum))
}

ViVldet <- function(u, XU, XF, SigmaM, tau2M, mle, SigmaB) {
  # Build blocks
  SMXFu <- tau2M * covar.sep(XU, XF, d = mle$d, g = 0)
  SMu <- tau2M * covar.sep(XF, d = mle$d, g = eps)
  
  # Build V from blocks
  V <- cbind(SigmaM, u * SMXFu)
  V <- rbind(V, cbind(u * t(SMXFu), u^2 * SMu + SigmaB))
  
  # Return inverse and determinant
  return(list(V = V,
              inv = solve(V), 
              ldet = as.numeric(determinant(V, log = TRUE)$modulus)))
}

predict_posterior <- function(X_new, XU, XF, yM, yF, uhat, mle, bhat, yMhat) {
  eps <- sqrt(.Machine$double.eps)
  
  # Compute covariance matrix for the computer model
  KM <- covar.sep(XU, d = mle$d, g = eps)
  tau2M <- drop(t(yM) %*% solve(KM) %*% yM / length(yM))
  SigmaM <- tau2M * KM
  
  # Compute covariance matrix for the discrepancy function
  KB <- covar.sep(XF, d = bhat$theta[1:ncol(XF)], g = bhat$theta[ncol(XF)+1])
  Ym <- uhat * predGPsep(yMhat, XF, lite = TRUE)$mean
  YmYm <- yF - Ym
  tau2B <- drop(t(YmYm) %*% solve(KB) %*% YmYm / length(YmYm))
  SigmaB <- tau2B * KB
  
  # Get inverse of V matrix
  VV <- ViVldet(uhat, XU, XF, SigmaM, tau2M, mle, SigmaB)
  Vinv <- VV$inv
  
  # Compute cross-covariance
  k_new_M <- tau2M * covar.sep(XU, X_new, d = mle$d, g = 0)
  k_new_F <- tau2M * covar.sep(XF, X_new, d = mle$d, g = 0)
  
  # Build joint covariance vector
  k_star_block1_corrected <- uhat * k_new_M
  k_delta_XF_Xnew <- tau2B * covar.sep(XF, X_new, d = bhat$theta[1:ncol(XF)], g = bhat$theta[ncol(XF)+1])
  k_star_block2_corrected <- uhat^2 * k_new_F + k_delta_XF_Xnew
  k_star <- rbind(k_star_block1_corrected, k_star_block2_corrected)
  
  # Compute auto-covariance matrix at X_new
  k_new_new <- uhat^2 * tau2M * covar.sep(X_new, d = mle$d, g = eps) + 
    tau2B * covar.sep(X_new, d = bhat$theta[1:ncol(XF)], g = bhat$theta[ncol(XF)+1])
  
  # Combine observed data
  y_obs <- c(yM, yF)
  
  # Compute posterior mean vector
  mu_post <- t(k_star) %*% Vinv %*% y_obs
  
  # Compute full posterior covariance matrix
  Sigma_post <- k_new_new - t(k_star) %*% Vinv %*% k_star
  
  # Ensure symmetric positive definite covariance matrix
  Sigma_post <- (Sigma_post + t(Sigma_post)) / 2
  diag(Sigma_post) <- pmax(diag(Sigma_post), eps)
  
  return(list(
    mean = mu_post,
    cov = Sigma_post,
    variance = diag(Sigma_post),
    lower_95 = mu_post - 1.96 * sqrt(diag(Sigma_post)),
    upper_95 = mu_post + 1.96 * sqrt(diag(Sigma_post))
  ))
}

########################### Main Code ##########################

set.seed(20256201)
eps <- .Machine$double.eps

setwd("~/Downloads")
for (i in 1:50) {
  filename <- paste0("p5_noise02_", i)
  
  ball2 <- read.csv(paste0("Calidata/2025618M_", filename, ".csv"))
  ball4 <- read.csv(paste0("Calidata/2025618F_", filename, ".csv"))
  
  ball2_0 <- ball2
  ball4_0 <- ball4
  
  # Standardize variables using combined dataset
  X_oricombine <- rbind(ball2_0[, 1:2], ball4_0[, 1:2])
  mean_all <- apply(X_oricombine, 2, mean)
  sd_all <- apply(X_oricombine, 2, sd)
  
  ball2[, 1:2] <- scale(ball2[, 1:2], center = mean_all, scale = sd_all)
  ball4[, 1:2] <- scale(ball4[, 1:2], center = mean_all, scale = sd_all)
  
  XU <- as.matrix(ball2[, 1:2])
  XF <- as.matrix(ball4[, 1:2])
  
  u_bound <- 10
  u <- seq(-u_bound, u_bound, length = 1000) # Calibration grid
  
  uhats_LOOsamples_summary <- list()
  uhats_LOOsamples <- list()
  uhat_summary <- list()
  
  for (yi in 1:1) {
    yM <- as.matrix(ball2[, (2 + yi)])
    yF <- as.matrix(ball4[, (2 + yi)])
    
    yMhat <- newGPsep(XU, yM, d = 0.1, g = 1e-7, dK = TRUE)
    mle <- mleGPsep(yMhat, tmin = 0.0001, tmax = 10)
    
    formals(bhat.fit)$da <- darg(d = list(mle = TRUE), X = XF)
    formals(bhat.fit)$ga <- garg(g = list(mle = TRUE), y = yF)
    
    unll <- rep(NA, length(u))
    for (j in 1:length(u)) {
      unll[j] <- calib(u[j], XF, yF, yMhat, bhat.fit)$nll
    }
    
    obj <- function(x, XF, yF, yMhat, fit) calib(x, XF, yF, yMhat, fit)$nll
    soln <- optimize(obj, lower = -u_bound, upper = u_bound, XF = XF, yF = yF, 
                     yMhat = yMhat, fit = bhat.fit, tol = .Machine$double.eps^0.5)
    uhat <- soln$minimum
    uhat_summary <- append(uhat_summary, uhat)
    print(paste("y", yi, "uhat", uhat))
    
    uhats <- q1 <- q2 <- m <- s2 <- rep(NA, nrow(XF))
    for (j in 1:nrow(XF)) {
      train <- XF[-j, , drop = FALSE]
      test <- XF[j, , drop = FALSE]
      cp <- calib.pred(test, train, yF[-j], yMhat, bhat.fit)
      m[j] <- cp$mean
      s2[j] <- cp$s2 
      q1[j] <- cp$q1
      q2[j] <- cp$q2 
      uhats[j] <- cp$uhat
    }
    
    uhats_LOOsamples_summary <- append(uhats_LOOsamples_summary, 
                                       list(c(summary(uhats), sd = sd(uhats), 
                                              quantile(uhats, 0.025), quantile(uhats, 0.975))))
    uhats_LOOsamples <- append(uhats_LOOsamples, list(uhats))
    print(paste("y", yi, "uhats LOO samples"))
    print(c(summary(uhats), quantile(uhats, 0.025), quantile(uhats, 0.975)))
    
    deleteGPsep(yMhat)
  }
  
  # Generate 100u samples for each y
  unew100 <- list()
  for (yi in 1:1) {
    kde <- density(uhats_LOOsamples[[yi]])
    unew100 <- append(unew100, list(sample(kde$x, size = 100, prob = kde$y, replace = FALSE)))
  }
  
  lhd_save <- list()
  obj_save <- list()
  results <- list()
  obj_min_x <- list()
  
  for (ui in 1:100) {
    uhat <- unew100[[1]][ui]
    obj_min_indx <- vector(mode = "numeric", length = 100)
    
    for (rep in 1:100) {
      # Generate new test X
      lhd_unit <- randomLHS(n = 500, k = 2)
      
      # Map lhd_unit to the original interval
      lhd_unit[, 1] <- lhd_unit[, 1] * (max(X_oricombine[, 1]) - min(X_oricombine[, 1])) + min(X_oricombine[, 1])
      lhd_unit[, 2] <- lhd_unit[, 2] * (max(X_oricombine[, 2]) - min(X_oricombine[, 2])) + min(X_oricombine[, 2])
      
      # Standardize using combined data statistics
      lhd_unit[, 1:2] <- scale(lhd_unit[, 1:2], center = mean_all, scale = sd_all)
      
      xnew_dm <- as.matrix(lhd_unit)
      lhd_save <- append(lhd_save, list(xnew_dm))
      
      obj <- vector(mode = "numeric", length = dim(xnew_dm)[1])
      XF_new <- rbind(XF, xnew_dm)
      
      yi <- 1
      yM <- as.matrix(ball2[, (2 + yi)])
      yF <- as.matrix(ball4[, (2 + yi)])
      
      yMhat <- newGPsep(XU, yM, d = 0.1, g = 1e-7, dK = TRUE)
      mle <- mleGPsep(yMhat, tmin = 0.0001, tmax = 10)
      
      formals(bhat.fit)$da <- darg(d = list(mle = TRUE), X = XF)
      formals(bhat.fit)$ga <- garg(g = list(mle = TRUE), y = yF)
      
      bhat <- calib(uhat, XF, yF, yMhat, bhat.fit, clean = FALSE)
      
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
      
      mu_ynew_y1 <- predictions$mean
      sigma2_ynew_y1 <- predictions$cov
      
      # Required library: mvtnorm (implied by rmvnorm usage)
      obj <- rmvnorm(1, mean = mu_ynew_y1, sigma = sigma2_ynew_y1)^2
      obj_save <- append(obj_save, list(obj))
      
      obj_min_indx[rep] <- which(obj == min(obj))
      obj_min_x <- append(obj_min_x, list(xnew_dm[obj_min_indx[rep], ]))
    }
    
    results <- append(results, list(obj_min_indx))
  }
  
  post_x <- t(matrix(unlist(obj_min_x), ncol = 10000, nrow = 2))
  x1 <- post_x[, 1]
  x2 <- post_x[, 2]
  
  x1_recov <- x1 * sd_all[1] + mean_all[1]
  x2_recov <- x2 * sd_all[2] + mean_all[2]
  
  sd(x1_recov)
  sd(x2_recov)
  
  med_x1 <- median(x1_recov)
  med_x2 <- median(x2_recov)
  
  save_path <- paste0('cali_2025618F50M200_', filename, '_u100rep100lhd500.RData')
  save.image(file = save_path)
}

# Plot results
par(mfrow = c(1, 2)) 
hist(x1_recov, freq = FALSE, main = "", 
     xlab = "x1", ylab = "density", xlim = c(-1, 0),
     breaks = seq(-1.0, 0, by = 0.02))
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
