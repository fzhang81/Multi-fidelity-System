##################################################
# Main Code
##################################################
set.seed(20250734)
eps <- .Machine$double.eps

##################################################
# Functions
##################################################

bhat.fit <- function(X, yF, Ym, da, ga, clean=TRUE) {
  bhat <- newGPsep(X, yF - Ym, d=da$start, g=ga$start, dK=TRUE)
  
  if(ga$mle) {
    cmle <- mleGPsep(bhat, param="both", tmin=c(da$min, ga$min),
                     tmax=c(da$max, ga$max), ab=c(da$ab, ga$ab))
  } else {
    cmle <- mleGPsep(bhat, tmin=da$min, tmax=da$max, ab=da$ab)
  }
  cmle$nll <- - llikGPsep(bhat, dab=da$ab, gab=ga$ab)
  
  if(clean) {
    deleteGPsep(bhat)
  } else { 
    cmle$gp <- bhat; cmle$gptype <- "sep" 
  }        
  
  return(cmle)
}

calib <- function(u, XF, yF, yMhat, fit, clean=TRUE) {
  XFu <- XF
  Ym <- u * predGPsep(yMhat, XFu, lite=TRUE)$mean
  cmle <- fit(XF, yF, Ym, clean=clean)
  return(cmle)
}

calib.pred <- function(XX, XF, yF, yMhat, fit) {
  soln <- optimize(obj, lower=-u_bound, upper=u_bound, XF=XF, yF=yF, 
                   yMhat=yMhat, fit=fit)
  bhat <- calib(soln$minimum, XF, yF, yMhat, fit, clean=FALSE)
  
  if(bhat$gptype == "sep") {
    pb <- predGPsep(bhat$gp, XX, lite=TRUE)
  } else {
    pb <- predGP(bhat$gp, XX, lite=TRUE)
  }
  pm <- predGPsep(yMhat, XX, lite=TRUE)
  m <- soln$minimum * pm$mean + pb$mean
  s2 <- soln$minimum^2 * pm$s2 + pb$s2 
  q1 <- qnorm(0.95, m, sqrt(s2))
  q2 <- qnorm(0.05, m, sqrt(s2))
  
  if(bhat$gptype == "sep") deleteGPsep(bhat$gp) else deleteGP(bhat$gp)
  
  return(list(mean=m, s2=s2, q1=q1, q2=q2, uhat=soln$minimum))
}



ViVldet <- function(u, XU, XF, SigmaM, tau2M, mle, SigmaB) {
  
  ## Build blocks
  SMXFu <- tau2M * covar.sep(XU, XF, d=mle$d, g=0)
  SMu <- tau2M * covar.sep(XF, d=mle$d, g=eps)
  
  ## Build V from blocks
  V <- cbind(SigmaM, u*SMXFu)
  V <- rbind(V, cbind(u*t(SMXFu), u^2*SMu + SigmaB))

  ## Return inverse and determinant
  return(list(V=V,
              inv=solve(V), 
              ldet=as.numeric(determinant(V, log=TRUE)$modulus)))
}

predict_posterior <- function(X_new, XU, XF, yM, yF, uhat, mle, bhat, yMhat) {
  eps <- sqrt(.Machine$double.eps)
  
  # Calculate covariance matrices
  KM <- covar.sep(XU, d=mle$d, g=eps)
  tau2M <- drop(t(yM) %*% solve(KM) %*% yM / length(yM))
  SigmaM <- tau2M * KM

  KB <- covar.sep(XF, d=bhat$theta[1:ncol(XF)], g=bhat$theta[ncol(XF)+1])
  Ym <- uhat * predGPsep(yMhat, XF, lite=TRUE)$mean
  YmYm <- yF - Ym
  tau2B <- drop(t(YmYm) %*% solve(KB) %*% YmYm / length(YmYm))
  SigmaB <- tau2B * KB
  
  # Invert V matrix
  VV <- ViVldet(uhat, XU, XF, SigmaM, tau2M, mle, SigmaB)
  Vinv <- VV$inv
 
  # Calculate cross-covariance
  k_new_M <- tau2M * covar.sep(XU, X_new, d=mle$d, g=0)
  k_new_F <- tau2M * covar.sep(XF, X_new, d=mle$d, g=0)
  
  # Build joint covariance vector
  k_star_block1_corrected <- uhat * k_new_M
  k_delta_XF_Xnew <- tau2B * covar.sep(XF, X_new, d=bhat$theta[1:ncol(XF)], g=bhat$theta[ncol(XF)+1])
  k_star_block2_corrected <- uhat^2 * k_new_F + k_delta_XF_Xnew
  k_star <- rbind(k_star_block1_corrected, k_star_block2_corrected)
  
  # Calculate auto-covariance matrix at X_new
  k_new_new <- uhat^2 * tau2M * covar.sep(X_new, d=mle$d, g=eps) + 
    tau2B * covar.sep(X_new, d=bhat$theta[1:ncol(XF)], g=bhat$theta[ncol(XF)+1])

  # Combine observation data
  y_obs <- c(yM, yF)

  # Posterior mean and covariance
  mu_post <- t(k_star) %*% Vinv %*% y_obs
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

# Update with nested data filenames
filenameF <- "Curing_33uF_noise01_lhd50_nested"  # High fidelity (50 points, nested)
filenameM <- "Curing_44DoCS_noise0_lhd200"      # Low fidelity (200 points)

ball4 <- read.csv(paste0(filenameF, ".csv"))
ball2 <- read.csv(paste0(filenameM, ".csv"))

ball2_0 <- ball2
ball4_0 <- ball4

# Combine data for standardization
X_oricombine <- rbind(ball2_0[, 1:4], ball4_0[, 1:4])
mean_all <- apply(X_oricombine, 2, mean)
sd_all <- apply(X_oricombine, 2, sd)

ball2[, 1:4] <- scale(ball2[, 1:4], center = mean_all, scale = sd_all)
ball4[, 1:4] <- scale(ball4[, 1:4], center = mean_all, scale = sd_all)

XU <- as.matrix(ball2[, 1:4])
XF <- as.matrix(ball4[, 1:4])

u_bound <- 100
u <- seq(-u_bound, u_bound, length=1000)

uhats_LOOsamples_summary <- list()
uhats_LOOsamples <- list()
uhat_summary <- list()

# Extract the 5th column for y observations
for (yi in 1:1) {
  yM <- as.matrix(ball2[, (4+yi)])
  yF <- as.matrix(ball4[, (4+yi)])
  
  yMhat <- newGPsep(XU, yM, d=0.1, g=1e-7, dK=TRUE)
  mle <- mleGPsep(yMhat, tmin=0.0001, tmax=10)
  
  formals(bhat.fit)$da <- darg(d=list(mle=TRUE), X=XF)
  formals(bhat.fit)$ga <- garg(g=list(mle=TRUE), y=yF)
  
  unll <- rep(NA, length(u))
  for(i in 1:length(u)) {
    unll[i] <- calib(u[i], XF, yF, yMhat, bhat.fit)$nll
  }

  obj <- function(x, XF, yF, yMhat, fit) calib(x, XF, yF, yMhat, fit)$nll
  soln <- optimize(obj, lower=-u_bound, upper=u_bound, XF=XF, yF=yF, 
                   yMhat=yMhat, fit=bhat.fit, tol = .Machine$double.eps^0.5)
  
  uhat <- soln$minimum
  uhat_summary <- append(uhat_summary, uhat)
  print(paste("y", yi, "uhat", uhat))
 
  uhats <- q1 <- q2 <- m <- s2 <- rep(NA, nrow(XF))
  for(i in 1:nrow(XF)) {
    train <- XF[-i, , drop=FALSE]
    test <- XF[i, , drop=FALSE]
    cp <- calib.pred(test, train, yF[-i], yMhat, bhat.fit)
    m[i] <- cp$mean
    s2[i] <- cp$s2 
    q1[i] <- cp$q1
    q2[i] <- cp$q2 
    uhats[i] <- cp$uhat
  }

  uhats_LOOsamples_summary <- append(uhats_LOOsamples_summary, 
                                     list(c(summary(uhats), sd=sd(uhats), quantile(uhats,0.025), quantile(uhats,0.975))))
  uhats_LOOsamples <- append(uhats_LOOsamples, list(uhats))
  print(paste("y", yi, "uhats LOO samples"))
  print(c(summary(uhats), quantile(uhats,0.025), quantile(uhats,0.975)))
 
  deleteGPsep(yMhat)
}

# Generate 100u samples
unew100 <- list()
for (yi in 1:1) {
  kde <- density(uhats_LOOsamples[[yi]])
  unew100 <- append(unew100, list(sample(kde$x, size=100, prob=kde$y, replace=FALSE)))
}

lhd_save <- list()
obj_save <- list()
results <- list()
obj_min_x <- list()

# Optimization over Latin Hypercube samples
for (ui in 1:100) {
  uhat <- unew100[[1]][ui]
  obj_min_indx <- vector(mode="numeric", length=100)

  for (rep in 1:100) {
    lhd_unit <- randomLHS(n = 200, k = 4)

    lhd_unit[, 1] <- lhd_unit[, 1] * (max(X_oricombine[, 1]) - min(X_oricombine[, 1])) + min(X_oricombine[, 1])
    lhd_unit[, 2] <- lhd_unit[, 2] * (max(X_oricombine[, 2]) - min(X_oricombine[, 2])) + min(X_oricombine[, 2])
    lhd_unit[, 3] <- lhd_unit[, 3] * (max(X_oricombine[, 3]) - min(X_oricombine[, 3])) + min(X_oricombine[, 3])
    lhd_unit[, 4] <- lhd_unit[, 4] * (max(X_oricombine[, 4]) - min(X_oricombine[, 4])) + min(X_oricombine[, 4])
 
    lhd_unit[, 1:4] <- scale(lhd_unit[, 1:4], center = mean_all, scale = sd_all)
    
    xnew_dm <- as.matrix(lhd_unit)
    lhd_save <- append(lhd_save, list(xnew_dm))
    
    XF_new <- rbind(XF, xnew_dm)

    yi <- 1
    yM <- as.matrix(ball2[, (4+yi)])
    yF <- as.matrix(ball4[, (4+yi)])
    
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

    mu_ynew_y1 <- predictions$mean
    sigma2_ynew_y1 <- predictions$cov

    obj <- rmvnorm(1, mean=mu_ynew_y1, sigma=sigma2_ynew_y1)^2
    obj_save <- append(obj_save, list(obj))

    obj_min_indx[rep] <- which(obj == min(obj))
    obj_min_x <- append(obj_min_x, list(xnew_dm[obj_min_indx[rep], ]))
  }
 
  results <- append(results, list(obj_min_indx))
}

post_x <- t(matrix(unlist(obj_min_x), ncol = 10000, nrow = 4))
x1 <- post_x[, 1]
x2 <- post_x[, 2]
x3 <- post_x[, 3]
x4 <- post_x[, 4]

x1_recov <- x1 * sd_all[1] + mean_all[1]
x2_recov <- x2 * sd_all[2] + mean_all[2]
x3_recov <- x3 * sd_all[3] + mean_all[3]
x4_recov <- x4 * sd_all[4] + mean_all[4]

med_x1 <- median(x1_recov)
med_x2 <- median(x2_recov)
med_x3 <- median(x3_recov)
med_x4 <- median(x4_recov)

# Save workspace
save_path <- paste0('Results_NestedData_Calibration.RData')
save.image(file=save_path)

##################################################
# Plot Results
##################################################

par(mfrow=c(2,2), mar=c(4.5, 4.5, 1, 1)) 

# Top-left plot
hist(x1_recov, freq=FALSE, main="", 
     xlab=expression(T[1]~"("*degree*C*")"), ylab="density", xlim=c(125, 180),
     breaks=seq(125, 180, by=2))
lines(density(x1_recov), col="red", lwd=1)
abline(v=med_x1, col="red", lwd=2, lty=2) 
text(x = med_x1 + 0.01, y = max(density(x1_recov)$y)*0.9, labels = round(med_x1, 2), col = "red", pos = 4)

# Top-right plot (ylab removed)
hist(x2_recov, freq=FALSE, main="",
     xlab=expression(T[2]~"("*degree*C*")"), ylab="", xlim=c(150, 180),
     breaks=seq(150, 180, by=1))
lines(density(x2_recov), col="red", lwd=1)
abline(v=med_x2, col="red", lwd=2, lty=2) 
text(x = med_x2 + 0.01, y = max(density(x2_recov)$y) * 0.9, labels = round(med_x2, 2), col = "red", pos = 4)

# Bottom-left plot
hist(x3_recov, freq=FALSE, main="",
     xlab=expression(t[1]~"(min)"), ylab="density", xlim=c(600, 6600),
     breaks=seq(600, 6600, by=200))
lines(density(x3_recov), col="red", lwd=1)
abline(v=med_x3, col="red", lwd=2, lty=2) 
text(x = med_x3 + 200, y = max(density(x3_recov)$y) * 0.9, labels = round(med_x3, 2), col = "red", pos = 4)

# Bottom-right plot (ylab removed)
hist(x4_recov, freq=FALSE, main="",
     xlab=expression(t[2]~"(min)"), ylab="", xlim=c(7240, 12000),
     breaks=seq(7200, 12000, by=100))
lines(density(x4_recov), col="red", lwd=1)
abline(v=med_x4, col="red", lwd=2, lty=2) 
text(x = med_x4 - 1100, y = max(density(x4_recov)$y) * 0.9, labels = round(med_x4, 2), col = "red", pos = 4)

# Restore default par settings to avoid impacting future plots
par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
