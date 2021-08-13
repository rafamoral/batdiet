################################################################################
## Kruszynski et al. - Pest suppression by bats in a human-modified landscape ##
################################################################################

## loading packages
require(hnp)
require(lme4)
require(ggplot2)
library(bivrp)
library(mvtnorm)
library(MixSIAR)
library(tidyverse)

################################################################################
## loading bat data
bats <- read.csv("bats.csv", h = T)
bats2 <- read.csv("bats_stacked.csv", h = T)
bats2$area <- bats$area
bats2$diet <- bats$diet
bats2

levels(bats$area)[1] <- "narrow"
levels(bats2$area)[1] <- "narrow"

ggplot(data = bats2, mapping = aes(x = X13C, y = X15N, colour = diet)) +
  theme_bw() + facet_wrap(~ tissue) +
  geom_point(cex = 1)

bats3 <- bats2
bats3$Area <- bats3$area
bats3$Area <- as.factor(bats3$Area)
levels(bats3$Area) <- c("Narrow","Open")
bats3$tissue <- as.factor(bats3$tissue)
levels(bats3$tissue) <- c("Fur","Liver","Wing membrane","Stomach content")
#png("fig2.tiff", w = 6, h = 4, res = 800, units = "in")
ggplot(data = subset(bats3, diet != "blood") %>%
         filter(tissue == "Fur"),
       mapping = aes(x = X13C, y = X15N, colour = Area, pch = Area)) +
  theme_bw() +
  geom_point(cex = 1.7, alpha = .85) +
  geom_point(data = subset(bats3, diet != "blood" & area == "Narrow") %>%
               filter(tissue == "fur"),
             shape = 1, cex = 1.7, colour = "black", alpha = .85) +
  geom_point(data = subset(bats3, diet != "blood" & area == "Open") %>%
               filter(tissue == "fur"),
             shape = 2, cex = 1.7, colour = "black", alpha = .85) +
  xlab(expression(paste(delta^{13}, "C (\u2030)", sep = ""))) +
  ylab(expression(paste(delta^{15}, "N (\u2030)", sep = "")))
#dev.off()

## fitting models for each tissue
fitCfur <- lmer(X13C.fur ~ area + (1 | date) + (1 | Species), data = bats)
hnp(fitCfur, verb = T, print = T)
summary(fitCfur)
fitCfur2 <- lmer(X13C.fur ~ 1 + (1 | date) + (1 | Species), data = bats)
anova(fitCfur, fitCfur2)

fitCmembrane <- lmer(X13C.membrane ~ area + (1 | date) + (1 | Species), data = subset(bats, !is.na(X13C.membrane)))
hnp(fitCmembrane, verb = T, print = T)
summary(fitCmembrane)
fitCmembrane2 <- lmer(X13C.membrane ~ 1 + (1 | date) + (1 | Species), data = subset(bats, !is.na(X13C.membrane)))
anova(fitCmembrane, fitCmembrane2)

fitCliver <- lmer(X13C.liver ~ area + (1 | date) + (1 | Species), data = subset(bats, !is.na(X13C.liver)))
hnp(fitCliver, verb = T, print = T)
summary(fitCliver)
fitCliver2 <- lmer(X13C.liver ~ 1 + (1 | date) + (1 | Species), data = subset(bats, !is.na(X13C.liver)))
anova(fitCliver, fitCliver2)

fitCstomach.content <- lmer(X13C.stomach.content ~ area + (1 | date) + (1 | Species), data = subset(bats, !is.na(X13C.stomach.content)))
hnp(fitCstomach.content, verb = T, print = T)
summary(fitCstomach.content)
fitCstomach.content2 <- lmer(X13C.stomach.content ~ 1 + (1 | date) + (1 | Species), data = subset(bats, !is.na(X13C.stomach.content)))
anova(fitCstomach.content, fitCstomach.content2)

## relationship between 13C liver and 13C fur
bats$Area <- bats$area
levels(bats$Area) <- c("Narrow","Open")
#png("fig3.tiff", w = 6, h = 4, res = 800, units = "in")
ggplot(data = bats, mapping = aes(x = X13C.fur, y = X13C.liver)) + theme_bw() +
  geom_point(aes(colour = Area)) + geom_smooth(method = "lm", col = 1) +
  xlab(expression(paste(delta^{13},C[fur]," (\u2030)",sep=""))) +
  ylab(expression(paste(delta^{13},C[liver]," (\u2030)",sep="")))
#dev.off()

fit1 <- lmer(X13C.liver ~ X13C.fur * area + (1 |date) + (1|Species), data = bats)
fit2 <- lmer(X13C.liver ~ X13C.fur + area + (1 |date) + (1|Species), data = bats)
fit3 <- lmer(X13C.liver ~ X13C.fur + (1 |date) + (1|Species), data = bats)
fit4 <- lmer(X13C.liver ~ 1 + (1 |date) + (1|Species), data = bats)
anova(fit1, fit2, fit3, fit4)

summary(fit3)

## test for shift along the axis
anova(manova(cbind(X13C.liver, X13C.fur) ~ 1, data = bats))

################################################################################
## loading insects data
ins <- read.table("insects.txt", header = TRUE)

ins_s <- split(ins[,6:7], ins$Local)
sapply(ins_s, colMeans)

ggplot(data = ins, mapping = aes(x = X15N, y = X13C, col = Local)) +
  theme_bw() +
  geom_point()

## code for fitting bivariate models
bivnormfit <- function(Y, X, covariance) {
  n <- nrow(X)
  p <- ncol(X)
  y <- cbind(Y[1:n],Y[(n+1):(2*n)])
  XtXinv <- solve(crossprod(X, X))
  beta.hat <- XtXinv %*% crossprod(X, y)
  mu.hat <- X%*%beta.hat
  sigma.hat <- 1/n * t(y - mu.hat) %*% (y - mu.hat)
  if(!covariance) sigma.hat <- diag(diag(sigma.hat))
  cov.betas <- sigma.hat %x% XtXinv
  se.s1 <- sqrt(2*sigma.hat[1]^2/(n-p+1))
  se.s2 <- sqrt(2*sigma.hat[4]^2/(n-p+1))
  if(!covariance) se.s12 <- NA else {
    rho <- sigma.hat[2]/sqrt(sigma.hat[1]*sigma.hat[4])
    se.s12 <- sqrt((1+rho^2)*sigma.hat[1]*sigma.hat[4]/(n-p+1))
  }
  se.betas <- sqrt(diag(cov.betas))
  se.sigma <- c(se.s1, se.s2, se.s12)
  coefs <- c(beta.hat, sigma.hat[1], sigma.hat[4], sigma.hat[2])
  #names(coefs) <- c("beta1.0", "beta1.1", "beta2.0", "beta2.1", "sig1", "sig2", "sig12")
  fitted <- c(mu.hat)
  resid <- Y - fitted
  Sig1 <- diag(rep(sigma.hat[1]), n)
  Sig2 <- diag(rep(sigma.hat[4]), n)
  Sig12 <- diag(rep(sigma.hat[2]), n)
  V <- rbind(cbind(Sig1, Sig12),
             cbind(Sig12, Sig2))
  llik <- dmvnorm(Y, c(mu.hat), V, log = TRUE)
  ret <- list("coefs" = coefs, "covariance" = covariance, "n" = n, 
              "X" = X, "fitted" = fitted, "resid" = resid, "loglik" = llik,
              "Y" = Y, "se" = c(se.betas, se.sigma))
  class(ret) <- "bivnormfit"
  return(ret)
}

## fitting models
X0 <- model.matrix(~ 1, data = ins)
X1 <- model.matrix(~ 0 + Local, data = ins)
Y <- c(ins[,6],ins[,7])

fit0 <- bivnormfit(Y, X0, covariance = TRUE)
fit1 <- bivnormfit(Y, X1, covariance = TRUE)
2*(fit1$loglik - fit0$loglik)
pchisq(22.45, 6, lower = FALSE)

## comparing all with pasture
X2 <- model.matrix(~ 0 + I(Local == "Pasture"), data = ins)
fit2 <- bivnormfit(Y, X2, covariance = TRUE)
2*(fit1$loglik - fit2$loglik)
pchisq(8.176994, 4, lower = FALSE) ## means that only pasture differs from every other ecosystem

## diagnostics
dfun <- function(obj) {
  r <- obj$resid
  n <- obj$n
  return(list(r[1:n], r[(n+1):(2*n)]))
}
sfun <- function(obj) {
  n <- obj$n
  fitted <- obj$fitted
  sig1 <- obj$coefs[9]
  sig2 <- obj$coefs[10]
  if(obj$covariance) sig12 <- obj$coefs[11] else sig12 <- 0
  Sig1 <- diag(rep(sig1), n)
  Sig2 <- diag(rep(sig2), n)
  Sig12 <- diag(rep(sig12), n)
  V <- rbind(cbind(Sig1, Sig12),
             cbind(Sig12, Sig2))
  mu1 <- as.numeric(obj$X %*% obj$coefs[1:4])
  mu2 <- as.numeric(obj$X %*% obj$coefs[5:8])
  Y <- as.numeric(rmvnorm(1, c(mu1, mu2), V))
  return(list(Y[1:n], Y[(n+1):(2*n)], "X" = obj$X, 
              "covariance" = obj$covariance))
}
ffun <- function(new.obj) {
  Ynew <- c(new.obj[[1]], new.obj[[2]])
  bivnormfit(Ynew, new.obj$X, new.obj$covariance)
}

bivrp(fit1, diagfun = dfun, simfun = sfun, fitfun = ffun, verb = TRUE,
      closest.angle = FALSE, density.bw = "nrd0")

################################################################################
## Bayesian mixture model
#### MixSIAR Script Bats simples ##############
## required JAGS-4.x.y.exe for R 3.3.0 or later
################################################

######  STEP 1 :Loading data ################################################
# Load the mixture/CONSUMER DATA ####################
mix1 <- load_mix_data(filename="insect_open_consumers.csv", 
                     iso_names=c("d13C","d15N"), 
                     factors=NULL, 
                     fac_random=NULL, 
                     fac_nested=NULL, 
                     cont_effects=NULL)

mix2 <- load_mix_data(filename="frugiv_consumers.csv",
                     iso_names=c("d13C","d15N"),
                     factors=NULL,
                     fac_random=NULL,
                     fac_nested=NULL,
                     cont_effects=NULL)

mix3 <- load_mix_data(filename="nectariv_consumers.csv",
                     iso_names=c("d13C","d15N"),
                     factors=NULL,
                     fac_random=NULL,
                     fac_nested=NULL,
                     cont_effects=NULL)

# Load the SOURCE DATA ####################

source1 <- load_source_data(filename="insect_sources.csv",
                           source_factors=NULL,
                           conc_dep=FALSE,
                           data_type="means",
                           mix1)

source2 <- load_source_data(filename="frugiv_sources.csv",
                           source_factors=NULL,
                           conc_dep=FALSE,
                           data_type="means",
                           mix2)

source3 <- load_source_data(filename="nectariv_sources.csv",
                           source_factors=NULL,
                           conc_dep=FALSE,
                           data_type="means",
                           mix3)

# Load DISCRIMINATION DATA ####################

discr1 <- load_discr_data(filename="insect_discrimination.csv", mix1)

discr2 <- load_discr_data(filename="frugiv_discrimination.csv", mix2)

discr3 <- load_discr_data(filename="nectariv_discrimination.csv", mix3)

#### Step 2: Exploring raw data   ############################################

#Make an isospace plot
plot_data(filename="isospace_plot1", plot_save_pdf=TRUE, plot_save_png=FALSE, mix1,source1,discr1)
plot_data(filename="isospace_plot2", plot_save_pdf=TRUE, plot_save_png=FALSE, mix2,source2,discr2)
plot_data(filename="isospace_plot3", plot_save_pdf=TRUE, plot_save_png=FALSE, mix3,source3,discr3)

# all sources must have the same probability, calculated as p = 1/n sources

#### Step 3: Write JAGS model file ###############################################

# Define model structure and write JAGS model file
model_filename <- "MixSIAR_model_kw_uninf.txt"   # Name of the JAGS model file
resid_err <- TRUE
process_err <- TRUE

#### Step 4: Run model ##########################################################
# Run the JAGS model ("very long" took ~5 min)
write_JAGS_model(model_filename, resid_err, process_err, mix1, source1)
#jags.uninf1 <- run_model(run="test",mix1,source1,discr1,model_filename,alpha.prior = 1, resid_err, process_err)
jags.uninf1 <- run_model(run="very long",mix1,source1,discr1,model_filename,alpha.prior = 1, resid_err, process_err)

write_JAGS_model(model_filename, resid_err, process_err, mix2, source2)
#jags.uninf2 <- run_model(run="test",mix2,source2,discr2,model_filename,alpha.prior = 1, resid_err, process_err)
jags.uninf2 <- run_model(run="very long",mix2,source2,discr2,model_filename,alpha.prior = 1, resid_err, process_err)

write_JAGS_model(model_filename, resid_err, process_err, mix3, source3)
#jags.uninf3 <- run_model(run="test",mix3,source3,discr3,model_filename,alpha.prior = 1, resid_err, process_err)
jags.uninf3 <- run_model(run="very long",mix3,source3,discr3,model_filename,alpha.prior = 1, resid_err, process_err)

# Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.uninf1, mix1, source1)
output_JAGS(jags.uninf2, mix2, source2)
output_JAGS(jags.uninf3, mix3, source3)

# analyse output settings
output_options <- list(summary_save = TRUE,
                       summary_name = "summary_statistics",
                       sup_post = FALSE,
                       plot_post_save_pdf = TRUE,
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,
                       plot_pairs_save_pdf = TRUE,
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,
                       geweke = TRUE,
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,
                       plot_post_save_png = FALSE,
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE)


#   Process diagnostics, summary stats, and posterior plots
output_JAGS(jags.uninf1, mix1, source1, output_options)
output_JAGS(jags.uninf2, mix2, source2, output_options)
output_JAGS(jags.uninf3, mix3, source3, output_options)