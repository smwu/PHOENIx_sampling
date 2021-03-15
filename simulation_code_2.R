#================================================#
# Simulation Code for evaluating sampling designs 
# for a nested case-control study within a CRT
#
# Updated 3/15/2021
#================================================#


library(SimCorMultRes)
library(evd)
library(dplyr)
library(geepack)
library(parallel)
library(foreach)
library(doParallel)
library(ggplot2)
library(tidyr)
library(gridExtra)

# https://cran.r-project.org/web/packages/SimCorMultRes/SimCorMultRes.pdf 
# https://cran.r-project.org/web/packages/SimCorMultRes/vignettes/SimCorMultRes.html
# https://journal.r-project.org/archive/2016/RJ-2016-034/index.html

#============================= Simulations ====================================

simulate_binary_data <- function(K, nk, p_baseline, OR, p_exposed, gamma) {
  # Function to simulate data with binary outcome and binary covariate
  # K: number of clusters
  # nk: size per cluster
  # p_baseline: P(Y=1|X=0)
  # OR: odds ratio
  # p_exposed: P(X=1)
  # gamma: correlation between the latent variables
  
  # sample size
  sample_size <- K
  # cluster size
  cluster_size <- nk
  # intercept
  beta_intercepts <- log(p_baseline/(1-p_baseline))
  # regression parameter associated with the covariate
  beta_coefficients <- log(OR)
  # individual-level covariate with proportions that may differ acros clusters
  x <- rbinom(cluster_size*sample_size, 1, prob=p_exposed)
  
  library(bindata)
  for (i in 1:sample_size) {
    rmvbin(cluster_size[i], margprob=p_exposed, bincorr)
  }
  x <- rmvbin()
  
  # # correlation matrix for the NORTA method
  # latent_correlation_matrix <- toeplitz(c(1, gamma, gamma, gamma))
  # # simulation of clustered binary responses
  # simulated_binary_dataset <- rbin(clsize = cluster_size, intercepts = beta_intercepts, 
  #                                  betas = beta_coefficients, xformula = ~x, cor.matrix = latent_correlation_matrix, 
  #                                  link = "probit")
  
  # simulation of correlated latent continuous variables
  simulated_latent_variables1 <- rmvevd(sample_size, dep = sqrt(1 - gamma), model = "log", 
                                        d = cluster_size)
  simulated_latent_variables2 <- rmvevd(sample_size, dep = sqrt(1 - gamma), model = "log", 
                                        d = cluster_size)
  simulated_latent_variables <- simulated_latent_variables1 - simulated_latent_variables2
  
  # simulation of clustered binary responses
  simulated_binary_dataset <- rbin(clsize = cluster_size, intercepts = beta_intercepts, 
                                   betas = beta_coefficients, xformula = ~x, rlatent = simulated_latent_variables)
  sim_data <- simulated_binary_dataset$simdata
  
  colnames(sim_data) <- c("y", "x", "cluster", "id")
  # order by cluster so that geeglm recognizes the strucutre
  sim_data <- group_by(sim_data, cluster)
  
  # write.csv(simulated_binary_dataset, file = "simdata_bin.csv", row.names = FALSE)
  return(sim_data)
}


gen_sub_des1 <- function(data, ratio) {
  # Function to generate subcohort for design 1: stratified one-stage cluster sampling
  # Returns subsetted data with weights as a new column
  # data: full dataset to sample from
  # ratio: control-to-case ratio
  
  # cluster status = 1 if cluster contains at least one case; 0 otherwise
  cluster_status <- data %>% group_by(cluster) %>% summarise(max = max(y), .groups = "keep")
  data <- merge(data, cluster_status, by="cluster")
  
  # get case clusters and label with type = 1
  case_clusters <- filter(data, max==1) %>% mutate(type = 1)
  n_cases <- sum(case_clusters$y==1)  # number of cases
  controls_in_case <- sum(case_clusters$y==0)  # number of controls in case clusters
  
  # get control clusters and label with type = 0
  control_clusters <- filter(data, max==0) %>% mutate(type = 0)
  # SRS of control cluster indices
  srs <- sample(x=unique(control_clusters$cluster), size=max(0, (ratio*n_cases-controls_in_case)/2))
  sub_controls <- filter(control_clusters, cluster %in% srs)  # get selected control clusters
  # combine data from case clusters and selected control clusters
  sub_cohort <- bind_rows(case_clusters, sub_controls) %>% group_by(cluster)
  
  w1 <- 1 
  w0 <- length(unique(control_clusters$cluster)) / length(srs)  # C0/c0
  
  sub_cohort <- mutate(sub_cohort, weights = ifelse(type==1, w1, w0))
  output <- sub_cohort %>% group_by(cluster)
  
  # fit IPWGEE and get robust variance for beta1
  fit <- geeglm(y ~ x, id=cluster, data=output, weights = weights, 
                family=binomial, scale.fix=TRUE, corstr="independence")
  var_slope <- summary(fit)$coefficients[2,2]
  
  return(list("var_slope" = var_slope, "estimate" = summary(fit)$coefficients[2,1], "output" = as.data.frame(output)))
}

gen_sub_des2 <- function(data, ratio) {
  # Function to generate subcohort for design 2: stratified two-stage cluster sampling
  # Returns subsetted data with weights as a new column
  # data: full dataset to sample from
  # ratio: control-to-case ratio
  
  # cluster status = 1 if cluster contains at least one case; 0 otherwise
  cluster_status <- data %>% group_by(cluster) %>% summarise(max = max(y), .groups = "keep")
  data <- merge(data, cluster_status, by="cluster")
  
  # get case clusters and label with type = 1
  case_clusters <- filter(data, max==1) %>% mutate(type = 1)
  n_cases <- sum(case_clusters$y==1)  # number of cases
  controls_in_case <- sum(case_clusters$y==0)  # number of controls in case clusters
  
  # get control clusters and label with type = 0
  control_clusters <- filter(data, max==0) %>% mutate(type = 0)
  # SRS of control cluster indices
  srs <- sample(x=unique(control_clusters$cluster), size=max(0, ratio*n_cases-controls_in_case))
  sub_controls <- filter(control_clusters, cluster %in% srs)  # get selected control clusters
  # combine data from case clusters and selected control clusters
  sub_cohort <- bind_rows(case_clusters, sub_controls) %>% group_by(cluster)
  
  w1 <- 1 
  w0 <- 2*length(unique(control_clusters$cluster)) / length(srs)  # 2C0/c0
  
  sub_cohort <- mutate(sub_cohort, weights = ifelse(type==1, w1, w0))
  output <- sub_cohort %>% group_by(cluster)
  
  # fit IPWGEE and get robust variance for beta1
  fit <- geeglm(y ~ x, id=cluster, data=output, weights = weights, 
                family=binomial, scale.fix=TRUE, corstr="independence")
  var_slope <- summary(fit)$coefficients[2,2]
  
  return(list("var_slope" = var_slope, "estimate" = summary(fit)$coefficients[2,1], "output" = as.data.frame(output)))
}


gen_sub_des3 <- function(data, ratio) {
  # Function to generate subcohort for design 3: dispersed stratified two-stage cluster sampling
  # Returns subsetted data with weights as a new column
  # data: full dataset to sample from
  # ratio: control-to-case ratio
  
  # cluster status = 1 if cluster contains at least one case; 0 otherwise
  cluster_status <- data %>% group_by(cluster) %>% summarise(max = max(y), .groups = "keep")
  data <- merge(data, cluster_status, by="cluster")
  
  # get cases
  cases <- filter(data, y==1) %>% mutate(type = 1)
  
  # get control clusters and label with type = 0
  control_clusters <- filter(data, max==0) %>% mutate(type = 0)
  # SRS of control cluster indices
  srs <- sample(x=unique(control_clusters$cluster), size=max(0, ratio*nrow(cases)))
  sub_controls <- filter(control_clusters, cluster %in% srs)  # get selected control clusters
  # combine data from cases and selected control clusters
  sub_cohort <- bind_rows(cases, sub_controls) %>% group_by(cluster)
  
  w1 <- 1 
  w0 <- 2*length(unique(control_clusters$cluster)) / length(srs)  # 2C0/c0
  
  sub_cohort <- mutate(sub_cohort, weights = ifelse(type==1, w1, w0))
  output <- sub_cohort %>% group_by(cluster)
  
  # fit IPWGEE and get robust variance for beta1
  fit <- geeglm(y ~ x, id=cluster, data=output, weights = weights, 
                family=binomial, scale.fix=TRUE, corstr="independence")
  var_slope <- summary(fit)$coefficients[2,2]
  
  return(list("var_slope" = var_slope, "estimate" = summary(fit)$coefficients[2,1], "output" = as.data.frame(output)))
}

gen_sub_des4 <- function(data, ratio) {
  # Function to generate subcohort for design 4: standard case-control
  # Returns subsetted data with weights as a new column
  # data: full dataset to sample from
  # ratio: control-to-case ratio
  cases <- filter(data, y==1)
  controls <- filter(data, y==0)
  srs <- sample(x=1:nrow(controls), size=min(nrow(controls), ratio*nrow(cases)))
  sub_controls <- controls[srs,]  # select sub-cohort of controls
  w1 <- 1  #N1/n1
  w0 <- nrow(controls)/nrow(sub_controls)  # N0/n0
  output <- rbind(cases, sub_controls) %>% mutate(weights = ifelse(y==1, 1, w0)) %>% group_by(cluster)
  
  # fit IPWGEE and get robust variance for beta1
  fit <- geeglm(y ~ x, id=cluster, data=output, weights = weights, 
                family=binomial, scale.fix=TRUE, corstr="independence")
  var_slope <- summary(fit)$coefficients[2,2]
  
  return(list("var_slope" = var_slope, "estimate" = summary(fit)$coefficients[2,1], "output" = as.data.frame(output)))
}

gen_sub_des5 <- function(data, ratio) {
  # Function to generate subcohort for Design 5: case-control with clustering
  # Returns subsetted data with weights as a new column
  # data: full dataset to sample from
  # ratio: control-to-case ratio
  cases <- filter(data, y==1)
  srs <- sample(x=unique(data$cluster), size=ratio*nrow(cases))  # SRS of cluster indices
  sub_controls <- filter(data, cluster %in% srs) %>%   # get selected clusters
    filter(y==0) %>%
    group_by(cluster) %>%
    group_modify(~ sample_n(.x, size=1))  # for each selected cluster, choose 1 control
  sub_cohort <- bind_rows(cases, sub_controls) %>% group_by(cluster)
  
  w1 <- 1  # N1/n1
  w10 <- 2*length(unique(data$cluster)) / length(srs)  # 2C/c
  w0 <- w10/2  # C/c
  
  type <- tapply(sub_cohort$y, sub_cohort$cluster, function(y) {ifelse(length(y)>1, 2, ifelse(y==1, 1, 0))})
  unique_clusters <- unique(sub_cohort$cluster)
  temp <- data.frame(cbind(unique_clusters, type))
  colnames(temp) <- c("cluster", "type")
  sub_cohort <- merge(sub_cohort, temp, by="cluster") %>%
    mutate(weights = ifelse(type==2, w10, ifelse(type==1, w1, w0)))
  
  output <- sub_cohort %>% group_by(cluster)
  
  # fit IPWGEE and get robust variance for beta1
  fit <- geeglm(y ~ x, id=cluster, data=output, weights = weights, 
                family=binomial, scale.fix=TRUE, corstr="independence")
  var_slope <- summary(fit)$coefficients[2,2]
  
  return(list("var_slope" = var_slope, "estimate" = summary(fit)$coefficients[2,1], "output" = as.data.frame(output)))
}



temp_trials <- function(gamma, MC_sims) {
  # Function to obtain simulation estimates and correlation in the full cohort (no weighting required)
  # gamma: pairwise correlation between latent e_ki^B random variables
  # MC_sims: number of monte carlo simulations
  temp <- as.data.frame(matrix(nrow=MC_sims, ncol=3))
  for (i in 1:MC_sims) {
    # simulate data
    data <- simulate_binary_data(K=734, nk=2, p_baseline=0.025, OR=2.05, p_exposed=0.2, gamma=gamma)
    fit <- summary(geeglm(y ~ x, id=cluster, data=data, family=binomial, scale.fix=TRUE, corstr="exchangeable"))
    temp[i, ] <- c(fit$corr[1,1], fit$coefficients[1,1], fit$coefficients[2,1])
  }
  apply(temp, 2, mean)
}

# Examine the relationship between gamma correlation and alpha correlation
gamma <- seq(0, 0.6, by=0.02)
temp_results <- as.data.frame(matrix(nrow=length(gamma), ncol=3))
colnames(temp_results) <- c("alpha", "beta0", "beta1")

for (i in 1:length(gamma)) {
  temp_results[i,] <- temp_trials(gamma[i], 100)
}

print(temp_results)
## When gamma = 0, alpha = 0.0023, P(Y=1|X=0) = 0.025, P(Y=1|X=1) = 1.97
## When gamma = 0.3 [16], alpha = 0.1496, P(Y=1|X=0) = 0.025, P(Y=1|X=1) = 2
## When gamma = 0.5 [26], alpha = 0.2971, P(Y=1|X=0) = 0.0254, P(Y=1|X=1) = 2.5



#================ Examine simulated data ===========================================================

library(dplyr)
library(ICCbin)
data <- simulate_binary_data(K=734, nk=2, p_baseline=0.025, OR=2.05, p_exposed=0.2, gamma=0.3)

get_data_chars <- function(data) {
  num_per_cluster <- data %>% count(cluster)
  nk_table <- table(num_per_cluster$n)
  nk_1 <- sum(num_per_cluster$n == 1)
  nk_2 <- sum(num_per_cluster$n == 2)
  mean_nk <- mean(num_per_cluster$n)
  median_nk <- median(num_per_cluster$n)
  xy_table <- addmargins(xtabs(~ x + y, data=data))
  prop_table <- xy_table / xy_table[3,3]
  p_case <- prop_table[3,2]
  p_ctrl <- prop_table[3,1]
  p_case_cond_exp <- prop_table[2,2]/prop_table[2,3]
  p_case_cond_unexp <- prop_table[1,2]/prop_table[1,3]
  OR <- (prop_table[2,2]/prop_table[2,1])/(prop_table[1,2]/prop_table[1,1])
  # paste0("P(Y=1): ", p_case, ", P(Y=0): ", p_ctrl, ", P(Y=1|X=1): ", 
  #        p_case_cond_exp, ", P(Y=1|X=0): ", p_case_cond_unexp, ", OR: ", OR)
  case_distr <- data %>% group_by(cluster) %>% summarise(clust_sum = sum(y), .groups='drop')
  K_case_table <- table(case_distr$clust_sum)
  # Use moment estimate from unbiased estimating equation
  # Potentially also use: Fleiss-Cuzick or Resampling method estimate
  # ICC <- iccbin(cid = cluster, y = y, data = data, method = "ub", ci.type="fc")
  temp <- data %>% group_by(cluster) %>% summarise(n = n(), y = sum(y), .groups='drop')
  ICC <- aod::iccbin(n, y, data=temp, method="C")@rho
  
  # P(X=1) = prop_table[2,3]
  sim_chars <- c(prop_table[2,3], p_case, p_case_cond_exp, p_case_cond_unexp, OR, ICC,  
                 K_case_table[1], K_case_table[2], sum(K_case_table[3:length(K_case_table)]),
                 nk_1, nk_2, mean_nk, median_nk)
  
  return(sim_chars)
}


fit_full <- geeglm(y ~ x, id=cluster, data=data, family=binomial, scale.fix=TRUE, corstr="exchangeable")
summary(fit_full)$corr

# chars_df <- as.data.frame(matrix(NA, nrow=13, ncol=5))

# Sanity check design characteristics after M simulations 
MC <- 100

# Initialize dataframes for the various designs
df_full <- as.data.frame(matrix(NA, nrow=13, ncol=MC))
df_Des1 <- as.data.frame(matrix(NA, nrow=13, ncol=MC))
df_Des2 <- as.data.frame(matrix(NA, nrow=13, ncol=MC))
df_Des3 <- as.data.frame(matrix(NA, nrow=13, ncol=MC))
df_Des4 <- as.data.frame(matrix(NA, nrow=13, ncol=MC))
df_Des5 <- as.data.frame(matrix(NA, nrow=13, ncol=MC))

# Set simulation parameters based on PHOENIx trial
K <- 734; nk <- 2; p_baseline <- 0.025; OR <- 2.05; p_exposed <- 0.2; gamma <- 0.3

# Run simulations
for (i in 1:MC) {
  data <- simulate_binary_data(K=K, nk=nk, p_baseline=p_baseline, OR=OR, p_exposed=p_exposed, gamma=gamma)
  df_full[,i] <- get_data_chars(data)
}

chars_df <- as.data.frame(cbind(apply(df_full, 1, mean))) 
colnames(chars_df) <- c("Full", "Des 1", "Des 2", "Des 3", "Des 4", "Des 5")
rownames(chars_df) <- c("P(X=1)", "P(Y=1)", "P(Y=1|X=1)", "P(Y=1|X=0)", "OR", "ICC", 
                        "K_case_0", "K_case_1", "K_case>2", "nk=1", "nk=2", "mean nk", "median nk")



library(mvtnorm)
x_it <- -3
cov <- pmvnorm(upper=c(0.2*x_it, 0.2*x_it), 
               corr=matrix(c(1,0.9, 0.9, 1), nrow=2, byrow=TRUE))[1] - (pnorm(0.2*x_it))^2
cor <- cov/(0.2*x_it*(1-0.2*x_it))
  
# Running the code: simulate_binary_data(K=734, nk=2, p_baseline=0.025, OR=2.05, p_exposed=0.2, gamma=0.3),
# We have P(Y=1)=0.033, P(Y=1|X=1)=0.04, P(Y=1|X=0)=0.03, and OR=1.39. The OR should be higher, but this may 
# be because the dichotomized model attenuates the OR



### Run the simulations on a local cluster

cl <- makeCluster(4)
registerDoParallel(cl)

# Results for alpha = 0.15
MC_sims <- 1000
ptm <- proc.time()
oper <- foreach (icount(MC_sims), .combine=rbind, .packages=c("SimCorMultRes", "evd", "dplyr", "geepack")) %dopar% {
  
  # Simulate data
  data <- simulate_binary_data(K=734, nk=2, p_baseline=0.025, OR=2.05, p_exposed=0.2, gamma=0.3)
  
  # Get variance estimates
  des5 <- gen_sub_des5(data, ratio=2)  # 80 clusters, 22 exposed, weights 17.1, 34.1
  des1 <- gen_sub_des1(data, ratio=2)  # 64 clusters, 35 exposed, weights 30.1
  des2 <- gen_sub_des2(data, ratio=2)  # 88 clusters, 44 exposed, weights 29.5
  des3 <- gen_sub_des3(data, ratio=2)  # 127 clusters, 42 exposed, weights 16.1
  des4 <- gen_sub_des4(data, ratio=2)  # 123 clusters, 37 exposed, weights 16.6
  
  # Get sampling design characteristics
  sampling_chars <- list(des5$output, des1$output, des2$output, des3$output, des4$output)
  n_clus <- sapply(sampling_chars, function(x) length(unique(x$cluster)))
  max_wt <- sapply(sampling_chars, function(x) max(x$weights))
  n_cases <- sapply(sampling_chars, function(x) sum(x$y))
  n_exp <- sapply(sampling_chars, function(x) sum(x$x))
  
  result <- c(des5[[1]], des1[[1]], des2[[1]], des3[[1]], des4[[1]], 
              des5[[2]], des1[[2]], des2[[2]], des3[[2]], des4[[2]], n_clus, max_wt, n_cases, n_exp)
  return(result)
}
colnames(oper) <- rep(paste0("Design ", 0:4), 6)
results <- apply(oper, 2, mean)
proc.time() - ptm  

rm(oper)

### Results for alpha = 0
ptm <- proc.time()
oper <- foreach (icount(MC_sims), .combine=rbind, .packages=c("SimCorMultRes", "evd", "dplyr", "geepack")) %dopar% {
  
  # Simulate data
  data <- simulate_binary_data(K=734, nk=2, p_baseline=0.025, OR=2.05, p_exposed=0.2, gamma=0)
  
  # Get variance estimates
  des5 <- gen_sub_des5(data, ratio=2)  # 80 clusters, 22 exposed, weights 17.1, 34.1
  des1 <- gen_sub_des1(data, ratio=2)  # 64 clusters, 35 exposed, weights 30.1
  des2 <- gen_sub_des2(data, ratio=2)  # 88 clusters, 44 exposed, weights 29.5
  des3 <- gen_sub_des3(data, ratio=2)  # 127 clusters, 42 exposed, weights 16.1
  des4 <- gen_sub_des4(data, ratio=2)  # 123 clusters, 37 exposed, weights 16.6
  result <- c(des5[[1]], des1[[1]], des2[[1]], des3[[1]], des4[[1]],
              des5[[2]], des1[[2]], des2[[2]], des3[[2]], des4[[2]])
  return(result)
}
results_alpha0 <- apply(oper, 2, mean)
proc.time() - ptm

rm(oper)

### Results for alpha = 0.3
ptm <- proc.time()
oper <- foreach (icount(MC_sims), .combine=rbind, .packages=c("SimCorMultRes", "evd", "dplyr", "geepack")) %dopar% {
  # simulate data
  data <- simulate_binary_data(K=734, nk=2, p_baseline=0.025, OR=2.05, p_exposed=0.2, gamma=0.5)
  
  # get variance estimates
  des5 <- gen_sub_des5(data, ratio=2)  # 80 clusters, 22 exposed, weights 17.1, 34.1
  des1 <- gen_sub_des1(data, ratio=2)  # 64 clusters, 35 exposed, weights 30.1
  des2 <- gen_sub_des2(data, ratio=2)  # 88 clusters, 44 exposed, weights 29.5
  des3 <- gen_sub_des3(data, ratio=2)  # 127 clusters, 42 exposed, weights 16.1
  des4 <- gen_sub_des4(data, ratio=2)  # 123 clusters, 37 exposed, weights 16.6
  result <- c(des5[[1]], des1[[1]], des2[[1]], des3[[1]], des4[[1]],
              des5[[2]], des1[[2]], des2[[2]], des3[[2]], des4[[2]])
  return(result)
}
results_alpha3 <- apply(oper, 2, mean)
proc.time() - ptm
stopCluster(cl)


### Save and display all results
results_alpha15 <- results[1:10]

results_all <- t(data.frame(results_alpha0, results_alpha15, results_alpha3))
rownames(results_all) <- paste0("$\\alpha$ = ", c(0, 0.15, 0.3))
colnames(results_all) <- rep(paste0("Design ", 0:4),2)

save(results_all, display, file="245_simulations.RData")

# Get information on simulated samples and weights for alpha = 0.15
display <- data.frame("Num_Clusters" = results[11:15],
                      "Max_Weight" = results[16:20],
                      "Num_Cases" = results[21:25],
                      "Num_Exposed" = results[26:30])
rownames(display) <- paste("Design", 0:4)












#========================= Numerical Analysis ===========================================

# calcVariance() computes the robust sandwich variance for the weighted GEE estimator 
# with binary covariate and outcome. All variables are constant across clusters, 
# either by design or in expectation.
#
# Args:
#   w0: sampling weight for controls
#   w1: sampling weight for cases
#   u0: risk of outcome given unexposed, P(Y=1|X=0)
#   u1: risk of outcome given exposed, P(Y=1|X=1)
#   alpha: within-cluster dependence
#   K: number of sample clusters
#   n: size of sample clusters
#   p00: proportion with X=0 and Y=0 per cluster
#   p01: proportion with X=0 and Y=1 per cluster
#   p10: proportion with X=1 and Y=0 per cluster
#   p11: proportion with X=1 and Y=1 per cluster
# 
# Returns a data frame containing: 
#   V(beta1): variance of beta1 estimator
#   V(beta2): variance of beta2 estimator
#   Cov: covariance between beta1 and beta2
#   beta1: log odds when X=0
#   beta2: log odds ratio
#   u0: baseline risk
#   u1: risk when X=1
#   v0: variance among unexposed (X=0)
#   v1: variance among exposed (X=1)
#   w0: sampling weight for controls
#   w1: sampling weight for cases
#   p0: proportion that is a control; vector
#   p1: proportion that is a case; vector
#   n: cluster sizes; vector
#   alpha: intracluster correlation
#   K: number of clusters
#   s1,s2,...,s8: components of Sigma_0 in the sandwich variance
calcVariance <- function(w0, w1=1, u0, u1, alpha, K, n, p00, p01, p10, p11) {
  
  # Calculate coefficients of interest
  beta1 <- log(u0/(1-u0))
  beta2 <- log(u0*(1-u1)/(u1*(1-u0)))
  
  # Calculate variance of Y for controls and cases
  v0 <- u0*(1-u0)
  v1 <- u1*(1-u1)
  
  # initialize variance estimator matrices
  Sigma1 <- matrix(0, nrow=2, ncol=2)
  Sigma0 <- matrix(0, nrow=2, ncol=2)
  # initialize variance component vectors
  s1 <- s2 <- s3 <- s4 <- s5 <- s6 <- s7 <- s8 <- a <- b <- c <- d <- numeric(K)
  
  # iterate over clusters
  for (k in 1:K) {
    # bread matrix
    P1_k <- matrix(c(p00[k]*v0*w0 + p01[k]*v0*w1 + p10[k]*v1*w0 + p11[k]*v1*w1, 
                     p10[k]*v1*w0 + p11[k]*v1*w1,
                     p10[k]*v1*w0 + p11[k]*v1*w1,
                     p10[k]*v1*w0 + p11[k]*v1*w1), nrow = 2, ncol = 2, byrow = TRUE)
    
    # meat matrix
    s1[k] <- w0^2*(1-alpha) + n[k]*(p00[k]*alpha*w0^2 + p01[k]*alpha*w0*w1 + 
                                      p10[k]*alpha*w0^2*sqrt(v1/v0) + p11[k]*alpha*w0*w1*sqrt(v1/v0))
    s2[k] <- w1^2*(1-alpha) + n[k]*(p01[k]*alpha*w1^2 + p00[k]*alpha*w0*w1 + 
                                      p11[k]*alpha*w1^2*sqrt(v1/v0) + p10[k]*alpha*w0*w1*sqrt(v1/v0))
    s3[k] <- w0^2*(1-alpha) + n[k]*(p10[k]*alpha*w0^2 + p11[k]*alpha*w0*w1 + 
                                      p00[k]*alpha*w0^2*sqrt(v0/v1) + p01[k]*alpha*w0*w1*sqrt(v0/v1))
    s4[k] <- w1^2*(1-alpha) + n[k]*(p11[k]*alpha*w1^2 + p10[k]*alpha*w0*w1 + 
                                      p01[k]*alpha*w1^2*sqrt(v0/v1) + p00[k]*alpha*w0*w1*sqrt(v0/v1))
    s5[k] <- n[k]*(p10[k]*alpha*w0^2*sqrt(v1/v0) + p11[k]*alpha*w0*w1*sqrt(v1/v0))
    s6[k] <- n[k]*(p11[k]*alpha*w1^2*sqrt(v1/v0) + p10[k]*alpha*w0*w1*sqrt(v1/v0))
    s7[k] <- w0^2*(1-alpha) + n[k]*(p10[k]*alpha*w0^2 + p11[k]*alpha*w0*w1)
    s8[k] <- w1^2*(1-alpha) + n[k]*(p11[k]*alpha*w1^2 + p10[k]*alpha*w0*w1)
    
    a[k] <- v0*(p00[k]*s1[k] + p01[k]*s2[k]) + v1*(p10[k]*s3[k] + p11[k]*s4[k])
    b[k] <- v1*(p10[k]*s3[k] + p11[k]*s4[k])
    c[k] <- v0*(p00[k]*s5[k] + p01[k]*s6[k]) + v1*(p10[k]*s7[k] + p11[k]*s8[k])
    d[k] <- v1*(p10[k]*s7[k] + p11[k]*s8[k])
    
    P0_k <- matrix(c(a[k], b[k], c[k], d[k]), nrow = 2, ncol = 2, byrow = TRUE)
    
    # Add to previous matrices
    Sigma1 <- Sigma1 + n[k]*P1_k
    Sigma0 <- Sigma0 + n[k]*P0_k
  }
  
  # Calculate robust sandwich variance matrix for coefficients of interest
  variance <- solve(Sigma1) %*% Sigma0 %*% solve(Sigma1)
  
  
  # Create output data frame for variance
  output <- as.data.frame(t(data.frame(variance[1,1], variance[2,2], variance[1,2])))
  
  rownames(output) <- c("Var($\\hat\\beta_1$)", "Var($\\hat\\beta_2$)", 
                        "Cov($\\hat\\beta_1, \\hat\\beta_2$)")
  colnames(output) <- "Value"
  
  # Create output data frame for 's' components
  s_comb <- as.data.frame(cbind(unique(s1), unique(s2), unique(s3), unique(s4), 
                                unique(s5), unique(s6), unique(s7), unique(s8)))
  colnames(s_comb) <- c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8")
  
  # Create output data frame for Sigma_0k components
  Sigma0_comp <- as.data.frame(cbind(unique(a), unique(b), unique(c), unique(d)))
  colnames(Sigma0_comp) <- c("a_k", "b_k", "c_k", "d_k")
  
  # Combine and return all outputs
  full_output <- list("Var" = round(output, 3), "Sigma1" = Sigma1, "Sigma0" = Sigma0, 
                      "s_unique" = round(s_comb, 3), "Sigma0_components" = round(Sigma0_comp,3))
  
  return(full_output)
  
}


# Args:
#   weights: vector of sampling weights, each entry for a cluster
#   u0: risk of outcome given unexposed, P(Y=1|X=0)
#   u1: risk of outcome given exposed, P(Y=1|X=1)
#   alpha: within-cluster dependence
#   K: number of sample clusters
#   n: size of sample clusters
#   p00: proportion with X=0 and Y=0 per cluster
#   p01: proportion with X=0 and Y=1 per cluster
#   p10: proportion with X=1 and Y=0 per cluster
#   p11: proportion with X=1 and Y=1 per cluster
# 
# Returns a data frame containing: 
#   V(beta1): variance of beta1 estimator
#   V(beta2): variance of beta2 estimator
#   Cov: covariance between beta1 and beta2
calcVariance2 <- function(weights, u0, u1, alpha, K, n, p00, p01, p10, p11) {
  
  # Calculate coefficients of interest
  beta1 <- log(u0/(1-u0))
  beta2 <- log(u0*(1-u1)/(u1*(1-u0)))
  
  # Calculate variance of Y for controls and cases
  v0 <- u0*(1-u0)
  v1 <- u1*(1-u1)
  
  # initialize variance estimator matrices
  Sigma1 <- matrix(0, nrow=2, ncol=2)
  Sigma0 <- matrix(0, nrow=2, ncol=2)
  
  # initialize variance component vectors
  p0 <- (p00+p01)
  p1 <- (p10+p11)
  a <- b <- c <- numeric(K)
  
  # iterate over clusters
  for (k in 1:K) {
    # bread matrix
    P1_k <- weights[k]*matrix(c(p0[k]*v0 + p1[k]*v1, p1[k]*v1, p1[k]*v1, p1[k]*v1), 
                              nrow = 2, ncol = 2, byrow = TRUE)
    
    # meat matrix
    a[k] <- weights[k]^2*(p1[k]*v1*(1+(n[k]*p1[k]-1)*alpha) + 2*n[k]*p0[k]*p1[k]*alpha*sqrt(v0*v1) + 
                            p0[k]*v0*(1+(n[k]*p0[k]-1)*alpha))
    b[k] <- weights[k]^2*p1[k]*v1*(1+(n[k]*p1[k]-1)*alpha) + n[k]*p0[k]*p1[k]*alpha*sqrt(v0*v1)
    c[k] <- weights[k]^2*p1[k]*v1*(1+(n[k]*p1[k]-1)*alpha)
    P0_k <- matrix(c(a[k], b[k], b[k], c[k]), nrow = 2, ncol = 2, byrow = TRUE)
    
    Sigma1 <- Sigma1 + n[k]*P1_k
    Sigma0 <- Sigma0 + n[k]*P0_k
  }
  
  # Calculate robust sandwich variance matrix for coefficients of interest
  variance <- solve(Sigma1) %*% Sigma0 %*% solve(Sigma1)
  
  # Create output data frame for variance
  output <- as.data.frame(t(data.frame(variance[1,1], variance[2,2], variance[1,2])))
  
  rownames(output) <- c("Var($\\hat\\beta_1$)", "Var($\\hat\\beta_2$)", 
                        "Cov($\\hat\\beta_1, \\hat\\beta_2$)")
  colnames(output) <- "Value"
  
  # Create output data frame for Sigma_0k components
  Sigma0_comp <- as.data.frame(cbind(unique(a), unique(b), unique(c)))
  colnames(Sigma0_comp) <- c("a_k", "b_k", "c_k")
  
  # Combine and return all outputs
  full_output <- list("Var" = round(output, 3), "Sigma1" = Sigma1, "Sigma0" = Sigma0, 
                      "Sigma0_components" = round(Sigma0_comp, 3))
  
  return(full_output)
}


# calcVarianceJoint() calculates the sandwich variance when both individual-level and cluster-level weights may be present,
# depending on the cluster type.
# 'weights': numeric vector that specifies weights for all clusters, with weights[k] equal to 0 if individual-level weights are 
#            required for cluster k. 
# 'w0' and 'w1' are the individual-level weights
# cluster: logical vector that specifies whether weights for a given cluster k are individual or cluster level.
# If cluster=TRUE, then 'weights' are used. Else, 'w0' and 'w1' are used
calcVarianceJoint <- function(w0, w1=1, weights, cluster, u0, u1, alpha, K, n, p00, p01, p10, p11) {
  
  # Calculate variance of Y for controls and cases
  v0 <- u0*(1-u0)
  v1 <- u1*(1-u1)
  
  # initialize variance estimator matrices
  Sigma1 <- matrix(0, nrow=2, ncol=2)
  Sigma0 <- matrix(0, nrow=2, ncol=2)
  # initialize variance component vectors
  s1 <- s2 <- s3 <- s4 <- s5 <- s6 <- s7 <- s8 <- a <- b <- c <- d <- numeric(K)
  
  # initialize variance component vectors
  p0 <- (p00+p01)
  p1 <- (p10+p11)
  a <- b <- c <- numeric(K)
  
  
  # iterate over clusters
  for (k in 1:K) {
    
    if (cluster[k]) { # cluster-level weights
      # bread matrix
      P1_k <- weights[k]*matrix(c(p0[k]*v0 + p1[k]*v1, p1[k]*v1, p1[k]*v1, p1[k]*v1), 
                                nrow = 2, ncol = 2, byrow = TRUE)
      
      # meat matrix
      a[k] <- weights[k]^2*(p1[k]*v1*(1+(n[k]*p1[k]-1)*alpha) + 2*n[k]*p0[k]*p1[k]*alpha*sqrt(v0*v1) + 
                              p0[k]*v0*(1+(n[k]*p0[k]-1)*alpha))
      b[k] <- weights[k]^2*p1[k]*v1*(1+(n[k]*p1[k]-1)*alpha) + n[k]*p0[k]*p1[k]*alpha*sqrt(v0*v1)
      c[k] <- weights[k]^2*p1[k]*v1*(1+(n[k]*p1[k]-1)*alpha)
      P0_k <- matrix(c(a[k], b[k], b[k], c[k]), nrow = 2, ncol = 2, byrow = TRUE)
      
    } else { # individual-level weights
      # bread matrix
      P1_k <- matrix(c(p00[k]*v0*w0 + p01[k]*v0*w1 + p10[k]*v1*w0 + p11[k]*v1*w1, 
                       p10[k]*v1*w0 + p11[k]*v1*w1,
                       p10[k]*v1*w0 + p11[k]*v1*w1,
                       p10[k]*v1*w0 + p11[k]*v1*w1), nrow = 2, ncol = 2, byrow = TRUE)
      
      # meat matrix
      s1[k] <- w0^2*(1-alpha) + n[k]*(p00[k]*alpha*w0^2 + p01[k]*alpha*w0*w1 + 
                                        p10[k]*alpha*w0^2*sqrt(v1/v0) + p11[k]*alpha*w0*w1*sqrt(v1/v0))
      s2[k] <- w1^2*(1-alpha) + n[k]*(p01[k]*alpha*w1^2 + p00[k]*alpha*w0*w1 + 
                                        p11[k]*alpha*w1^2*sqrt(v1/v0) + p10[k]*alpha*w0*w1*sqrt(v1/v0))
      s3[k] <- w0^2*(1-alpha) + n[k]*(p10[k]*alpha*w0^2 + p11[k]*alpha*w0*w1 + 
                                        p00[k]*alpha*w0^2*sqrt(v0/v1) + p01[k]*alpha*w0*w1*sqrt(v0/v1))
      s4[k] <- w1^2*(1-alpha) + n[k]*(p11[k]*alpha*w1^2 + p10[k]*alpha*w0*w1 + 
                                        p01[k]*alpha*w1^2*sqrt(v0/v1) + p00[k]*alpha*w0*w1*sqrt(v0/v1))
      s5[k] <- n[k]*(p10[k]*alpha*w0^2*sqrt(v1/v0) + p11[k]*alpha*w0*w1*sqrt(v1/v0))
      s6[k] <- n[k]*(p11[k]*alpha*w1^2*sqrt(v1/v0) + p10[k]*alpha*w0*w1*sqrt(v1/v0))
      s7[k] <- w0^2*(1-alpha) + n[k]*(p10[k]*alpha*w0^2 + p11[k]*alpha*w0*w1)
      s8[k] <- w1^2*(1-alpha) + n[k]*(p11[k]*alpha*w1^2 + p10[k]*alpha*w0*w1)
      
      a[k] <- v0*(p00[k]*s1[k] + p01[k]*s2[k]) + v1*(p10[k]*s3[k] + p11[k]*s4[k])
      b[k] <- v1*(p10[k]*s3[k] + p11[k]*s4[k])
      c[k] <- v0*(p00[k]*s5[k] + p01[k]*s6[k]) + v1*(p10[k]*s7[k] + p11[k]*s8[k])
      d[k] <- v1*(p10[k]*s7[k] + p11[k]*s8[k])
      
      P0_k <- matrix(c(a[k], b[k], c[k], d[k]), nrow = 2, ncol = 2, byrow = TRUE)
      
    }
    # Add to previous matrices
    Sigma1 <- Sigma1 + n[k]*P1_k
    Sigma0 <- Sigma0 + n[k]*P0_k
  }
  
  # Calculate robust sandwich variance matrix for coefficients of interest
  variance <- solve(Sigma1) %*% Sigma0 %*% solve(Sigma1)
  
  # Create output data frame for variance
  output <- as.data.frame(t(data.frame(variance[1,1], variance[2,2], variance[1,2])))
  
  rownames(output) <- c("Var($\\hat\\beta_1$)", "Var($\\hat\\beta_2$)", 
                        "Cov($\\hat\\beta_1, \\hat\\beta_2$)")
  colnames(output) <- "Value"
  
  ## Combine and return all outputs
  full_output <- list("Var" = round(output, 3), "Sigma1" = Sigma1, "Sigma0" = Sigma0)
  
  return(full_output)
}


alphas <- c(0, 0.15, 0.3)
variances <- as.data.frame(matrix(0, nrow = length(alphas), ncol = 5))
colnames(variances) <- c("Design 5", "Design 1", "Design 2", "Design 3", "Design 4")
rownames(variances) <- paste0("$\\alpha$ = ", alphas)
n1 <- 44

for (i in 1:length(alphas)) {
  alpha <- alphas[i]
  # Design 5 Scenario 1
  # Order: double case, single case, single control
  d5_s1 <- calcVarianceJoint(w0=734/88, w1=1, weights=c(rep(0, 5), rep(1, 39), rep(734/88*2, 83)), 
                             cluster=c(rep(FALSE, 5), rep(TRUE, 122)), u0=0.025, u1=0.05, alpha=alpha,
                             K=127, n=c(rep(2, 5), rep(1, 122)),
                             p00=c(rep(0.42, 5), rep(0, 39), rep(0.84, 83)), p01=c(rep(0.36, 5), rep(0.72, 39), rep(0, 83)),
                             p10=c(rep(0.08, 5), rep(0, 39), rep(0.16, 83)), p11=c(rep(0.14, 5), rep(0.28, 39), rep(0, 83)))
  
  # Design 1 Scenario 1
  # Order: double case, double control
  d1_s1 <- calcVariance2(weights=c(rep(1, n1), rep(690/22, n1/2)), u0=0.025, u1=0.05, alpha=alpha, 
                         K=3/2*n1, n=rep(2, 3/2*n1), 
                         p00=c(rep(0.42, n1), rep(0.84, n1/2)), p01=c(rep(0.36, n1), rep(0, n1/2)), 
                         p10=c(rep(0.08, n1), rep(0.16, n1/2)), p11=c(rep(0.14, n1), rep(0, n1/2)))
  
  # Design 2 Scenario 1
  # Order: double case, single control
  d2_s1 <- calcVariance2(weights=c(rep(1, n1), rep(690/44*2, n1)), u0=0.025, u1=0.05, alpha=alpha, 
                         K=2*n1, n=c(rep(2, n1), rep(1, n1)), 
                         p00=c(rep(0.42, n1), rep(0.84, n1)), p01=c(rep(0.36, n1), rep(0, n1)), 
                         p10=c(rep(0.08, n1), rep(0.16, n1)), p11=c(rep(0.14, n1), rep(0, n1)))
  
  # Design 3 Scenario 1
  # Order: single case, single control
  d3_s1 <- calcVariance2(weights=c(rep(1, n1), rep(690/88*2, 2*n1)), u0=0.025, u1=0.05, alpha=alpha, 
                         K=3*n1, n=c(rep(1, n1), rep(1, 2*n1)), 
                         p00=c(rep(0, n1), rep(0.84, 2*n1)), p01=c(rep(0.72, n1), rep(0, 2*n1)), 
                         p10=c(rep(0, n1), rep(0.16, 2*n1)), p11=c(rep(0.28, n1), rep(0, 2*n1)))
  
  # Design 4 Scenario 1
  # Order: double case, single case, double control, single control
  d4_s1 <- calcVariance(w0=1424/88, w1=1, u0=0.025, u1=0.05, alpha=alpha,
                        K=123, n=c(rep(2, 3), rep(1, 41), rep(2, 6), rep(1, 73)),
                        p00=c(rep(0.42, 3), rep(0, 41), rep(0.84, 6+73)), p01=c(rep(0.36, 3), rep(0.72, 41), rep(0, 6+73)),
                        p10=c(rep(0.08, 3), rep(0, 41), rep(0.16, 6+73)), p11=c(rep(0.14, 3), rep(0.28, 41), rep(0, 6+73)))
  
  variances[i, 1] <- d5_s1[[1]][2,1]  # Design 5
  variances[i, 2] <- d1_s1[[1]][2,1]  # Design 1
  variances[i, 3] <- d2_s1[[1]][2,1]  # Design 2
  variances[i, 4] <- d3_s1[[1]][2,1]  # Design 3
  variances[i, 5] <- d4_s1[[1]][2,1]  # Design 4
}

kable(variances, "latex", booktabs = TRUE, escape = FALSE, digits = 3, longtable = TRUE) %>%
  kable_styling(latex_options = "striped") %>%
  row_spec(c(3,6), bold = TRUE)

save(variances, file="245_numerical.RData")









#================================== Comparison of Results =======================================

# Load in simulations
load("C:/Users/Lang/Documents/Research/TB_Control_sampling/245_simulations.RData")
# Load in numerical analysis
load("C:/Users/Lang/Documents/Research/TB_Control_sampling/245_numerical.RData")

### Plot estimates on OR scale to check for bias 
estimates <- as.data.frame(exp(results_all[,6:10]))
estimates$Alpha <- c(0, 0.15, 0.3)
plot_data <- gather(estimates, "Design", "Estimate", 1:5)
ggplot(plot_data, aes(x=Alpha, y=Estimate, color=Design)) + geom_line() + 
  geom_hline(yintercept= 2.05, linetype = "dashed") + geom_point()

### Plot simulation results
results_all <- as.data.frame(results_all[, 1:5])
results_all$Alpha <- c(0, 0.15, 0.3)
plot_output <- gather(results_all, "Design", "Variance", 1:5)
plot_output$Type <- "Simulations"
plot2 <- ggplot(plot_output, aes(x=Alpha, y=Variance, color=Design)) + geom_line() + geom_point() +
  labs(title="Simulations") + theme(plot.title = element_text(hjust = 0.5))

### Plot numerical analysis results
variances$Alpha <- c(0, 0.15, 0.3)
plot_variances <- gather(variances, "Design", "Variance", 1:5)
plot3 <- ggplot(plot_variances, aes(x=Alpha, y=Variance, color=Design)) + geom_line() + geom_point() +
  labs(title="Numerical Analysis") + theme(plot.title = element_text(hjust = 0.5))

grid.arrange(plot3, plot2, ncol=2)

### Plot combined results
variances$Type <- "Numerical"
variances <- gather(variances, "Design", "Variance", 1:5)
plot_both <- rbind(plot_output, variances)
ggplot(plot_both, aes(x=Alpha, y=Variance, color=Design)) + 
  geom_line(aes(linetype=Type, color=Design)) + geom_point()
