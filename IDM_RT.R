
#  test whether your RTDist installation is capable 
# of running simulations on the CPU
source('sourceAllRTDist.R')

testWrappers()


settings<-list()
settings$nWalkers<-500000 # number of simulations
settings$nGPUs<-0
settings$gpuIds<-0
settings<- checkSettings(settings)

#LDM
#LDMPars=list()
#LDMPars$a=0.08
#LDMPars$nStimuli=2
#LDMPars$v=matrix(c(0,0.1),1,2)
#LDMPars= checkLDMPars(LDMPars)
#out = LDMDist(LDMPars,settings) # "distributions" "extra" "ok"  "used"  
#plot(matrix(1:3000,3000, 2*LDMPars$nStimuli), out$distributions,'l')

# IDM
IDMPars<-list()
IDMPars$nDim <- 2
IDMPars$nStimuli <- 1
IDMPars$N <- matrix(c(1000,1000),2,1)
IDMPars$W <- matrix(c(52500,8400,8400, 52500),2,2)
IDMPars$h <- matrix(c(0.4,0.4,0.4, 0.4),2,2)
IDMPars$Theta <- matrix(c(51450,51450),2,1)
#IDMPars$B <- matrix(c(3000,3000,3000,3000),2,2) # depends on number of stimuli 
IDMPars$B <- matrix(c(3300,2700),2,1) # depends on number of stimuli 

IDMPars$spontaneousTime <- 1
IDMPars$dynamics <- 1
IDMPars$D <- matrix(c(0.05,0.05),2,1)  # original
IDMPars$sigma <- matrix(c(0.02,0.02),2,1)
IDMPars$deltat <- 0.001
  


IDMPars<- checkIDMPars(IDMPars)
out <- IDMDist(IDMPars,settings) # "distributions" "extra" "ok"  "used"  
plot(matrix(1:3000,3000, 2*IDMPars$nStimuli), out$distributions,'l',
     xlab = "RT bin index", ylab = "frequency")


# Plot
plot(1:3000, out$distributions[,1], col="red", type = "l",
     xlab = "RT bin index", ylab = "frequency")
lines(1:3000, out$distributions[,2], col="blue", type = "l")
legend(1500, 8000, legend=c("Response 1", "Response 2"),
       col=c("red", "blue"), lty=1:1, cex=0.8)

# truncated
plot(1:500, out$distributions[1:500,1], col="red", type = "l",
     xlab = "RT bin index", ylab = "frequency")
lines(1:500, out$distributions[1:500,2], col="blue", type = "l")
legend(300, 8000, legend=c("Response 1", "Response 2"),
       col=c("red", "blue"), lty=1:1, cex=0.8)

rt_dist <- get_pdf_RT()
plot(matrix(1:3000,3000, 2*IDMPars$nStimuli), rt_dist/length(sample_rt_c),'l',
     xlab = "RT bin index", ylab = "frequency")


# quantiles:
q <- c(0.1, 0.3, 0.5, 0.7, 0.9)

sample_rt_c <- c()
for (i in 1:3000){
  sample_rt_c <- c(sample_rt_c, rep(i, (rt_dist[,1][i])))
}

sample_rt_e <- c()
for (i in 1:3000){
  sample_rt_e <- c(sample_rt_e, rep(i, (rt_dist[,2][i])))
}



### get quantile from rt distribution
get_quantile  <- function(pdf){
  sample_rt <- c()
  for (i in 1:3000){
    sample_rt <- c(sample_rt, rep(i, (pdf[i])))
  }
  quantile <- quantile(sample_rt, probs = q)
  return(quantile)
}

q_c <- as.numeric(get_quantile(rt_dist[,1]))
q_e <- as.numeric(get_quantile(rt_dist[,2]))


plot(ecdf(sample_rt_c))




cdf_c <- pdf_cdf(rt_dist[,1])
cdf_e <- pdf_cdf(rt_dist[,2])

plot(cdf_c, type = "l")

q_expected_c <- cdf_c[q_c]
q_expected_e <- cdf_e[q_e]



data_RT <- simul_IDM(nsim=500, prob_scaled_vec, B1 = 3100, B2 = 2900, Ter = 0.3, h = 0.4)
rt_dist <- get_pdf_RT(C=0.5, h=0.4, Ter=0.3)

prop.table(table(data_RT$Response))
get_chisquare_stats(data_RT, rt_dist, nsim = 1000)



#######
fit_IDM <- function(data, param_init = c(0.5, 0.3, 0.4), iter = 10){ # param = c(C, Ter, h)
  
}

param_mat <- matrix(rep(param_init,4),4,3,byrow=TRUE) + matrix(runif(12,0,0.1),4,3,byrow=TRUE)
chisquare_test(C = 0.5, h = 0.4, Ter = 0.3)


### For DDM, try to compare DDM model first to see whether the code is correct
# I 
data=simdiffT(N=1000,a=1,mv=1,sv=0.3,ter=0.3)
real_data <- data.frame('RT' = data$rt, 'Response' = data$x+1)

data=simdiffT(N=500000,a=1,mv=1,sv=0.3,ter=0.3)
simul_data <- data.frame('RT' = data$rt, 'Response' = data$x+1)

c_cdf <- ecdf(simul_data[simul_data$Response==2,]$RT)
e_cdf <- ecdf(simul_data[simul_data$Response==1,]$RT)

get_chisquare_stats(real_data, cbind(c_cdf(seq(0.001, 3, 0.001)), e_cdf(seq(0.001, 3, 0.001))), 
                    nsim = 1000)
chis <- c()
for (a in seq(0.1,2,0.1)){
  data=simdiffT(N=500000,a=a,mv=1,sv=0.3,ter=0.3)
  simul_data <- data.frame('RT' = data$rt, 'Response' = data$x+1)
  c_cdf <- ecdf(simul_data[simul_data$Response==2,]$RT)
  e_cdf <- ecdf(simul_data[simul_data$Response==1,]$RT)
  
  chi <-get_chisquare_stats(real_data, cbind(c_cdf(seq(0.001, 3, 0.001)), e_cdf(seq(0.001, 3, 0.001))), 
                      nsim = 1000)
  print(chi)
  chis <- c(chis,chi)
}

plot(seq(0.1,2,0.1), chis, xlab = "a", ylab = "Chisquare")

chis <- c()
for (v in seq(0.1,2,0.1)){
  data=simdiffT(N=500000,a=1,mv=v,sv=0.3,ter=0.3)
  simul_data <- data.frame('RT' = data$rt, 'Response' = data$x+1)
  c_cdf <- ecdf(simul_data[simul_data$Response==2,]$RT)
  e_cdf <- ecdf(simul_data[simul_data$Response==1,]$RT)
  
  chi <-get_chisquare_stats(real_data, cbind(c_cdf(seq(0.001, 3, 0.001)), e_cdf(seq(0.001, 3, 0.001))), 
                            nsim = 1000)
  print(chi)
  chis <- c(chis,chi)
}

plot(seq(0.1,2,0.1), chis, xlab = "v", ylab = "Chisquare")

chis <- c()
for (ter in seq(0,0.5,0.05)){
  data=simdiffT(N=500000,a=1,mv=1,sv=0.3,ter=ter)
  simul_data <- data.frame('RT' = data$rt, 'Response' = data$x+1)
  c_cdf <- ecdf(simul_data[simul_data$Response==2,]$RT)
  e_cdf <- ecdf(simul_data[simul_data$Response==1,]$RT)
  
  chi <-get_chisquare_stats(real_data, cbind(c_cdf(seq(0.001, 3, 0.001)), e_cdf(seq(0.001, 3, 0.001))), 
                            nsim = 1000)
  print(chi)
  chis <- c(chis,chi)
}

plot(seq(0,0.5,0.05), chis, xlab = "Ter", ylab = "Chisquare")


# see how chisquare stats changes when change rt data
data=simdiffT(N=500000,a=1,mv=1,sv=0.3,ter=0.3)
simul_data <- data.frame('RT' = data$rt, 'Response' = data$x+1)

c_cdf <- epdf(simul_data[simul_data$Response==2,]$RT)
e_cdf <- epdf(simul_data[simul_data$Response==1,]$RT)

timefit(simul_data[simul_data$Response==2,]$RT)
simul_data <- cbind(c_cdf(seq(0.001, 3, 0.001)), e_cdf(seq(0.001, 3, 0.001)))

chis <- c()
iter <- 0
for (a in seq(0.2,2,0.1)){
  data=simdiffT(N=1000,a=a,mv=1,sv=0.3,ter=0.3)
  real_data <- data.frame('RT' = data$rt, 'Response' = data$x+1)
  chi <- get_chisquare_stats(real_data, simul_data, 
                      nsim = 1000)
  
  chis <- c(chis, chi)
  iter <- iter + 1
  cat("iteration: ", iter, ". ")
}

plot(seq(0.2,2,0.1),chis, xlab = "drift rate", ylab = "chisqure")

hist(simul_data[simul_data$Response==2,]$RT)


# `DDM try optim function to recover the parameters

get_chisquare_stats_DDM <- function(param, data_RT, nsim = 1000){

  q = c(0.1, 0.3, 0.5, 0.7, 0.9)
  p = c(.1,.2,.2,.2,.2,.1)
  a <- param[1]
  v <- param[2]
  sv <- param[3]
  ter <- param[4]
  # empirical RT quantile
  data_RT_q_c <- as.numeric(quantile(data_RT[data_RT$Response == 2, "RT"], probs = q))
  data_RT_q_e <- as.numeric(quantile(data_RT[data_RT$Response == 1, "RT"], probs = q))
  
  data=simdiffT(N=1000,a=a,mv=v,sv=sv,ter=ter)
  simul_data <- data.frame('RT' = data$rt, 'Response' = data$x+1)
  
  cdf_c <- ecdf(simul_data[simul_data$Response==2,]$RT)
  cdf_e <- ecdf(simul_data[simul_data$Response==1,]$RT)
  
  #cdf_c <- rt_dist[,1]
  #cdf_e <- rt_dist[,2]
  
  # expected cumulative probability
  expected_cp_c <- cdf_c(data_RT_q_c*1000) # from s to index of ms
  expected_cp_e <- cdf_e(data_RT_q_c*1000)
  
  expected_cp_c <- c(0, expected_cp_c, 1)
  expected_cp_e <- c(0, expected_cp_e, 1)
  
  expected_freq_c <- numeric(6)
  expected_freq_e <- numeric(6)
  for (i in 1:6){
    expected_freq_c[i] <- expected_cp_c[i+1] - expected_cp_c[i]
    expected_freq_e[i] <- expected_cp_e[i+1] - expected_cp_e[i]
  }
  expected_freq_c <- expected_freq_c*nsim
  expected_freq_e <- expected_freq_e*nsim
  
  chisquare_c <- sum((expected_freq_c-p*nsim)^2/p*nsim)
  chisquare_e <- sum((expected_freq_e-p*nsim)^2/p*nsim)
  
  log(chisquare_c + chisquare_e)
}
data <- simdiffT(N=1000,a=2,mv=1,sv=0.3,ter=0.3)
real_data <- data.frame('RT' = data$rt, 'Response' = data$x+1)

optim(param = c(0.2, 0.2, 0.1, 0.4), get_chisquare_stats_DDM(data_RT = real_data), 
      lower = 0, method = "L-BFGS-B")

get_chisquare_stats_DDM(param = c(0.2, 0.2, 0.1, 0.4), data_RT=real_data)
