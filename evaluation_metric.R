require(graphics)
source('helper_functions.R')
N<-2000
W_pos<-52500
W_neg<-8400
B1<-3300
B2<-2700
beta<-1/24
theta<-51450
# practice optim
fr <- function(x) {   ## Rosenbrock Banana function
  x1 <- x[1]
  x2 <- x[2]
  return(x1^2+2*x1 +1 + x2^2-4*x2 +1)
}

optim(c(1,5), fr)
(res <- optim(c(-1.2,1), fr, grr, method = "BFGS"))

############

data_RT <- simul_IDM(1000, prob_scaled_vec = prob_scaled_vec,
                     B1 = 3100, B2 = 2900, Ter = 0.3, h = 0.4)
data_RT_cg <- simul_IDM(1000, prob_scaled_vec = prob_scaled_vec,
                     B1 = 3100, B2 = 2900, Ter = 0.3, h = 0.4, sigma = 0.01, mode = "coarse-grained")

# quantiles:
q <- c(0.1, 0.3, 0.5, 0.7, 0.9)

# observed data:
(data_RT_q_c <- quantile(data_RT[data_RT$Response == 2, "RT"], probs = q))
(data_RT_q_e <- quantile(data_RT[data_RT$Response == 1, "RT"], probs = q))
(data_RT_cg_q_c <- quantile(data_RT_cg[data_RT_cg$Response == 2, "RT"], probs = q))
(data_RT_cg_q_e <- quantile(data_RT_cg[data_RT_cg$Response == 1, "RT"], probs = q))


plot(data_RT_q_c, q*prop.table(table(data_RT$Response))[2], pch = 2, 
     ylab = "Cumulative Probability", xlab = "Response Time (sec)", main = "Diffusion",
     xlim = c(0.3,0.6), ylim = c(0,1))
points(data_RT_q_e, q*prop.table(table(data_RT$Response))[1], pch = 2, col = "red")

points(data_RT_cg_q_c, q*prop.table(table(data_RT_cg$Response))[2], pch = 1)
points(data_RT_cg_q_e, q*prop.table(table(data_RT_cg$Response))[1], pch = 1, col = "red")

legend("topright", legend = c("coarse-grained (correct)", "Euler-Maruyama(correct)", "coarse-grained(error)", "Euler-Maruyama(error)"), 
       pch = c(2,2,1,1),  col = c("black","red","black","red"))




### chisqures to fit IDM Stijn's data to IDM Stijn's data.

data_RT <- get_pdf_RT(C = 0.2, h=0.4, Ter=0.3, Bns = 2500, Bs = 1000, mode = 1)

chis <- c()
accuracies <- c()
Cs <- seq(0,0.4,0.04)
for (C in Cs){
  simul_data <- get_pdf_RT(C = C, h=0.4, Ter=0.3, mode = 1)
  chi <- get_chisquare_stats_IDM(data_RT,simul_data)
  chis <- c(chis, chi)
  #accuracy <- sum(simul_data[,1])/500000
  #accuracies <- c(accuracies, accuracy)
  #print(chi)
  #print(accuracy)
}
plot(Cs, chis, xlab = "C", ylab = "Log Chisquare")
#plot(Cs, accuracies, xlab = "C", ylab = "Accuracy")

chis <- c()
Ters <- seq(0,0.6,0.05)
for (Ter in Ters){
  simul_data <- get_pdf_RT(C = 0.2, h=0.4, Ter=Ter, mode = 1)
  chi <- get_chisquare_stats_IDM(data_RT,simul_data)
  chis <- c(chis, chi)
  print(chi)
}
plot(Ters, chis, xlab = "Ter", ylab = "Log Chisquare")

chis <- c()
hs <- seq(0.2,0.6,0.04)
for (h in hs){
  simul_data <- get_pdf_RT(C = 0.2, h=h, Ter=0.3, mode = 1)
  chi <- get_chisquare_stats_IDM(data_RT,simul_data)
  chis <- c(chis, chi)
  print(chi)
}
plot(hs, chis, xlab = "h", ylab = "Log Chisquare")


### chisqures to fit IDM simulated data to IDM real data.
data_RT <- simul_IDM(2000, prob_scaled_vec = prob_scaled_vec,
                     C = 0.2, Ter = 0.3, h = 0.4)#, mode = "coarse-grained")
length(data_RT[data_RT$Response == 2, "RT"])/2000

chis <- c()
Cs <- seq(0,0.4,0.04)

for (C in Cs){
  simul_data <- get_pdf_RT(C = C, h=0.4, Ter=0.3, mode = 1)
  chi <- get_chisquare_stats(data_RT,simul_data)
  chis <- c(chis, chi)
}
plot(Cs, chis, xlab = "C", ylab = "Log Chisquare")


chis <- c()
Ters <- seq(0,0.6,0.06)
for (Ter in Ters){
  simul_data <- get_pdf_RT(C = 0.2, h=0.4, Ter=Ter, mode = 1)
  chi <- get_chisquare_stats(data_RT, simul_data)
  chis <- c(chis, chi)
}
plot(Ters, chis, xlab = "Ter", ylab = "Log Chisquare")

chis <- c()
hs <- seq(0.2,0.6,0.02)
for (h in hs){
  simul_data <- get_pdf_RT(C = 0.2, h=h, Ter=0.3, mode = 1)
  chi <- get_chisquare_stats(data_RT,simul_data)
  chis <- c(chis, chi)
}
plot(hs, chis, xlab = "h", ylab = "Log Chisquare")



#### to see whether quantile likelihood works for IDM simulation and fitting
data_RT <- simul_IDM(100000, prob_scaled_vec = prob_scaled_vec,
                     C = 0.2, Ter = 0.3, h = 0.4)#, mode = "coarse-grained")
prop.table(table(data_RT$Response))

log_ls <- c()
Cs <- seq(0,0.4,0.02)
for (C in Cs){
  pdf_RT <- get_pdf_RT(C = C, h=0.4, Ter=0.3, mode = 0)
  log_l <- IDM_likelihood(data_RT, pdf_RT)
  log_ls <- c(log_ls, log_l)
}
plot(Cs, log_ls, xlab = "C", ylab = "log likelihood")

log_ls <- c()
Ters <- seq(0,0.6,0.03)
for (Ter in Ters){
  pdf_RT <- get_pdf_RT(C = 0.5, h=0.4, Ter=Ter, mode = 0)
  log_l <- IDM_likelihood(data_RT, pdf_RT)
  log_ls <- c(log_ls, log_l)
}
plot(Ters, log_ls, xlab = "Ter", ylab = "log likelihood")

log_ls <- c()
hs <- seq(0.2,0.6,0.04)
for (h in hs){
  pdf_RT <- get_pdf_RT(C = 0.5, h=h, Ter=0.3, mode = 0)
  log_l <- IDM_likelihood(data_RT, pdf_RT)
  log_ls <- c(log_ls, log_l)
}
plot(hs, log_ls, xlab = "h", ylab = "log likelihood")



#### to see whether quantile likelihood works for DDM simulation and fitting
data=simdiffT(N=1000,a=1,mv=1,sv=0.3,ter=0.3)
real_data <- data.frame('RT' = data$rt, 'Response' = data$x+1)


chis <- c()
for (a in seq(0.1,2,0.1)){
  data=simdiffT(N=500000,a=a,mv=1,sv=0.3,ter=0.3)
  simul_data <- data.frame('RT' = data$rt, 'Response' = data$x+1)
  c_cdf <- ecdf(simul_data[simul_data$Response==2,]$RT)
  e_cdf <- ecdf(simul_data[simul_data$Response==1,]$RT)
  c_pdf <- cdf_pdf(c_cdf(seq(0.001, 3, 0.001)))
  e_pdf <- cdf_pdf(e_cdf(seq(0.001, 3, 0.001)))
  
  chi <-IDM_likelihood(real_data, cbind(c_pdf, e_pdf))
  print(chi)
  chis <- c(chis,chi)
}

plot(seq(0.1,2,0.1), chis, xlab = "a", ylab = "Log Likelihood")


chis <- c()
for (v in seq(0.1,2,0.2)){
  data=simdiffT(N=500000,a=1,mv=v,sv=0.3,ter=0.3)
  simul_data <- data.frame('RT' = data$rt, 'Response' = data$x+1)
  c_cdf <- ecdf(simul_data[simul_data$Response==2,]$RT)
  e_cdf <- ecdf(simul_data[simul_data$Response==1,]$RT)
  c_pdf <- cdf_pdf(c_cdf(seq(0.001, 3, 0.001)))
  e_pdf <- cdf_pdf(e_cdf(seq(0.001, 3, 0.001)))
  
  chi <-IDM_likelihood(real_data, cbind(c_pdf, e_pdf))
  print(chi)
  chis <- c(chis,chi)
}

plot(seq(0.1,2,0.2), chis, xlab = "v", ylab = "Log Likelihood")

chis <- c()
for (ter in seq(0,0.5,0.02)){
  data=simdiffT(N=500000,a=1,mv=1,sv=0.3,ter=ter)
  simul_data <- data.frame('RT' = data$rt, 'Response' = data$x+1)
  c_cdf <- ecdf(simul_data[simul_data$Response==2,]$RT)
  e_cdf <- ecdf(simul_data[simul_data$Response==1,]$RT)
  c_pdf <- cdf_pdf(c_cdf(seq(0.001, 3, 0.001)))
  e_pdf <- cdf_pdf(e_cdf(seq(0.001, 3, 0.001)))
  
  chi <-IDM_likelihood(real_data, cbind(c_pdf, e_pdf))
  print(chi)
  chis <- c(chis,chi)
}

plot(seq(0,0.5,0.02), chis, xlab = "Ter", ylab = "Log Likelihood")
