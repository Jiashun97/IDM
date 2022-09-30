source('helper_functions.R')






# grid uniform initialization of C, Ter, and h
nSim <- 100
a <-c()
v <- c()
t0 <- c()
Cs <- runif(100,0,0.4)
Ters <- runif(100,0,0.6)
hs <- runif(100,0.2,0.6)

plot(0,0, xlim = c(0,0.4), ylim = c(0,6), col = "white",
     xlab = "Stimulus Distinctness", ylab = "Drift Rate")

for (i in 1:nSim){
  cat("simulation trial", i, ". ")
  start_time <- Sys.time()
  DDM_fit_results <- IDM_DDM(C=Cs[i], Ter=Ters[i], h=hs[i])
  DDM_fit_params <- DDM_fit_results$par
  if(DDM_fit_results$convergence == 0){
    points(Cs[i], DDM_fit_params["v"]) # plot one point to see online pattern
    
    a <- c(a, as.numeric(DDM_fit_params["a"]))
    v <- c(v, as.numeric(DDM_fit_params["v"]))
    t0 <- c(t0, as.numeric(DDM_fit_params["t0"]))
  }else if(DDM_fit_results$convergence != 0){
    print("Failed convergence")
  }
  print(Sys.time() - start_time)
}

cor(Cs, v)
plot(Cs, v)



{
par(mfrow = c(3, 3))
sample_size <- 100
plot(Cs[1:sample_size], a, xlim = c(0,0.4), xlab = "Stimulus Distinctness", ylab = "Boundary Seperation")
plot(Cs[1:sample_size], v, xlim = c(0,0.4), xlab = "Stimulus Distinctness", ylab = "Drift Rate")
plot(Cs[1:sample_size], t0, xlim = c(0,0.4), xlab = "Stimulus Distinctness", ylab = "Non-decision Time")

plot(Ters[1:sample_size], a, xlab = "Non-decision Time", ylab = "Boundary Seperation")
plot(Ters[1:sample_size], v, xlab = "Non-decision Time", ylab = "Drift Rate")
plot(Ters[1:sample_size], t0, xlab = "Non-decision Time", ylab = "Non-decision Time")

plot(hs[1:sample_size], a, xlab = "Detection Box Size", ylab = "Boundary Seperation")
plot(hs[1:sample_size], v, xlab = "Detection Box Size", ylab = "Drift Rate")
plot(hs[1:sample_size], t0, xlab = "Detection Box Size", ylab = "Non-decision Time")
}
  
plot(Cs[1:sample_size], v, xlab = "Stimulus Distinctness", ylab = "Drift Rate")
plot(Cs[1:sample_size],v/a, xlab = "Stimulus Distinctness", ylab = "Drift Rate")

plot(hs[1:sample_size], a, xlab = "Detection Box Size", ylab = "Boundary Seperation")
plot(1/hs[1:sample_size], a, xlab = "Detection Box Size", ylab = "Boundary Seperation")

cor(Cs[1:sample_size],v)
cor(Cs[1:sample_size],v/a)

results <- data.frame(Cs = Cs[1:sample_size], a = a, v = v, t0 = t0)
pcor(results)

simul_IDM(1000, prob_scaled_vec = prob_scaled_vec)





