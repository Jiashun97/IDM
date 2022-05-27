source('helper_functions.R')
library(ppcor)




## generate data from IDM and fit with DDM. 
# varying parameters in IDM to see changes in parameters in DDM

# Bns <- 2500, Bs <- 500, to make it symmetric over 3000



a <- c()
v <- c()
t0 <- c()
sz <- c()
st0 <- c()
sv <- c()
Cs <- seq(0.2,1,0.05)

iter <- 1
for (C in Cs){
  for (i in 1:1){
    DDM_fit_params <- IDM_DDM(C)
    #for (name in names(DDM_fit_params)){
    a_temp <- c(a, as.numeric(DDM_fit_params["a"]))
    v_temp <- c(v, as.numeric(DDM_fit_params["v"]))
    t0_temp <- c(t0, as.numeric(DDM_fit_params["t0"]))
    #sz_temp <- c(sz, as.numeric(DDM_fit_params["sz"]))
    #st0_temp <- c(st0, as.numeric(DDM_fit_params["st0"]))
    sv_temp <- c(sv, as.numeric(DDM_fit_params["sv"]))
    print(iter)
    iter <- iter + 1
  }
  a <- c(a, mean(a_temp))
  v <- c(v, mean(v_temp))
  t0 <- c(t0, mean(t0_temp))
  sz <- c(sz, mean(sz_temp))
  st0 <- c(st0, mean(st0_temp))
  sv <- c(sv, mean(sv_temp))
  
}

vary_C <- list("a" = a, "v" = v, "t0" = t0)
plot(Cs, a, xlab = "Stimulus Distinctness", ylab = "Boundary Seperation")
plot(Cs, v, xlab = "Stimulus Distinctness", ylab = "Drift Rate")
plot(Cs, t0, xlab = "Stimulus Distinctness", ylab = "Non-decision Time")
plot(Cs, sz)
plot(Cs, st0)
plot(Cs, sv)


#
cor(C, v)

#
# freeze other params? like a  decision_boundary-a, Ter-t0
# vary Ter

a <- c()
v <- c()
t0 <- c()
sz <- c()
st0 <- c()
sv <- c()
Ters <- seq(0.3,0.6,0.03)

iter <- 1
for (Ter in Ters){
  for (i in 1:1){
    DDM_fit_params <- IDM_DDM(Ter = Ter)
    #for (name in names(DDM_fit_params)){
    a_temp <- c(a, as.numeric(DDM_fit_params["a"]))
    v_temp <- c(v, as.numeric(DDM_fit_params["v"]))
    t0_temp <- c(t0, as.numeric(DDM_fit_params["t0"]))
    #sz_temp <- c(sz, as.numeric(DDM_fit_params["sz"]))
    #st0_temp <- c(st0, as.numeric(DDM_fit_params["st0"]))
    #sv_temp <- c(sv, as.numeric(DDM_fit_params["sv"]))
    print(iter)
    iter <- iter + 1
  }
  a <- c(a, mean(a_temp))
  v <- c(v, mean(v_temp))
  t0 <- c(t0, mean(t0_temp))
  sz <- c(sz, mean(sz_temp))
  st0 <- c(st0, mean(st0_temp))
  sv <- c(sv, mean(sv_temp))
  
}
vary_Ter <- list("a" = a, "v" = v, "t0" = t0)
plot(Ters, a, xlab = "Non-decision Time", ylab = "Boundary Seperation")
plot(Ters, v, xlab = "Non-decision Time", ylab = "Drift Rate")
plot(Ters, t0, xlab = "Non-decision Time", ylab = "Non-decision Time")
plot(Ters, sz)
plot(Ters, st0)
plot(Ters, sv)

# vary h

a <- c()
v <- c()
t0 <- c()
sz <- c()
st0 <- c()
sv <- c()
hs <- seq(0.3,0.7,0.05)

iter <- 1
for (h in hs){
  for (i in 1:1){
    DDM_fit_params <- IDM_DDM(h = h)
    #for (name in names(DDM_fit_params)){
    a_temp <- c(a, as.numeric(DDM_fit_params["a"]))
    v_temp <- c(v, as.numeric(DDM_fit_params["v"]))
    t0_temp <- c(t0, as.numeric(DDM_fit_params["t0"]))
    sz_temp <- c(sz, as.numeric(DDM_fit_params["sz"]))
    st0_temp <- c(st0, as.numeric(DDM_fit_params["st0"]))
    sv_temp <- c(sv, as.numeric(DDM_fit_params["sv"]))
    print(iter)
    iter <- iter + 1
  }
  a <- c(a, mean(a_temp))
  v <- c(v, mean(v_temp))
  t0 <- c(t0, mean(t0_temp))
  sz <- c(sz, mean(sz_temp))
  st0 <- c(st0, mean(st0_temp))
  sv <- c(sv, mean(sv_temp))
  
}

vary_h <- list("a" = a, "v" = v, "t0" = t0)

plot(hs, a, xlab = "Detection Box Size", ylab = "Boundary Seperation")
plot(hs, v, xlab = "Detection Box Size", ylab = "Drift Rate")
plot(hs, t0, xlab = "Detection Box Size", ylab = "Non-decision Time")
plot(hs, sz)
plot(hs, st0)
plot(hs, sv)

# plot all results together
par(mfrow = c(3, 3))

plot(Cs, vary_C[["a"]], xlab = "Stimulus Distinctness", ylab = "Boundary Seperation")
plot(Cs, vary_C[["v"]], xlab = "Stimulus Distinctness", ylab = "Drift Rate")
plot(Cs, vary_C[["t0"]], xlab = "Stimulus Distinctness", ylab = "Non-decision Time")

plot(Ters, vary_Ter[["a"]], xlab = "Non-decision Time", ylab = "Boundary Seperation")
plot(Ters, vary_Ter[["v"]], xlab = "Non-decision Time", ylab = "Drift Rate")
plot(Ters, vary_Ter[["t0"]], xlab = "Non-decision Time", ylab = "Non-decision Time")

plot(hs, vary_h[["a"]], xlab = "Detection Box Size", ylab = "Boundary Seperation")
plot(hs, vary_h[["v"]], xlab = "Detection Box Size", ylab = "Drift Rate")
plot(hs, vary_h[["t0"]], xlab = "Detection Box Size", ylab = "Non-decision Time")




# grid uniform initialization of C, Ter, and h
Cs <- c()
Ters <- c()
hs <- c()
a <- c()
v <- c()
t0 <- c()


plot(0,0, v, xlim = c(-1,1), ylim = c(-1,10), col = "white",
     xlab = "Stimulus Distinctness", ylab = "Drift Rate")

nSim <- 200
sim <- 0
for (i in 1:nSim){
  start_time <- Sys.time()
  C <- runif(1,-1,1)
  Ter <- runif(1,0.3,0.6)
  h <- runif(1,0.3,0.7)
  
  Cs <- c(Cs, C)
  Ters <- c(Ters, Ter)
  hs <- c(hs, h)
  
  DDM_fit_results <- IDM_DDM(C=C, Ter=Ter, h=h)
  DDM_fit_params <- DDM_fit_results$par
  if(DDM_fit_results$convergence == 0){
    points(C, DDM_fit_params["v"])
    
    #plot(Ters, vary_Ter[["t0"]], xlab = "Non-decision Time", ylab = "Non-decision Time")
    
    #plot(hs, vary_h[["a"]], xlab = "Detection Box Size", ylab = "Boundary Seperation")
    
    a <- c(a, as.numeric(DDM_fit_params["a"]))
    v <- c(v, as.numeric(DDM_fit_params["v"]))
    t0 <- c(t0, as.numeric(DDM_fit_params["t0"]))
  }else if(DDM_fit_results$convergence == 0){
    print("Failed convergence")
  }
  
  sim <- sim + 1
  cat("simulation trial", sim)
  print(Sys.time() - start_time)
}

par(mfrow = c(3, 3))
sample_size <- 38
plot(Cs[1:sample_size], a, xlab = "Stimulus Distinctness", ylab = "Boundary Seperation")
plot(Cs[1:sample_size], v, xlab = "Stimulus Distinctness", ylab = "Drift Rate")
plot(Cs[1:sample_size], t0, xlab = "Stimulus Distinctness", ylab = "Non-decision Time")

plot(Ters[1:sample_size], a, xlab = "Non-decision Time", ylab = "Boundary Seperation")
plot(Ters[1:sample_size], v, xlab = "Non-decision Time", ylab = "Drift Rate")
plot(Ters[1:sample_size], t0, xlab = "Non-decision Time", ylab = "Non-decision Time")

plot(hs[1:sample_size], a, xlab = "Detection Box Size", ylab = "Boundary Seperation")
plot(hs[1:sample_size], v, xlab = "Detection Box Size", ylab = "Drift Rate")
plot(hs[1:sample_size], t0, xlab = "Detection Box Size", ylab = "Non-decision Time")

  
plot(Cs[1:sample_size], v, xlab = "Stimulus Distinctness", ylab = "Drift Rate")
plot(Cs[1:sample_size],v/a, xlab = "Stimulus Distinctness", ylab = "Drift Rate")

plot(hs[1:sample_size], a, xlab = "Detection Box Size", ylab = "Boundary Seperation")
plot(1/hs[1:sample_size], a, xlab = "Detection Box Size", ylab = "Boundary Seperation")

cor(Cs[1:sample_size],v)
cor(Cs[1:sample_size],v/a)

results <- data.frame(Cs = Cs[1:sample_size], a = a, v = v, t0 = t0)
pcor(results)


## generate data from DDM and fit with IDM. 
# varying parameters in DDM to see changes in parameters in IDM





## try to make data generated from DDM similar to what it looks like in Urai et al. (2019).



## IDM for speed accuracy trade-off?