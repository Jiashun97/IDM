source('helper_functions.R')




## generate data from IDM and fit with DDM. 
# varying parameters in IDM to see changes in parameters in DDM

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






