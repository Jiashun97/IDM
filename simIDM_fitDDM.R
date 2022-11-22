source('helper_functions.R')






# grid uniform initialization of C, Ter, and h
nSim <- 100
a <-c()
v <- c()
t0 <- c()
Cs <- runif(100,0.05,0.4)
hs <- runif(100,0.2,0.8)
Ters <- runif(100,0.1,0.6)


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


# full axis labels
{
par(mfrow = c(3, 3))
sample_size <- 100
plot(Cs[1:sample_size], v, xlim = c(0,0.4), xlab = "Stimulus Distinctness", ylab = "Drift Rate")
plot(hs[1:sample_size], v, xlab = "Detection Box Size", ylab = "Drift Rate")
plot(Ters[1:sample_size], v, xlab = "Non-decision Time", ylab = "Drift Rate")


plot(Cs[1:sample_size], a, xlab = "Stimulus Distinctness", ylab = "Boundary Seperation")
plot(hs[1:sample_size], a, xlab = "Detection Box Size", ylab = "Boundary Seperation")
plot(Ters[1:sample_size], a, xlab = "Non-decision Time", ylab = "Boundary Seperation")

plot(Cs[1:sample_size], t0, xlim = c(0,0.4), xlab = "Stimulus Distinctness", ylab = "Non-decision Time")
plot(hs[1:sample_size], t0, xlab = "Detection Box Size", ylab = "Non-decision Time")
plot(Ters[1:sample_size], t0, xlab = "Non-decision Time", ylab = "Non-decision Time")
}

# marginal labels
{
  plot_mar <- 3
  par(mfrow = c(3, 3))
  par(mar = c(plot_mar, plot_mar, 0.5, 0.5))
  
  plot(Cs, v, xlim = c(0.05,0.4), xlab = "", ylim = c(0,3), ylab = "Drift Rate")
  mtext("Drift Rate", side = 2, line = 2, cex = 0.7)
  plot(hs, v, xlim = c(0.2,0.6), xlab = "", ylim = c(0,3), ylab = "")
  plot(Ters, v, xlim = c(0.1,0.6), xlab = "", ylim = c(0,3), ylab = "")
  
  plot(Cs, a, xlim = c(0.05,0.4), xlab = "", ylim = c(0.4,0.9), ylab = "Boundary Seperation")
  mtext("Boundary Seperation", side = 2, line = 2, cex = 0.7)
  plot(hs, a, xlim = c(0.2,0.6), xlab = "", ylim = c(0.4,0.9), ylab = "")
  plot(Ters, a, xlim = c(0.1,0.6), xlab = "", ylim = c(0.4,0.9), ylab = "")
  
  plot(Cs, t0, xlim = c(0.05,0.4), xlab = "", ylim = c(0.1,0.7), ylab = "Non-decision Time")
  mtext("Non-decision Time", side = 2, line = 2, cex = 0.7)
  mtext("Stimulus Distinctness", side = 1, line = 2, cex = 0.7)
  plot(hs, t0, xlim = c(0.2,0.6), xlab = "", ylim = c(0.1,0.7), ylab = "")
  mtext("Detection Box Size", side = 1, line = 2, cex = 0.7)
  plot(Ters, t0, xlim = c(0.1,0.6), xlab = "", ylim = c(0.1,0.7), ylab = "")
  mtext("Non-decision Time", side = 1, line = 2, cex = 0.7)
}
#bottom_text <- c("Drift Rate", "Drift Rate", "Drift Rate")
  

length(Cs)
length(v)
  
plot(Cs[1:sample_size], v, xlab = "Stimulus Distinctness", ylab = "Drift Rate")
plot(Cs[1:sample_size],v/a, xlab = "Stimulus Distinctness", ylab = "Drift Rate")

plot(hs[1:sample_size], a, xlab = "Detection Box Size", ylab = "Boundary Seperation")
plot(1/hs[1:sample_size], a, xlab = "Detection Box Size", ylab = "Boundary Seperation")

cor(Cs[1:sample_size],v)
cor(Cs[1:sample_size],v/a)

results <- data.frame(Cs = Cs, Ters = Ters, hs = hs, a = a, v = v, t0 = t0)
write.csv(results,'IDM_DDM_result.csv')

simul_IDM(1000, prob_scaled_vec = prob_scaled_vec)



df <- read.csv("/Users/wangjiashun/Documents/GitHub/IDM/IDM_DDM_result.csv")
Cs <- df$Cs
hs <- df$hs
Ters <- df$Ters
v <- df$v
a <- df$a
t0 <- df$t0

par(mfrow = c(3, 3))
plot(df$c, df$a, xlim = c(0,0.4), xlab = "Stimulus Distinctness", ylab = "Boundary Seperation")
plot(df$c, df$v, xlim = c(0,0.4), xlab = "Stimulus Distinctness", ylab = "Drift Rate")
plot(df$c, df$t0, xlab = "Stimulus Distinctness", ylab = "Non-decision Time")

plot(df$ter, df$a, xlab = "Non-decision Time", ylab = "Boundary Seperation")
plot(df$ter, df$v, xlab = "Non-decision Time", ylab = "Drift Rate")
plot(df$ter, df$t0, xlab = "Non-decision Time", ylab = "Non-decision Time")

plot(df$h, df$a, xlab = "Detection Box Size", ylab = "Boundary Seperation")
plot(df$h, df$v, xlab = "Detection Box Size", ylab = "Drift Rate")
plot(df$h, df$t0, xlab = "Detection Box Size", ylab = "Non-decision Time")

