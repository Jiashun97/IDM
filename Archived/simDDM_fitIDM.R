source('helper_functions.R')




## generate data from DDM and fit with IDM. 
# varying parameters in DDM to see changes in parameters in IDM

DDM_data <- simdiffT(N=10,a=1,mv=1,sv=0.3,ter=0.3)
DDM_data <- data.frame('RT' = DDM_data$rt, 'Response' = DDM_data$x+1)


data_RT <- simul_IDM(3000, prob_scaled_vec = prob_scaled_vec,
                     C = 0.2, h = 0.4, Ter = 0.3, mode = "coarse-grained")
data <- data_RT

#f <- function(data, par){
par <- c(0.2,0.3,0.4)
f <- function(par){
    
  C <- par[1]
  h <- par[2]
  Ter <- par[3]
  
  simul_data <- get_pdf_RT(C = C, h=h, Ter=Ter, Bns = 2500, Bs = 1000, mode = 0)
  
  #with(data, get_chisquare_stats(data, simul_data))
  with(data, -IDM_likelihood(data, simul_data))
  
}

f(par = c(0.2,0.4,0.3))
f(par = c(0.1,0.2,0.2))

## try to make data generated from DDM similar to what it looks like in Urai et al. (2019).

## use mode finding algorithm
#outDEoptim <- DEoptim(fn = f,
#                      lower=c(0, 0, 0),
#                      upper=c(1, 1, 1),
#                      DEoptim.control(itermax = 10)
#                      )
# print output information
#summary(outDEoptim)

# plot the best members
#plot(outDEoptim, type = 'b')

# best parameter
#outDEoptim$member$bestmemit[10,]
#f(par = outDEoptim$member$bestmemit[10,]) # likelihood


# Create a general prior distribution by specifying an arbitrary density function and a
# corresponding sampling function
density <- function(par){
  d1 <- dunif(par[1], 0, 0.5, log =TRUE)
  d2 <- dunif(par[2], 0, 0.7, log =TRUE)
  d3 <- dunif(par[3], 0, 0.7, log =TRUE)
  return(d1 + d2 + d3)
}

# The sampling is optional but recommended because the MCMCs can generate automatic starting
# conditions if this is provided
sampler <- function(n=1){
  d1 <- runif(n, 0, 0.5)
  d2 <- runif(n, 0, 0.7)
  d3 <- runif(n, 0, 0.7)
  return(cbind(d1,d2,d3))
}

prior <- createPrior(density = density, 
                     sampler = sampler,
                     lower = c(0, 0, 0), 
                     upper = c(0.5, 0.7, 0.7))

bayesianSetup = createBayesianSetup(likelihood = f, prior = prior)


ddm <- c()
idm <- c()
## run DRAM to fit IDM 
for (i in 1:30){

  mv <- runif(1, 0, 4)
  a <- runif(1, 0, 0.8)
  ter <- runif(1, 0, 0.7)
  
  data1 <- simdiffT(N=10000,a=a,mv=mv,sv=0.3,ter=ter)
  data <- data.frame('RT' = data1$rt, 'Response' = data1$x+1)
  
  
  outDEoptim <- DEoptim(fn = f,
                        lower=c(0, 0, 0),
                        upper=c(0.5, 0.7, 0.7),
                        DEoptim.control(itermax = 10))
  startValues <- as.vector(outDEoptim$member$bestmemit[10,]) # .5139145 0.8287143 0.2542243
  #startValues <- c(.5139145, 0.8287143, 0.2542243)
  settings <- list(nrChains= 1, adapt = T, DRlevels = 2,
                   iterations = 10000, 
                   gibbsProbabilities = NULL, temperingFunction = NULL,
                   optimize = T,  message = FALSE)
  
  applySettingsDefault(list(startValue = startValues, burnin = 10000))
  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
  #summary(out)
  #plot(out)
  ddm <- rbind(ddm, c(mv, a, ter))
  idm <- rbind(idm,  as.numeric(MAP(out)$parametersMAP))
  
}
plot(ddm[,1], idm[,1], xlab = "Drift Rate", ylab = "Stimulus Distinctness") # v - c
plot(ddm[,2], idm[,2], xlab = "Boundary Seperation", ylab = "Detection Box Size") # h - a
plot(ddm[,3], idm[,3], xlab = "Non-decision Time", ylab = "Non-decision Time") # ter

cor(ddm[,1], idm[,1])
cor(ddm[,2], idm[,2])
cor(ddm[,3], idm[,3])

length(ddm[,1])
length(idm[,1])

result <- data.frame(cbind(ddm, idm))

colnames(result) <- c("v", "a", "ter", "c","h", "ter")

write.csv(result, "fit_results_EM_right_range")

## try DRAM function
#dram_result <- DRAM(startValue = as.vector(outDEoptim$member$bestmemit[10,]),
#                    iterations = 100, nBI = 100, parmin = c(0, 0, 0),
#                    parmax = c(1, 1, 1), FUN = f)

## plot IDM with different parameters
sampler <- function(n=1){
  d1 <- runif(n, 0,0.4)
  d2 <- runif(n, 0,0.6)
  d3 <- runif(n, 0.2,0.6)
  return(cbind(d1,d2,d3))
}
f(par)
pars <- c()
lls <- c()

for (i in 1:1000){
  par <- as.numeric(sampler())
  pars <- rbind(pars, par)
  lls <- c(lls, f(par))
  print(length(lls))
}

fig <- plot_ly(
  x = as.numeric(pars[,1]), 
  y = as.numeric(pars[,3]), 
  z = lls, 
  type = "contour",
  contours = list(start = min(lls), end = max(lls), size = 1000)
)
fig
as.numeric(pars[,1])

## IDM for speed accuracy trade-off?

## DDM connection between potential function and starting point




