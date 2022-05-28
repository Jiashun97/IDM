
#  test whether your RTDist installation is capable 
# of running simulations on the CPU
source('sourceAllRTDist.R')
testWrappers()


settings=list()
settings$nWalkers=500000
settings$nGPUs=0
settings$gpuIds=0
settings= checkSettings(settings)

#LDM
LDMPars=list()
LDMPars$a=0.08
LDMPars$nStimuli=2
LDMPars$v=matrix(c(0,0.1),1,2)
LDMPars= checkLDMPars(LDMPars)
out = LDMDist(LDMPars,settings) # "distributions" "extra" "ok"  "used"  
plot(matrix(1:3000,3000, 2*LDMPars$nStimuli), out$distributions,'l')

# IDM
IDMPars=list()
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
IDMPars$D <- matrix(c(0.05,0.05),2,1)
IDMPars$sigma <- matrix(c(0.02,0.02),2,1)
IDMPars$deltat <- 0.001
  


IDMPars= checkIDMPars(IDMPars)
out = IDMDist(IDMPars,settings) # "distributions" "extra" "ok"  "used"  
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

