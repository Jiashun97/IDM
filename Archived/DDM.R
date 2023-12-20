require(rtdists)
library(diffIRT)

# Exp. 1; Wagenmakers, Ratcliff, Gomez, & McKoon (2008, JML)
data(speed_acc)   
# remove excluded trials:
speed_acc <- droplevels(speed_acc[!speed_acc$censor,]) 
# create numeric response variable where 1 is an error and 2 a correct response: 
speed_acc$corr <- with(speed_acc, as.numeric(stim_cat == response))+1 
# select data from participant 11, accuracy condition, non-word trials only
p11 <- speed_acc[speed_acc$id == 11 & 
                   speed_acc$condition == "accuracy" & 
                   speed_acc$stim_cat == "nonword",] 
# proportion of correct responses
prop.table(table(p11$corr))


ll_diffusion <- function(pars, rt, response){
  densities <- ddiffusion(rt, response=response, 
                          a=pars["a"], 
                          v=pars["v"], 
                          t0=pars["t0"], 
                          sz=pars["sz"], 
                          st0=pars["st0"],
                          sv=pars["sv"])
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}
start <- c(runif(2, 0.5, 3), 0.1, runif(3, 0, 0.5))
names(start) <- c("a", "v", "t0", "sz", "st0", "sv")
p11_diff <- nlminb(start, ll_diffusion, lower = 0, 
                   rt=p11$rt, response=p11$corr)
p11_diff

# quantiles:
q <- c(0.1, 0.3, 0.5, 0.7, 0.9)

# observed data:
(p11_q_c <- quantile(p11[p11$corr == 2, "rt"], probs = q))
(p11_q_e <- quantile(p11[p11$corr == 1, "rt"], probs = q))

### diffusion:
# same result as when using Inf, but faster:
(pred_prop_correct_diffusion <- pdiffusion(rt = 20,  response = "upper",
                                           a=p11_diff$par["a"], 
                                           v=p11_diff$par["v"], 
                                           t0=p11_diff$par["t0"], 
                                           sz=p11_diff$par["sz"], 
                                           st0=p11_diff$par["st0"], 
                                           sv=p11_diff$par["sv"]))  

(pred_correct_diffusion <- qdiffusion(q*pred_prop_correct_diffusion, 
                                      response = "upper",
                                      a=p11_diff$par["a"], 
                                      v=p11_diff$par["v"], 
                                      t0=p11_diff$par["t0"], 
                                      sz=p11_diff$par["sz"], 
                                      st0=p11_diff$par["st0"], 
                                      sv=p11_diff$par["sv"]))
(pred_error_diffusion <- qdiffusion(q*(1-pred_prop_correct_diffusion), 
                                    response = "lower",
                                    a=p11_diff$par["a"], 
                                    v=p11_diff$par["v"], 
                                    t0=p11_diff$par["t0"], 
                                    sz=p11_diff$par["sz"], 
                                    st0=p11_diff$par["st0"], 
                                    sv=p11_diff$par["sv"]))

plot(p11_q_c, q*prop.table(table(p11$corr))[2], pch = 2, ylim=c(0, 1), xlim = c(0.4, 1.3), ylab = "Cumulative Probability", xlab = "Response Time (sec)", main = "Diffusion")
points(p11_q_e, q*prop.table(table(p11$corr))[1], pch = 2)
lines(pred_correct_diffusion, q*pred_prop_correct_diffusion, type = "b")
lines(pred_error_diffusion, q*(1-pred_prop_correct_diffusion), type = "b")

rand_rts <- rdiffusion(1e5, a=p11_diff$par["a"], 
                       v=p11_diff$par["v"], 
                       t0=p11_diff$par["t0"], 
                       sz=p11_diff$par["sz"], 
                       st0=p11_diff$par["st0"], 
                       sv=p11_diff$par["sv"])
plot(ecdf(rand_rts[rand_rts$response == "upper","rt"]))
normalised_pdiffusion = function(rt,...) pdiffusion(rt,...)/pdiffusion(rt=Inf,...) 
curve(normalised_pdiffusion(x, response = "upper",
                            a=p11_diff$par["a"], 
                            v=p11_diff$par["v"], 
                            t0=p11_diff$par["t0"], 
                            sz=p11_diff$par["sz"], 
                            st0=p11_diff$par["st0"], 
                            sv=p11_diff$par["sv"]), 
      add=TRUE, col = "yellow", lty = 2)




## Experiment 2. fit the model with simulated data from IDM
# proportion of correct responses

IDM_data <- simdiff_results
prop.table(table(IDM_data$Response))


IDM_diffusion <- function(pars, rt, response){
  densities <- ddiffusion(rt, response=response, 
                          a=pars["a"], 
                          v=pars["v"], 
                          t0=pars["t0"], 
                          #sz=pars["sz"], 
                          #st0=pars["st0"],
                          sv=pars["sv"])
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}
start <- c(runif(2, 0.5, 3), 0.1, runif(1, 0, 0.5))
names(start) <- c("a", "v", "t0", "sv") # , "sz", "st0"
IDM_diff <- nlminb(start, IDM_diffusion, lower = 0, 
                   rt=IDM_data$RT, 
                   response=IDM_data$Response)
IDM_diff # 2.294046e+00 2.926467e+00 3.218440e-14 1.836911e-01 2.993832e-01 

# quantiles:
q <- c(0.1, 0.3, 0.5, 0.7, 0.9)

# observed data:
(IDM_q_c <- quantile(IDM_data[IDM_data$Response == 2, "RT"], probs = q))
(IDM_q_e <- quantile(IDM_data[IDM_data$Response == 1, "RT"], probs = q))


### diffusion:
# same result as when using Inf, but faster:
(pred_prop_correct_diffusion <- pdiffusion(rt = 20,  response = "upper",
                                           a=IDM_diff$par["a"], 
                                           v=IDM_diff$par["v"], 
                                           t0=IDM_diff$par["t0"], 
                                           sz=IDM_diff$par["sz"], 
                                           st0=IDM_diff$par["st0"], 
                                           sv=IDM_diff$par["sv"]))  

(pred_correct_diffusion <- qdiffusion(q*pred_prop_correct_diffusion, 
                                      response = "upper",
                                      a=IDM_diff$par["a"], 
                                      v=IDM_diff$par["v"], 
                                      t0=IDM_diff$par["t0"], 
                                      sz=IDM_diff$par["sz"], 
                                      st0=IDM_diff$par["st0"], 
                                      sv=IDM_diff$par["sv"]))
(pred_error_diffusion <- qdiffusion(q*(1-pred_prop_correct_diffusion), 
                                    response = "lower",
                                    a=IDM_diff$par["a"], 
                                    v=IDM_diff$par["v"], 
                                    t0=IDM_diff$par["t0"], 
                                    sz=IDM_diff$par["sz"], 
                                    st0=IDM_diff$par["st0"], 
                                    sv=IDM_diff$par["sv"]))

plot(IDM_q_c, q*prop.table(table(IDM_data$Response))[2], pch = 2, ylab = "Cumulative Probability", xlab = "Response Time (sec)", main = "Diffusion")
lines(pred_correct_diffusion, q*pred_prop_correct_diffusion, type = "b")

points(IDM_q_e, q*prop.table(table(IDM_data$Response))[1], pch = 2, col = "red")
lines(pred_error_diffusion, q*(1-pred_prop_correct_diffusion), type = "b", col = "red")

legend("topright", legend = c("data(correct)", "predictions(correct)", "data(error)", "predictions(error)"), 
       pch = c(2,1,2,1), lty = c(0,1,0,1), col = c("black","black","red","red"))







rand_rts <- rdiffusion(1e5, a=IDM_diff$par["a"], 
                       v=IDM_diff$par["v"], 
                       t0=IDM_diff$par["t0"], 
                       sz=IDM_diff$par["sz"], 
                       st0=IDM_diff$par["st0"], 
                       sv=IDM_diff$par["sv"])
plot(ecdf(rand_rts[rand_rts$response == "upper","rt"]))
normalised_pdiffusion = function(rt,...) pdiffusion(rt,...)/pdiffusion(rt=Inf,...) 
curve(normalised_pdiffusion(x, response = "upper",
                            a=IDM_diff$par["a"], 
                            v=IDM_diff$par["v"], 
                            t0=IDM_diff$par["t0"], 
                            sz=IDM_diff$par["sz"], 
                            st0=IDM_diff$par["st0"], 
                            sv=IDM_diff$par["sv"]), 
      add=TRUE, col = "yellow", lty = 2)


#### recover parameters from data simulated from DDM
data=simdiffT(N=10000,a=2,mv=1,sv=0.3,ter=0.3)
simdiff_results <- data.frame('RT' = data$rt, 'Response' = data$x+1)
prop.table(table(simdiff_results$Response))


simdiffT_diffusion <- function(pars, rt, response){
  densities <- ddiffusion(rt, response=response, 
                          a=pars["a"], 
                          v=pars["v"], 
                          t0=pars["t0"], 
                          #sz=pars["sz"], 
                          #st0=pars["st0"],
                          sv=pars["sv"])
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}
start <- c(runif(2, 0.5, 3), 0.1, runif(1, 0, 0.5))
names(start) <- c("a", "v", "t0", "sv") # , "sz", "st0"
simdiffT_diff <- nlminb(start, simdiffT_diffusion, lower = 0, 
                   rt=simdiff_results$RT, 
                   response=simdiff_results$Response)
simdiffT_diff







