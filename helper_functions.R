library(diffIRT)

free_energy <- function(y1,y2,N=2000,W_pos=52500,W_neg=8400,B1=0,B2=0,beta=1/24,theta=51450){
  sum_y <- -N/2*(y1*log(y1)+(1-y1)*log(1-y1) + y2*log(y2)+(1-y2)*log(1-y2)) 
  # no multiple lines!!!
  F <- -W_pos*y1^2 - W_pos*y2^2 + W_neg*y1*y2 - (B1-theta)*y1 - (B2-theta)*y2 - beta^(-1)*sum_y
  return(F)
}

## get the probability for all states 
N <- 1000 # here N is the number of neurons in one population
beta <- 1/24
prob <- matrix(0, nrow = N-1, ncol = N-1)
# first calculate Z
for(i in 1:N-1){
  for(j in 1:N-1){
    prob[i,j] <- exp(-beta*free_energy(i/N,j/N))
  }
}

Z <- sum(prob)
prob_scaled <- prob/Z # divided by Z to make it a valid pmf
prob_scaled_vec <- as.vector(prob_scaled)



dFdy1 <- function(y1,y2, W_pos, W_neg, B1, theta, beta, N){
  #print(log(y1/(1-y1)))
  result <- -W_pos*2*y1+W_neg*y2-(B1-theta)+beta^(-1)*0.5*N*(log(y1/(1-y1)))
  return(result)
}

dFdy2 <- function(y1,y2, W_pos, W_neg, B2, theta, beta, N){
  result <- -W_pos*2*y2+W_neg*y1-(B2-theta)+beta^(-1)*0.5*N*(log(y2/(1-y2)))
  return(result)
}


clamp <- function(n, minn=0, maxn=1){
  return(max(min(maxn, n), minn))
}

simul_IDM <- function(nsim, prob_scaled_vec, N = 2000, W_pos = 52500, W_neg = 8400, B1 = 3300, B2 = 3300, beta = 1/24, theta = 51450, D = 0.05, Ter = 0.3, h = 0.4){
  start_time <- Sys.time()
  result <- c()
  for(i in 1:nsim){
    sample_y <- sample_from_Py(N/2,prob_scaled_vec)
    y1 <- sample_y[1]
    y2 <- sample_y[2]
    dt <- 0.001
    iter <- 1
    wiener_pre1 <- 0
    wiener_pre2 <- 0
    location <- c(y1, y2)
    #print(((y1<0.6)&(y2<0.6))|((y1>0.4)&(y2>0.4)))
    while(((y1<(1-h))&(y2<(1-h)))|((y1>h)&(y2>h))){
      #dy1 <- -beta*D* dFdy1(y1,y2, W_pos, W_neg, B1, theta, beta, N) * dt + sqrt(2*D)*Wiener(n = 1, pts = 1, K = 1)
      wiener_post1 <- rnorm(1, 0, sqrt(iter*dt))
      wiener_post2 <- rnorm(1, 0, sqrt(iter*dt))
      dy1 <- -beta*D* dFdy1(y1,y2, W_pos, W_neg, B1, theta, beta, N) * dt + sqrt(2*D)*(wiener_post1 - wiener_pre1)
      dy2 <- -beta*D* dFdy2(y1,y2, W_pos, W_neg, B2, theta, beta, N) * dt + sqrt(2*D)*(wiener_post2 - wiener_pre2)
      
      wiener_pre1 <- wiener_post1
      wiener_pre2 <- wiener_post2
      
      y1 <- clamp(y1 + dy1) # clamp to [0,1]
      y2 <- clamp(y2 + dy2)
      location <- cbind(location, c(y1,y2))
      iter <- iter + 1
      
    if ((y1>(1-h))&(y2<h)){
      RT <- iter*dt + Ter
      if(C>=0){
        R <- 2 # correct response
      }else if(C<0){
        R <- 1 # incorrect response
      }
      
      result <- rbind(result, c(RT, R))
    }
    if ((y1<h)&(y2>(1-h))){
      RT <- iter*dt + Ter
      
      if(C>=0){
        R <- 1 # incorrect response
      }else if(C<0){
        R <- 2 # correct response
      }
      result <- rbind(result, c(RT, R))
    }
    }
  }
  result_df <- data.frame(result)
  colnames(result_df) <- c('RT','Response')
  print(Sys.time() - start_time)
  return(result_df)
}

sample_from_Py <- function(N,prob_scaled_vec){
  #print(length(prob_scaled_vec))
  #print(sum(prob_scaled_vec))
  index <- sample(x=1:((N-1)*(N-1)), size=1, replace = T, prob=prob_scaled_vec)
  quentient <- index%/%(N-1)
  remainder <- index%%(N-1)
  if (remainder != 0){
    original_index <- c(remainder, quentient+1)
  }else if(remainder == 0){
    original_index <- c((N-1), quentient)
  }
  return(original_index/N)
}


DDM_fit <- function(IDM_data){
  #start_time <- Sys.time()
  
  IDM_diffusion <- function(pars, rt, response){
    densities <- ddiffusion(rt, response=response, 
                            a=pars["a"],        # threshold separation.
                            v=pars["v"],        # drift rate.
                            t0=pars["t0"],      # non-decision time. 
                            #sz=pars["sz"],     # inter-trial-variability of starting point.
                            #st0=pars["st0"],   # inter-trial-variability of non-decision components.
                            sv=pars["sv"]       # inter-trial-variability of drift rate, 0<sv<2, default = 0.
                            )
    if (any(densities == 0)) return(1e6)
    return(-sum(log(densities)))
  }
  start <- c(runif(2, 0.5, 3), 0.1, runif(1, 0, 0.5))
  names(start) <- c("a", "v", "t0", "sv")#, "sz", "st0")
  IDM_diff <- nlminb(start, IDM_diffusion, lower = 0, 
                     rt=IDM_data$RT, response=IDM_data$Response)
  #print(Sys.time() - start_time)
  return(IDM_diff)
}
# simulate in IDM and fit with DDM
IDM_DDM <- function(C=0.2, Bns = 2800, Bs = 200, Ter = 0.3, h = 0.4){
  B1 <- Bs*(1+C) + Bns
  B2 <- Bs*(1-C) + Bns
  simul_results <- simul_IDM(1000, prob_scaled_vec = prob_scaled_vec, B1 = B1, B2 = B2, Ter = Ter, h = h)
  DDM_fit_results <- DDM_fit(simul_results)
  #DDM_fit_params <- DDM_fit_results$par
  #DDM_fit_convergence <- DDM_fit_results$convergence
  return(DDM_fit_results)
}

## simulation drift-diffusion model
simul_DDM_one_trial <- function(v=6, a=3, t0=0.3, c=1, sv=1){
  v <- rnorm(1, v, sv)
  ys <- c()
  y <- a/2
  dt <- 0.001
  iter <- 1
  wiener_pre <- 0
  
  while((y<a)&(y>0)){
    wiener_post <- rnorm(1, 0, sqrt(iter*dt))
    dy <- v* dt + c*(wiener_post - wiener_pre)
    print(c*(wiener_post - wiener_pre)/v*dt)
    
    wiener_pre <- wiener_post
    
    y <- y + dy 
    iter <- iter + 1
    ys <- c(ys, y)
    }
  return(ys)
}

#ys <- simul_DDM_one_trial()
#ys
#plot(ys, type = "l")


## many DDM trials
simul_DDM <- function(nsim, v=6, a=3, t0=0.3, c=1, sv=1){
  start_time <- Sys.time()
  result <- c()
  for(i in 1:nsim){
    v <- rnorm(1, v, sv)
    y <- a/2
    wiener_pre <- 0
    dt <- 0.001
    iter <- 1

    while((y<a)&(y>0)){
      wiener_post <- rnorm(1, 0, sqrt(iter*dt))
      dy <- v*dt + c*(wiener_post - wiener_pre)

      wiener_pre <- wiener_post

      y <- y + dy 
      iter <- iter + 1
      }
    if (y>=a){
      R <- 2 # correct response
    }
    if (y<=0){
      R <- 1
    }
    RT <- iter*dt + t0
    result <- rbind(result, c(RT, R))

  }
  result_df <- data.frame(result)
  #print(result)
  colnames(result_df) <- c('RT','Response')
  print(Sys.time() - start_time)
  return(result_df)
}

#simul_results_DDM <- simul_DDM(nsim = 1000, a=2, v=1, t0=0.3, sv=0.3)
#prop.table(table(simul_results_DDM$Response))

## simulate data accroding to the traditional diffusion model
# set.seed(1)


#data=simdiffT(N=100,a=2,mv=1,sv=0.3,ter=0.3)
#simdiff_results <- data.frame('RT' = data$rt, 'Response' = data$x+1)
#prop.table(table(simdiff_results$Response))
