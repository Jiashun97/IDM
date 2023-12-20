require(rtdists)
library(diffIRT)
library(condMVNorm)
library(MASS)
library(pracma)
library(DEoptim)
library(BayesianTools)


energy <- function(y1,y2,W_pos=52500,W_neg=8400,B1=0,B2=0,theta=51450){
  # no multiple lines!!!
  E <- -W_pos*y1^2 - W_pos*y2^2 + W_neg*y1*y2 - (B1-theta)*y1 - (B2-theta)*y2
  return(E)
}

dEdy1 <- function(y1,y2, W_pos, W_neg, B1, theta){
  result <- -W_pos*2*y1+W_neg*y2-(B1-theta)
  return(result)
}

dEdy2 <- function(y1,y2, W_pos, W_neg, B2, theta){
  result <- -W_pos*2*y2+W_neg*y1-(B2-theta)
  return(result)
}

free_energy <- function(y1,y2,N=2000,W_pos=52500,W_neg=8400,B1=0,B2=0,beta=1/24,theta=51450){
  sum_y <- -N/2*(y1*log(y1)+(1-y1)*log(1-y1) + y2*log(y2)+(1-y2)*log(1-y2)) 
  # no multiple lines!!!
  F <- -W_pos*y1^2 - W_pos*y2^2 + W_neg*y1*y2 - (B1-theta)*y1 - (B2-theta)*y2 - beta^(-1)*sum_y
  return(F)
}

dFdy1 <- function(y1,y2, W_pos, W_neg, B1, theta, beta, N){
  #print(log(y1/(1-y1)))
  result <- -W_pos*2*y1+W_neg*y2-(B1-theta)+beta^(-1)*0.5*N*(log(y1/(1-y1)))
  return(result)
}

dFdy2 <- function(y1,y2, W_pos, W_neg, B2, theta, beta, N){
  result <- -W_pos*2*y2+W_neg*y1-(B2-theta)+beta^(-1)*0.5*N*(log(y2/(1-y2)))
  return(result)
}

## get the probability for all states 
#N <- 1000 # here N is the number of neurons in one population
#beta <- 1/24
#prob <- matrix(0, nrow = N-1, ncol = N-1)
# first calculate Z
#for(i in 1:N-1){
#  for(j in 1:N-1){
#    prob[i,j] <- exp(-beta*free_energy(i/N,j/N))
#  }
#}

#Z <- sum(prob)
#prob_scaled <- prob/Z # divided by Z to make it a valid pmf
#prob_scaled_vec <- as.vector(prob_scaled)





clamp <- function(n, minn=0.0001, maxn=0.9999){
  return(max(min(maxn, n), minn))
}

# sigma = 0.001, dt <- 0.00001 -> D = 0.05
# sigma = 0.001,  dt <- 0.001 -> D = 5e-04
# sigma = 0.01,  dt <- 0.001 -> D = 0.05
get_res <- function(y1,y2,h,B1,B2){
  if ((y1>(1-h))&(y2<h)){
    if(B1>=B2){
      R <- 2 # correct response
    }else if(B1<B2){
      R <- 1 # incorrect response
    }
  }
  
  if ((y1<h)&(y2>(1-h))){
    if(B1>=B2){
      R <- 1 # incorrect response
    }else if(B1<B2){
      R <- 2 # correct response
    }
  }
  return(R)
}

continue <- function(y1, y2, h){
  overlap <- ((y1>(1-h))&(y2>(1-h)))&((y1<h)&(y2<h)) # in overlapped region?
  if (overlap){
    return(TRUE)
  }else{
    return(((y1<(1-h))&(y2<(1-h)))|((y1>h)&(y2>h)))
  }
}

get_state <- function(y1,y2,h){
  overlap <- ((y1>(1-h))&(y2>(1-h)))&((y1<h)&(y2<h))
  if(overlap){
    return(0)
  }
  
  if((y1>(1-h))&(y2<h)){
    state <- 1
  }else if((y1<h)&(y2>(1-h))){
    state <- 2
  }else{
    state <- 0
  }
  return(state)
}

get_state_DDM <- function(y,a){
  if(y > a){
    state <- 1
  }else if(y < 0){
    state <- 2
  }else{
    state <- 0
  }
  return(state)
}

mind_change <- function(s1,s2){
  if((s1==1)&(s2==2)){
    return(TRUE)
  }
  
  if((s1==2)&(s2==1)){
    return(TRUE)
  }
  return(FALSE)
}

simul_IDM_non_stop <- function(C = 0.2, N = 2000, W_pos = 52500, W_neg = 8400, Bns = 2500, Bs = 1000, # mode = "coarse-grained"
                                beta = 1/24, theta = 51450, Ter = 0.3, h = 0.4, D = 0.05, sigma = 0.01, dt = 0.001, 
                                init_loc = c(0.3,0.3), n_iter = 100){ # D = 0.05, 2/(1000*exp(1)) 
  B1 <- Bs*(1+C) + Bns
  B2 <- Bs*(1-C) + Bns
  #start_time <- Sys.time() 
  y1 <- init_loc[1] #sample_y[1]
  y2 <- init_loc[2] #sample_y[2]
  iter <- 1
  location <- c(y1, y2)
    
  while(iter<=n_iter){
    dW1 <- rnorm(1)*sqrt(dt) 
    dW2 <- rnorm(1)*sqrt(dt)
    
    dy1 <- -beta*D*dFdy1(y1,y2, W_pos, W_neg, B1, theta, beta, N)*dt + sqrt(2*D)*dW1
    dy2 <- -beta*D*dFdy2(y1,y2, W_pos, W_neg, B2, theta, beta, N)*dt + sqrt(2*D)*dW2
    
    y1 <- clamp(y1 + dy1) # clamp to [0,1]
    y2 <- clamp(y2 + dy2)
    location <- cbind(location, c(y1,y2))
    iter <- iter + 1
  }
  state <- get_state(y1,y2,h)
  return(list(loc = c(y1,y2), state = state))
}



simul_DDM_non_stop <- function(v, dt = 0.001, a = 0.8, init_loc = 0.4, n_iter = 2000, c = sqrt(0.1)){

  iter <- 1
  location <- init_loc
  y <- init_loc
  
  while(iter<=n_iter){
    dW <- rnorm(1)*sqrt(dt) 
    dy <- v*dt + c*dW
    y <- y + dy
    
    location <- c(location, y)
    iter <- iter + 1
  }
  state <- get_state_DDM(y,a)
  return(list(loc = y, state = state))
}


simul_DDM <- function(v, dt = 0.001, a = 0.8, init_loc = 0.4, c = sqrt(0.1)){
  
  iter <- 1
  location <- init_loc
  y <- init_loc
  
  while(get_state_DDM(y,a)==0){
    dW <- rnorm(1)*sqrt(dt) 
    dy <- v*dt + c*dW
    y <- y + dy
    
    location <- c(location, y)
    iter <- iter + 1
  }
  
  state <- get_state_DDM(y,a)
  return(list(loc = y, state = state, iter = iter, traj = location))
}

simul_DDM_post <- function(v, dt = 0.001, a = 0.8, post_iter = 500, init_loc = 0.4, c = sqrt(0.1)){
  
  iter <- 1
  location <- init_loc
  y <- init_loc
  
  while(iter <= post_iter){
    dW <- rnorm(1)*sqrt(dt) 
    dy <- v*dt + c*dW
    y <- y + dy
    
    location <- c(location, y)
    iter <- iter + 1
  }
  
  state <- get_state_DDM(y,a)
  return(list(loc = y, state = state))
}
# y>a => 1, y<0 => 2, otherwise 0
simul_DDM_com <- function(v, dt = 0.001, a = 0.8, post_iter = 500, init_loc = 0.4, c = sqrt(0.1)){
  init_result <- simul_DDM(v = v, a = a, init_loc = init_loc)
  final_result <- simul_DDM_post(v = v, a = a, post_iter = post_iter, init_loc = init_result$loc)
  
  uncertain <- as.numeric(final_result$state == 0) # whether final state is in the middle of two boundaries
  
  if((init_result$state==1)&(final_result$state==2)){ # whether decision changes
    change_of_mind <- 1
  }else if((init_result$state==2)&(final_result$state==1)){
    change_of_mind <- 1
  }else{
    change_of_mind <- 0
  }
  
  return(list(uncertain = uncertain, change = change_of_mind, iter = init_result$iter))
  
}


simul_IDM_one_trial <- function(C = 0.2, N = 2000, W_pos = 52500, W_neg = 8400, Bns = 2500, Bs = 1000, # mode = "coarse-grained"
                      beta = 1/24, theta = 51450, Ter = 0.3, h = 0.4, D = 0.05, sigma = 0.01, dt = 0.001, 
                      init_loc = c(0.3,0.3), mode = "coarse-grained"){ # D = 0.05, 2/(1000*exp(1))
  B1 <- Bs*(1+C) + Bns
  B2 <- Bs*(1-C) + Bns
  y1 <- init_loc[1] #sample_y[1]
  y2 <- init_loc[2] #sample_y[2]
  iter <- 1
  location <- c(y1, y2)
  if(mode == "Euler-Maruyama"){
    while(continue(y1=y1, y2=y2, h=h)){
      dW1 <- rnorm(1)*sqrt(dt) 
      dW2 <- rnorm(1)*sqrt(dt)
      
      dy1 <- -beta*D*dFdy1(y1,y2, W_pos, W_neg, B1, theta, beta, N)*dt + sqrt(2*D)*dW1#(wiener_post1 - wiener_pre1)
      dy2 <- -beta*D*dFdy2(y1,y2, W_pos, W_neg, B2, theta, beta, N)*dt + sqrt(2*D)*dW2#(wiener_post2 - wiener_pre2)
      
      y1 <- clamp(y1 + dy1) # clamp to [0,1]
      y2 <- clamp(y2 + dy2)
      location <- cbind(location, c(y1,y2))
      iter <- iter + 1
    }
  }else if(mode == "coarse-grained"){
    while(continue(y1=y1, y2=y2, h=h)){
      y_step <- mvrnorm(n = 1, mu=c(0,0), Sigma=matrix(c(sigma^2,0,0,sigma^2), 2,2))
      y_candidate1 <- clamp(y1+y_step[1])
      y_candidate2 <- clamp(y2+y_step[2])
      #print(y_step)
      df <- free_energy(y_candidate1,y_candidate2,B1=B1,B2=B2) - free_energy(y1,y2,B1=B1,B2=B2)
      # Matropolis acceptance rule
      if (df<=0){ # accept the move
        y1 <- y_candidate1 # clamp to [0,1]
        y2 <- y_candidate2
      }else{ # accept the move with a probability
        p_accept <- exp(-beta*df)
        if (p_accept >= runif(1)){
          y1 <- y_candidate1 # clamp to [0,1]
          y2 <- y_candidate2
          location <- cbind(location, c(y1,y2))
          }
        }
      iter <- iter + 1
    }
  }else{
    print("Wrong mode!")
  }
  RT <- clamp(iter*dt + Ter, 0, 3)
  R <- get_state(y1 = y1, y2 = y2, h = h)
  
  return(list("Trajectory" = location, "RT" = RT, "R" = R, "loc" = c(y1,y2), "iter" = iter))
}



simul_IDM_post_trial <- function(duration, init_loc = c(0.3,0.3), C = 0.2, N = 2000, W_pos = 52500, W_neg = 8400, 
                                 Bns = 2500, Bs = 1000, beta = 1/24, theta = 51450, Ter = 0.3, h = 0.4, 
                                 D = 0.05, sigma = 0.01, dt = 0.001, mode = "coarse-grained"){ # D = 0.05, 2/(1000*exp(1))
  B1 <- Bs*(1+C) + Bns
  B2 <- Bs*(1-C) + Bns
  #start_time <- Sys.time() 
  y1 <- init_loc[1]
  y2 <- init_loc[2]
  iter <- 1
  location <- c(y1, y2)
  if(mode == "Euler-Maruyama"){
    
    while(iter<=duration){
      dW1 <- rnorm(1)*sqrt(dt) 
      dW2 <- rnorm(1)*sqrt(dt)
      
      dy1 <- -beta*D*dFdy1(y1,y2, W_pos, W_neg, B1, theta, beta, N)*dt + sqrt(2*D)*dW1#(wiener_post1 - wiener_pre1)
      dy2 <- -beta*D*dFdy2(y1,y2, W_pos, W_neg, B2, theta, beta, N)*dt + sqrt(2*D)*dW2#(wiener_post2 - wiener_pre2)
      
      y1 <- clamp(y1 + dy1) # clamp to [0,1]
      y2 <- clamp(y2 + dy2)
      location <- cbind(location, c(y1,y2))
     # if(iter>100){
     #    return(list("Trajectory" = location))
     # }
      iter <- iter + 1
    }
  }else if(mode == "coarse-grained"){
    loc_com <- NA
    iter_post <- NA
    record <- TRUE
    s1 <- get_state(y1 = y1, y2 = y2, h = h)
    while(iter<=duration){
      y_step <- mvrnorm(n = 1, mu=c(0,0), Sigma=matrix(c(sigma^2,0,0,sigma^2), 2,2))
      y_candidate1 <- clamp(y1+y_step[1])
      y_candidate2 <- clamp(y2+y_step[2])
      #print(y_step)
      df <- free_energy(y_candidate1,y_candidate2,B1=B1,B2=B2) - free_energy(y1,y2,B1=B1,B2=B2)
      # Matropolis acceptance rule
      if (df<=0){ # accept the move
        y1 <- y_candidate1 # clamp to [0,1]
        y2 <- y_candidate2
      }else{ # accept the move with a probability
        p_accept <- exp(-beta*df)
        if (p_accept >= runif(1)){
          y1 <- y_candidate1 # clamp to [0,1]
          y2 <- y_candidate2
          location <- cbind(location, c(y1,y2))
          
        }
      }
      s2 <- get_state(y1 = y1, y2 = y2, h = h)
      if (mind_change(s1,s2)&record){
        #print(s1)
        #print(s2)
        loc_com <- c(y1,y2)
        iter_post <- iter
        record <- FALSE
      }
      iter <- iter + 1
    }
    #print(mind_change(s1,s2))
    #print(s1)
    #print(s2)
  }else{
    print("Wrong mode!")
  }
  R <- get_state(y1 = y1, y2 = y2, h = h)
  #print(c(y1,y2))
  return(list("Trajectory" = location, "loc" = c(y1,y2), "R" = R, "loc_com" = loc_com, iter_post = iter_post))
}


simul_IDM_com <- function(C = 0.1, h = 0.4, duration = 500){
  init_result <- simul_IDM_one_trial(C = C, h = h, mode = "coarse-grained")
  final_result <- simul_IDM_post_trial(duration = duration, h = h, init_loc = init_result$loc, mode = "coarse-grained")
  
  uncertain <- as.numeric(final_result$R == 0) # whether final state is in the middle of two boundaries
  
  if((init_result$R==1)&(final_result$R==2)){ # whether decision changes
    change_of_mind <- 1
  }else if((init_result$R==2)&(final_result$R==1)){
    change_of_mind <- 1
  }else{
    change_of_mind <- 0
  }
  
  return(list(RT = init_result$RT, uncertain = uncertain, change = change_of_mind, 
              iter = init_result$iter, loc_init = init_result$loc, loc_com = final_result$loc_com, iter_post = final_result$iter_post))
  
}

simul_continuous <- function(stim_t = 5000, rest_t = 2000, init_loc = c(0.3,0.3), C = 0.2, N = 2000, W_pos = 52500, W_neg = 8400, 
                             Bns = 2500, Bs = 1000, beta = 1/24, theta = 51450, Ter = 0.3, h = 0.4, 
                             D = 0.05, sigma = 0.01, dt = 0.001, mode = "Euler-Maruyama"){
  # stimulus onset until decision made
  result <- simul_IDM_one_trial(C = C, W_pos = W_pos)
  loc1 <- result$Trajectory
  #print(loc1)
  # stimulus onset after decision made
  loc2 <- simul_IDM_post_trial(duration = (stim_t-result$RT*1000), init_loc = loc1[,dim(loc1)[2]], C = C, W_pos = W_pos)$Trajectory

  # stimulus offset
  loc3 <- simul_IDM_post_trial(duration = rest_t, init_loc = loc2[,dim(loc2)[2]], Bns = 0, Bs = 0, W_pos = W_pos)$Trajectory
  
  locations <- cbind(loc1, loc2, loc3)
  
  return(list("Trajectory" = locations, "loc1" = loc1, "loc2" = loc2, "loc3" = loc3, "RT" = result$RT, "R" = result$R))
}

simul_continuous_multiple <- function(simulations = 100, stim_t = 5000, rest_t = 2000, init_loc = c(0.3,0.3), 
                                      C = 0.2, N = 2000, W_pos = 52500, W_neg = 8400, 
                             Bns = 2500, Bs = 1000, beta = 1/24, theta = 51450, Ter = 0.3, h = 0.4, 
                             D = 0.05, sigma = 0.01, dt = 0.001, mode = "Euler-Maruyama"){
  results <- c()
  loc3 <- c(0.7,0.3)
  #simulations <- 3
  #W_pos = 51500
  for(i in 1:simulations){
    # check whether last location still in the region
    if (!continue(loc3[1],loc3[2],h)){
      cat("simulation", i, "same as before, skipped")
      loc2 <- simul_IDM_post_trial(duration = (stim_t-result$RT*1000), init_loc = loc3, C = C, W_pos = W_pos)$Trajectory
      loc2 <- loc2[,dim(loc2)[2]]
      # stimulus offset
      loc3 <- simul_IDM_post_trial(duration = rest_t, init_loc = loc2, Bns = 0, Bs = 0, W_pos = W_pos)$Trajectory
      loc3 <- loc3[,dim(loc3)[2]]
      
      next
    }
    # stimulus onset until decision made
    result <- simul_IDM_one_trial(C = C, W_pos = W_pos, init_loc = loc3)
    loc1 <- result$Trajectory
    loc1 <- loc1[,dim(loc1)[2]]
    #print(loc1)
    # stimulus onset after decision made
    loc2 <- simul_IDM_post_trial(duration = (stim_t-result$RT*1000), init_loc = loc1, C = C, W_pos = W_pos)$Trajectory
    loc2 <- loc2[,dim(loc2)[2]]
    # stimulus offset
    loc3 <- simul_IDM_post_trial(duration = rest_t, init_loc = loc2, Bns = 0, Bs = 0, W_pos = W_pos)$Trajectory
    loc3 <- loc3[,dim(loc3)[2]]

    results <- rbind(results, c(result$RT, result$R))

  }
  results_df <- data.frame(results)
  colnames(results_df) <- c('RT','Response')
  return(results_df)
}

# result <- simul_continuous(W_pos = W_pos)
# plot_trajectory(result, W_pos)

plot_trajectory <- function(result, W_pos){
  plot(result$loc1[1,], result$loc1[2,], 
       type = "l", xlim = c(0,1), ylim = c(0,1),
       xlab = "y1", ylab = "y2")
  
  lines(result$loc2[1,], result$loc2[2,], type = "l",col = "red")
  lines(result$loc3[1,], result$loc3[2,], type = "l",col = "green")
  
  
  lines(rep(0.6,11), seq(0,0.4,0.04), type = "l")
  lines(rep(0.4,11), seq(0.6,1,0.04), type = "l")
  
  lines(seq(0,0.4,0.04), rep(0.6,11), type = "l")
  lines(seq(0.6,1,0.04), rep(0.4,11), type = "l")
  
  legend(0.45, 1, legend=c("Stim. on Pre-resp.", "Stim. on Post-resp.", "Stim. off"),
         col=c("black", "red", "green"), lty=1:1, cex=0.8)
  title(main = paste("Self-excitation is", W_pos))
}

simul_IDM <- function(nsim, prob_scaled_vec, C, N = 2000, W_pos = 52500, W_neg = 8400, Bns = 2500, Bs = 500, # Bns = 1500, Bs = 1000
                      beta = 1/24, theta = 51450, Ter = 0.3, h = 0.4, D = 0.05, sigma = 0.01, dt = 0.001, mode = "Euler-Maruyama"){ # D = 0.05, 2/(1000*exp(1))
  
  B1 <- Bs*(1+C) + Bns
  B2 <- Bs*(1-C) + Bns
  #start_time <- Sys.time() 
  result <- c()
  iters <- c()
  for(i in 1:nsim){
    #sample_y <- sample_from_Py(N/2,prob_scaled_vec)
    y1 <- 0.3 #sample_y[1]
    y2 <- 0.3 #sample_y[2]
    iter <- 1
    #location <- c(y1, y2)
    if(mode == "Euler-Maruyama"){
      
      while(continue(y1=y1, y2=y2, h=h)){
        dW1 <- rnorm(1)*sqrt(dt) 
        dW2 <- rnorm(1)*sqrt(dt)
        
        dy1 <- -beta*D*dFdy1(y1,y2, W_pos, W_neg, B1, theta, beta, N)*dt + sqrt(2*D)*dW1#(wiener_post1 - wiener_pre1)
        dy2 <- -beta*D*dFdy2(y1,y2, W_pos, W_neg, B2, theta, beta, N)*dt + sqrt(2*D)*dW2#(wiener_post2 - wiener_pre2)
        
        y1 <- clamp(y1 + dy1) # clamp to [0,1]
        y2 <- clamp(y2 + dy2)
        #location <- cbind(location, c(y1,y2))
        iter <- iter + 1
      }
    }else if(mode == "coarse-grained"){
      while(continue(y1=y1, y2=y2, h=h)){
        y_step <- mvrnorm(n = 1, mu=c(0,0), Sigma=matrix(c(sigma^2,0,0,sigma^2), 2,2))
        y_candidate1 <- clamp(y1+y_step[1])
        y_candidate2 <- clamp(y2+y_step[2])
        #print(y_step)
        df <- free_energy(y_candidate1,y_candidate2,B1=B1,B2=B2) - free_energy(y1,y2,B1=B1,B2=B2)
        # Matropolis acceptance rule
        if (df<=0){ # accept the move
          y1 <- y_candidate1 # clamp to [0,1]
          y2 <- y_candidate2
        }else{ # accept the move with a probability
          p_accept <- exp(-beta*df)
          if (p_accept >= runif(1)){
            y1 <- y_candidate1 # clamp to [0,1]
            y2 <- y_candidate2
          }
          
        }
        #location <- cbind(location, c(y1,y2))
        
        iter <- iter + 1
      }
    }else{
      print("Wrong mode!")
    }
    RT <- clamp(iter*dt + Ter, 0, 3)
    R <- get_res(y1 = y1, y2 = y2, h = h, B1 = B1, B2 = B2)
    result <- rbind(result, c(RT, R))
    iters <- c(iters, iter)
    #print(i)
    #print(iter)
    if ((i%%1000) == 0){
      cat("simulation ",i," has finished. ")
      cat(iter, "iterations.")
    }
  }
  result_df <- data.frame(result)
  colnames(result_df) <- c('RT','Response')
  #print(Sys.time() - start_time)
  cat("Averaged teration is ", mean(iters), ". ")
  
  return(result_df)
}

#data_RT <- result_df

#data_RT <- simul_IDM(10, prob_scaled_vec = prob_scaled_vec,
                     #B1 = 3100, B2 = 2900, Ter = 0.3, h = 0.4, dt = 0.001)

#prop.table(table(data_RT$Response))

#log_ls <- c()
#for (C in seq(0.1,1,0.1)){
#  log_l <- IDM_likelihood(data_RT, C = C)
#  log_ls <- c(log_ls, log_l)
#}
#plot(seq(0.1,1,0.1), log_ls, xlab = "C", ylab = "log likelihood")


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

### simulate rt distribution from IDM
get_pdf_RT <- function(C = 0.5, h=0.4, Ter=0.3, Bns = 2500, Bs = 1000, mode = 1){ # 0 continuous, 1 coarse-grain
  B1 <- Bs*(1+C) + Bns
  B2 <- Bs*(1-C) + Bns
  #print(B1)
  #print(B2)
  print(C)
  print(h)
  print(Ter)
  
  settings<-list()
  settings$nWalkers<-100000 # number of simulations
  settings$nGPUs<-0
  settings$gpuIds<-0
  settings$seed <- 0
  settings$dt <- 0.001
  settings$nBins <- 3000
  settings$binWidth <- 0.001
  #settings= checkSettings(settings)
  
  
  IDMPars<-list()
  IDMPars$startPos <- matrix(c(0.3,0.3),2,1)
  IDMPars$beta <- 1/24
  IDMPars$boxShape <- 0
  IDMPars$nDim <- 2
  IDMPars$nStimuli <- 1
  IDMPars$N <- matrix(c(1000,1000),2,1)
  IDMPars$W <- matrix(c(52500,8400,8400, 52500),2,2)
  IDMPars$h <- matrix(c(h,h,h,h),2,2)
  IDMPars$Theta <- matrix(c(51450,51450),2,1)
  IDMPars$B <- matrix(c(B1,B2),2,1) # depends on number of stimuli 
  IDMPars$Ter <- Ter
  IDMPars$sTer <- 0
  IDMPars$spontaneousTime <- 0
  IDMPars$dynamics <- mode
  IDMPars$D <- matrix(c(0.05,0.05),2,1)
  IDMPars$sigma <- matrix(c(0.01,0.01),2,1)
  IDMPars$deltat <- 0.001
  IDMPars$profile <- 1
  # check
  #IDMPars= checkIDMPars(IDMPars)
  out = IDMDist(IDMPars,settings) # "distributions" "extra" "ok"  "used" 
  return(out$distributions)
}


pdf_cdf <- function(pdf){
  cdf <- c(pdf[1])
  for (i in 2:3000){
    cdf <- c(cdf, cdf[i-1] + pdf[i])
    #print(cdf)
  }
  return (cdf/sum(pdf))
}

cdf_pdf <- function(cdf){
  pdf <- c(cdf[1])
  for (i in 2:3000){
    pdf <- c(pdf, cdf[i] - cdf[i-1])
    #print(cdf)
  }
  return (pdf/sum(pdf))
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
  start <- c(runif(2, 0, 3), 0.1, runif(1, 0, 0.5))
  names(start) <- c("a", "v", "t0", "sv")#, "sz", "st0")
  IDM_diff <- nlminb(start, IDM_diffusion, lower = 0, 
                     rt=IDM_data$RT, response=IDM_data$Response)
  #print(Sys.time() - start_time)
  return(IDM_diff)
}
# simulate in IDM and fit with DDM
IDM_DDM <- function(C = 0.2, Bns = 2500, Bs = 1000, Ter = 0.3, h = 0.4, n_trials = 1000, D=0.05, mode = "coarse-grained"){
  simul_results <- simul_IDM(n_trials, prob_scaled_vec = prob_scaled_vec, C = C, Ter = Ter, h = h, D = D, mode = mode)
  DDM_fit_results <- DDM_fit(simul_results)
  #DDM_fit_params <- DDM_fit_results$par
  #DDM_fit_convergence <- DDM_fit_results$convergence
  return(DDM_fit_results)
}

## simulation drift-diffusion model
#simul_DDM_one_trial <- function(v=6, a=3, t0=0.3, c=1, sv=1){
#  v <- rnorm(1, v, sv)
#  ys <- c()
#  y <- a/2
#  dt <- 0.001
#  iter <- 1
#  wiener_pre <- 0
  
#  while((y<a)&(y>0)){
#    wiener_post <- rnorm(1, 0, sqrt(iter*dt))
#    dy <- v* dt + c*(wiener_post - wiener_pre)
#    print(c*(wiener_post - wiener_pre)/v*dt)
    
#    wiener_pre <- wiener_post
    
#    y <- y + dy 
#    iter <- iter + 1
#    ys <- c(ys, y)
#    }
#  return(ys)
#}

#ys <- simul_DDM_one_trial()
#ys
#plot(ys, type = "l")


## many DDM trials
#simul_DDM <- function(nsim, v=6, a=3, t0=0.3, c=1, sv=1){
#  start_time <- Sys.time()
#  result <- c()
#  for(i in 1:nsim){
#    v <- rnorm(1, v, sv)
#    y <- a/2
#    wiener_pre <- 0
#    dt <- 0.001
#    iter <- 1
#
#    while((y<a)&(y>0)){
#      print(iter)
#      wiener_post <- rnorm(1, 0, sqrt(iter*dt))
#      dy <- v*dt + c*(wiener_post - wiener_pre)#

#      wiener_pre <- wiener_post

#      y <- y + dy 
#      iter <- iter + 1
#      }
#    if (y>=a){
#      R <- 2 # correct response
#    }
#    if (y<=0){
#      R <- 1
#    }
#    RT <- iter*dt + t0
#    result <- rbind(result, c(RT, R))

#  }
#  result_df <- data.frame(result)
#  #print(result)
#  colnames(result_df) <- c('RT','Response')
#  #print(Sys.time() - start_time)
#  return(result_df)
#}

#simul_results_DDM <- simul_DDM(nsim = 1000, a=2, v=1, t0=0.3, sv=0.3)
#prop.table(table(simul_results_DDM$Response))

## simulate data accroding to the traditional diffusion model
# set.seed(1)




IDM_fit <- function(IDM_data){
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
  start <- c(runif(1, 0.5, 2),runif(1, 0, 4), runif(1, 0, 1), runif(1, 0, 0.5))
  names(start) <- c("a", "v", "t0", "sv")#, "sz", "st0")
  IDM_diff <- nlminb(start, IDM_diffusion, lower = 0, 
                     rt=IDM_data$RT, response=IDM_data$Response)
  #print(Sys.time() - start_time)
  return(IDM_diff)
}

#pdf_RT <- get_pdf_RT(C = 0.5, h=0.4, Ter=0.3, Bns = 2800, Bs = 200)

IDM_likelihood <- function(data_RT, pdf_RT){
  
  probs <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  q <- as.numeric(quantile(data_RT[data_RT$Response == 2, "RT"], probs = probs))*1000 # s to ms
  n <- hist(data_RT[data_RT$Response == 2, "RT"], br = c(0, q/1000, 10), plot = FALSE)$counts
  m <- hist(data_RT[data_RT$Response == 1, "RT"], br = c(0, q/1000, 10), plot = FALSE)$counts

  pdf_c <- (pdf_RT[,1]+0.01)/sum(pdf_RT[,1]+0.01) # make it numerically more stable
  pdf_e <- (pdf_RT[,2]+0.01)/sum(pdf_RT[,2]+0.01)
    
  cdf_c <- pdf_cdf(pdf_c)
  cdf_e <- pdf_cdf(pdf_e)
  
  # correct
  log_p1_c <- (n[1]-1)*log(cdf_c[q[1]])
  log_p2_c <- log(pdf_c[q[1]]) + (n[2]-1)*log(cdf_c[q[2]] - cdf_c[q[1]])
  log_p3_c <- log(pdf_c[q[2]]) + (n[3]-1)*log(cdf_c[q[3]] - cdf_c[q[2]])
  log_p4_c <- log(pdf_c[q[3]]) + (n[4]-1)*log(cdf_c[q[4]] - cdf_c[q[3]])
  log_p5_c <- log(pdf_c[q[4]]) + (n[5]-1)*log(cdf_c[q[5]] - cdf_c[q[4]])
  log_p6_c <- log(pdf_c[q[5]]) + (n[6])*log(1 - cdf_c[q[5]])
  log_c <- log_p1_c + log_p2_c + log_p3_c + log_p4_c + log_p5_c + log_p6_c
  
  # error
  log_p1_e <- m[1]*log(cdf_e[q[1]])
  log_p2_e <- m[2]*log(cdf_e[q[2]] - cdf_e[q[1]])
  log_p3_e <- m[3]*log(cdf_e[q[3]] - cdf_e[q[2]])
  log_p4_e <- m[4]*log(cdf_e[q[4]] - cdf_e[q[3]])
  log_p5_e <- m[5]*log(1 - cdf_e[q[5]])
  log_e <- log_p1_e + log_p2_e + log_p3_e + log_p4_e + log_p5_e
  
  # log likelihood
  L <- (log_c + log_e)
  return(L)
}

pdf_sample <- function(pdf){
  sample <- c()
  for (i in 1:3000){
    sample <- c(sample, rep(i, (pdf[i])))
  }
  return(sample)
}

##### the other way around


#log_ls <- c()
#for (C in seq(0,1,0.1)){
#  print(C)
#  data_RT <- simul_IDM(1000, prob_scaled_vec = prob_scaled_vec,
#                       B1 = 200*(1+C) + 2800, B2 = 200*(1-C) + 2800, Ter = 0.3, h = 0.4)
#  log_l <- IDM_likelihood(data_RT)
#  log_ls <- c(log_ls, log_l)
#}

#plot(seq(0,1,0.1), log_ls, xlab = "C", ylab = "log likelihood")





#rt_dist <- get_pdf_RT(C = 0, h=0.4, Ter=0.3, Bns = 2500, Bs = 1000)

get_chisquare_stats <- function(data_RT, rt_dist){
  q <-  c(0.1, 0.3, 0.5, 0.7, 0.9)
  p <-  c(.1,.2,.2,.2,.2,.1)
  nsim_c <- length(data_RT[data_RT$Response == 2, "RT"])
  nsim_e <- length(data_RT[data_RT$Response == 1, "RT"])
  # empirical RT quantile
  data_RT_q_c <- as.numeric(quantile(data_RT[data_RT$Response == 2, "RT"], probs = q))
  data_RT_q_e <- as.numeric(quantile(data_RT[data_RT$Response == 1, "RT"], probs = q))
  
  cdf_c <- pdf_cdf(rt_dist[,1])
  cdf_e <- pdf_cdf(rt_dist[,2])
  
  # expected cumulative probability
  expected_cp_c <- cdf_c[data_RT_q_c*1000] # from s to index of ms
  expected_cp_e <- cdf_e[data_RT_q_e*1000]
  
  expected_cp_c <- c(0, expected_cp_c, 1)
  expected_cp_e <- c(0, expected_cp_e, 1)
  
  expected_freq_c <- numeric(6)
  expected_freq_e <- numeric(6)
  for (i in 1:6){
    expected_freq_c[i] <- expected_cp_c[i+1] - expected_cp_c[i]
    expected_freq_e[i] <- expected_cp_e[i+1] - expected_cp_e[i]
  }
  expected_freq_c <- expected_freq_c*nsim_c
  expected_freq_e <- expected_freq_e*nsim_e
  if (prod(expected_freq_c) == 0){
    expected_freq_c <- expected_freq_c + rep(5,6)
    print("There is some bin containing 0 in correct responses")
  } 
  
  if(prod(expected_freq_e) == 0){
    expected_freq_e <- expected_freq_e + rep(5,6)
    print("There is some bin containing 0 in error responses")
  }
  chisquare_c <- sum((expected_freq_c-p*nsim_c)^2/expected_freq_c)
  chisquare_e <- sum((expected_freq_e-p*nsim_e)^2/expected_freq_e)
  
  return(log(chisquare_c + chisquare_e))
}
#simul_data <- get_pdf_RT(C = 0.2, h=0.4, Ter=0.8, Bns = 2500, Bs = 1000)


#rt_dist <- get_pdf_RT(C = 0, h=0.4, Ter=0.3, Bns = 2500, Bs = 1000)

get_chisquare_stats <- function(data_RT, rt_dist){
  q <-  c(0.1, 0.3, 0.5, 0.7, 0.9)
  p <-  c(.1,.2,.2,.2,.2,.1)
  nsim_c <- length(data_RT[data_RT$Response == 2, "RT"])
  nsim_e <- length(data_RT[data_RT$Response == 1, "RT"])
  
  # empirical RT quantile
  data_RT_q_c <- as.numeric(quantile(data_RT[data_RT$Response == 2, "RT"], probs = q))
  data_RT_q_e <- as.numeric(quantile(data_RT[data_RT$Response == 1, "RT"], probs = q))
  pdf_c <- rt_dist[,1]
  pdf_e <- rt_dist[,2]
  data_RT_q_c <- c(-1, data_RT_q_c*1000)
  data_RT_q_e <- c(-1, data_RT_q_e*1000)
  expected_freq_c <- numeric(6)
  expected_freq_e <- numeric(6)
  for (i in 1:5){ 
    #print(data_RT_q_c[i]+1)
    #print(data_RT_q_c[i+1])
    expected_freq_c[i] <- sum(pdf_c[(data_RT_q_c[i]+1):data_RT_q_c[i+1]])/sum(rt_dist[,1])*nsim_c
    expected_freq_e[i] <- sum(pdf_e[(data_RT_q_e[i]+1):data_RT_q_e[i+1]])/sum(rt_dist[,2])*nsim_e
  }
  expected_freq_c[6] <- sum(pdf_c[(data_RT_q_c[6]+1):3000])/sum(rt_dist[,1])*nsim_c
  expected_freq_e[6] <- sum(pdf_e[(data_RT_q_e[6]+1):3000])/sum(rt_dist[,2])*nsim_e
  
  if (prod(expected_freq_c) == 0){
    expected_freq_c <- expected_freq_c + rep(5,6)
    print("There is some bin containing 0 in correct responses")
  } 
  
  if(prod(expected_freq_e) == 0){
    expected_freq_e <- expected_freq_e + rep(5,6)
    print("There is some bin containing 0 in error responses")
  }
  chisquare_c <- sum((expected_freq_c-p*nsim_c)^2/expected_freq_c)
  chisquare_e <- sum((expected_freq_e-p*nsim_e)^2/expected_freq_e)
  
  return(log(chisquare_c + chisquare_e))
}

#rt_dist <- simul_data

get_chisquare_stats_IDM <- function(data_RT, rt_dist){
  q = c(0.1, 0.3, 0.5, 0.7, 0.9)
  p = c(.1,.2,.2,.2,.2,.1)
  # number of correct and incorrect responses
  nsim_c <- sum(data_RT[,1])
  nsim_e <- sum(data_RT[,2])
  # empirical RT quantile
  data_RT_q_c <- as.numeric(quantile(pdf_sample(data_RT[,1]), probs = q))
  data_RT_q_e <- as.numeric(quantile(pdf_sample(data_RT[,2]), probs = q))
  
  cdf_c <- pdf_cdf(rt_dist[,1])
  cdf_e <- pdf_cdf(rt_dist[,2])
  
  # expected cumulative probability
  expected_cp_c <- cdf_c[data_RT_q_c]
  expected_cp_e <- cdf_e[data_RT_q_e]
  
  expected_cp_c <- c(0, expected_cp_c, 1)
  expected_cp_e <- c(0, expected_cp_e, 1)
  
  expected_freq_c <- numeric(6)
  expected_freq_e <- numeric(6)
  for (i in 1:6){
    expected_freq_c[i] <- expected_cp_c[i+1] - expected_cp_c[i]
    expected_freq_e[i] <- expected_cp_e[i+1] - expected_cp_e[i]
  }
  expected_freq_c <- expected_freq_c*nsim_c
  expected_freq_e <- expected_freq_e*nsim_e
  
  if (prod(expected_freq_c) == 0){
    expected_freq_c <- expected_freq_c + rep(5,6)
    print("There is some bin containing 0 in correct responses")
  } 
  
  if(prod(expected_freq_e) == 0){
    expected_freq_e <- expected_freq_e + rep(5,6)
    print("There is some bin containing 0 in error responses")
  }
  
  chisquare_c <- sum((expected_freq_c-p*nsim_c)^2/expected_freq_c)
  chisquare_e <- sum((expected_freq_e-p*nsim_e)^2/expected_freq_e)
  #log(chisquare_c + chisquare_e)
  
  return(log(chisquare_c + chisquare_e))
}

