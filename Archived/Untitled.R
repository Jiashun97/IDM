


#create data frame
nsim_c <- length(data_RT[data_RT$Response == 2, "RT"])
nsim_e <- length(data_RT[data_RT$Response == 1, "RT"])
data_RT_q_c <- as.numeric(quantile(data_RT[data_RT$Response == 2, "RT"], probs = q))
data_RT_q_e <- as.numeric(quantile(data_RT[data_RT$Response == 1, "RT"], probs = q))

df <- data.frame(nsim_c=nsim_c,
                 nsim_e=nsim_e,
                 data_RT_q_c=data_RT_q_c,
                 data_RT_q_e=data_RT_q_e)

data <- c(nsim_c, nsim_e, data_RT_q_c, data_RT_q_e)

get_chisquare_stats <- function(data, rt_dist){
  q = c(0.1, 0.3, 0.5, 0.7, 0.9)
  p = c(.1,.2,.2,.2,.2,.1)
  
  nsim_c <- data[1]
  nsim_e <- data[2]
  # empirical RT quantile
  data_RT_q_c <- data[3:7]
  data_RT_q_e <- data[8:12]
  
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
  chisquare_c <- sum((expected_freq_c-p*nsim)^2/expected_freq_c)
  chisquare_e <- sum((expected_freq_e-p*nsim)^2/expected_freq_e)
  
  return(log(chisquare_c + chisquare_e))
}





#create data frame
df <- data.frame(x=c(1, 3, 3, 5, 6, 7, 9, 12),
                 y=c(4, 5, 8, 6, 9, 10, 13, 17))

#define function to minimize residual sum of squares
min_residuals <- function(data, par) {
  with(data, sum((par[1] + par[2] * x - y)^2))
}

#find coefficients of linear regression model
optim(par=c(0, 1), fn=min_residuals, data=df)

library(RWiener)
### another way to simulate wiener process
rwiener(n = 100, alpha = 2,tau = 0.3,beta = 0.5,delta = 0.5)

### brownian path simulation
randn('state',100)           % set the state of randn
T = 1; N = 500; dt = T/N;
dW = zeros(1,N);             % preallocate arrays ...
W = zeros(1,N);              % for efficiency

dW(1) = sqrt(dt)*randn;      % first approximation outside the loop ...
W(1) = dW(1);                % since W(0) = 0 is not allowed
for j = 2:N
dW(j) = sqrt(dt)*randn;   % general increment
W(j) = W(j-1) + dW(j); 
end

plot([0:dt:T],[0,W],'r-')    % plot W against t
xlabel('t','FontSize',16) 
ylabel('W(t)','FontSize',16,'Rotation',0)


set.seed(1)
T <- 1
N <- 500
dt <- T/N
dW <- numeric(N)
W <- numeric(N)

# first iteration
dW[1] <- sqrt(dt)*rnorm(1)      # first approximation outside the loop ...
W[1] <- dW[1]

for (i in 2:N){
  dW[i] <- sqrt(dt)*rnorm(1)
  W[i] <- W[i-1] + dW[i]
}

plot(1:N,W)


simul_IDM <- function(nsim, prob_scaled_vec, C, N = 2000, W_pos = 52500, W_neg = 8400, Bns = 2500, Bs = 1000, # mode = "coarse-grained"
                      beta = 1/24, theta = 51450, Ter = 0.3, h = 0.4, D = 0.05, sigma = 0.01, dt = 0.001, mode = "Euler-Maruyama"){ # D = 0.05, 2/(1000*exp(1))
  
  B1 <- Bs*(1+C) + Bns
  B2 <- Bs*(1-C) + Bns
  #start_time <- Sys.time() 
  result <- c()
  iters <- c()
  for(i in 1:nsim){
    sample_y <- sample_from_Py(N/2,prob_scaled_vec)
    y1 <- sample_y[1]
    y2 <- sample_y[2]
    iter <- 1
    #location <- c(y1, y2)
    if(mode == "Euler-Maruyama"){
      w_pre1 <- 0
      w_pre2 <- 0
      w_post1 <- rnorm(1,0,sqrt(dt))
      w_post2 <- rnorm(1,0,sqrt(dt))
      while(((y1<(1-h))&(y2<(1-h)))|((y1>h)&(y2>h))){
        #dE1 <- dEdy1(y1,y2, W_pos, W_neg, B1, theta)
        #dE2 <- dEdy2(y1,y2, W_pos, W_neg, B2, theta)
        #D1 <- (1-sign(dE1)*(1-2*y1))/(N^2*dt)  
        #D2 <- (1-sign(dE2)*(1-2*y2))/(N^2*dt)
        #cat("D1", D1)
        
        dy1 <- -beta*D*dFdy1(y1,y2, W_pos, W_neg, B1, theta, beta, N)*dt + sqrt(2*D)*(w_post1-w_pre1)#(wiener_post1 - wiener_pre1)
        dy2 <- -beta*D*dFdy2(y1,y2, W_pos, W_neg, B2, theta, beta, N)*dt + sqrt(2*D)*(w_post2-w_pre2)#(wiener_post2 - wiener_pre2)
        
        w_pre1 <- w_post1
        w_pre2 <- w_post2
        w_post1 <- rcmvnorm(n=1, mean=c(0,0), sigma=matrix(c(iter*dt,iter*dt,iter*dt,(iter+1)*dt),2,2), 
                            dependent.ind = 2, given.ind = 1, X = w_pre1, method="chol")
        w_post2 <- rcmvnorm(n=1, mean=c(0,0), sigma=matrix(c(iter*dt,iter*dt,iter*dt,(iter+1)*dt),2,2), 
                            dependent.ind = 2, given.ind = 1, X = w_pre2, method="chol")
        
        y1 <- clamp(y1 + dy1) # clamp to [0,1]
        y2 <- clamp(y2 + dy2)
        #location <- cbind(location, c(y1,y2))
        iter <- iter + 1
      }
    }else if(mode == "coarse-grained"){
      while(((y1<(1-h))&(y2<(1-h)))|((y1>h)&(y2>h))){
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
    
    
    if ((y1>(1-h))&(y2<h)){
      if(B1>=B2){
        R <- 2 # correct response
      }else if(B1<B2){
        R <- 1 # incorrect response
      }
      result <- rbind(result, c(RT, R))
    }
    
    if ((y1<h)&(y2>(1-h))){
      if(B1>=B2){
        R <- 1 # incorrect response
      }else if(B1<B2){
        R <- 2 # correct response
      }
      result <- rbind(result, c(RT, R))
    }
    iters <- c(iters, iter)
    #print(i)
    #print(iter)
    if ((i%%100) == 0){
      #cat("simulation ",i," has finished. ")
      #cat(iter, "iterations.")
    }
  }
  result_df <- data.frame(result)
  colnames(result_df) <- c('RT','Response')
  #print(Sys.time() - start_time)
  cat("Averaged teration is ", mean(iters), ". ")
  
  return(result_df)
}




stim_t <- 5000 
  rest_t <- 2000 
  init_loc <- c(0.30.3) 
  C <- 0.2 
  N <- 2000 
  W_pos <- 52500 
  W_neg <- 8400 
  Bns <- 2500 
  Bs <- 1000 
  beta <- 1/24 
  theta <- 51450 
  Ter <- 0.3 
  h <- 0.4 
  D <- 0.05 
  sigma <- 0.01 
  dt <- 0.001 
  mode <- "Euler-Maruyama"

