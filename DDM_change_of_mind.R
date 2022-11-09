library(ggplot2)
## simulation

plot_result <- matrix(data = NA, nrow = 3, ncol = 21)

a <- 0.4
post_iters <- seq(500,1500,500)

for (j in 1:3){
  print(j)
  post_iter <- post_iters[j]
  
  n_sim <- 1000
  init_loc <- a*0.5
  vs <- seq(0,2,0.1) 
  props_uncertain <- c()
  props_change <- c()
  
  for(v in vs){
    
    #v <- v
    
    uncertain <- c()
    change <- c()
    iter <- c()
    for (i in 1: n_sim){
      result <- simul_DDM_com(v = v, post_iter = post_iter, a = a, init_loc = init_loc)
      uncertain <- c(uncertain, result$uncertain)
      change <- c(change, result$change)
      iter <- c(iter, result$iter)
    }
    
    props_uncertain <- c(props_uncertain, sum(uncertain)/n_sim)
    props_change <- c(props_change, sum(change)/n_sim)
    
  }
  plot_result[j,] <- props_change
}


df <- data.frame(x = vs,
                 val = c(plot_result[1,], plot_result[2,], plot_result[3,]), 
                 variable = rep(c("0.5 s", "1 s", "1.5 s"), each=21))
plot <- ggplot(data = df, aes(x=x, y=val, colour=variable)) + 
  geom_line() + 
  #scale_fill_discrete(breaks=c("500 ms", "1000 ms", "2000 ms"))
  xlab("Drfit Rate") + 
  ylab("P(Change of Mind)") +
  labs(colour = "Post-decision Period")
plot

#plot(vs, props_change, type = "l", xlab = "Drift Rate", ylab = "P(Change of Mind)")

#table(uncertain)
#table(change)
