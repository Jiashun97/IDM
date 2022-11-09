library(ggplot2)

## simulation
Cs <- seq(0,0.3,0.015) 
plot_result <- matrix(data = NA, nrow = 3, ncol = length(Cs))

h <- 0.5
durations <- seq(500,1500,500)

for (j in 1:3){
  print(j)
  duration <- durations[j]
  
  n_sim <- 5000 #######
  props_uncertain <- c()
  props_change <- c()
  
  for(C in Cs){
    cat('C = ', C)
    uncertain <- c()
    change <- c()
    iter <- c()
    for (i in 1: n_sim){
      result <- simul_IDM_com(C = C, h = h, duration = duration)
      uncertain <- c(uncertain, result$uncertain)
      change <- c(change, result$change)
      iter <- c(iter, result$iter)
    }
    
    props_uncertain <- c(props_uncertain, sum(uncertain)/n_sim)
    props_change <- c(props_change, sum(change)/n_sim)
    
  }
  plot_result[j,] <- props_change
}


df <- data.frame(x = Cs,
                 val = c(plot_result[1,], plot_result[2,], plot_result[3,]), 
                 variable = rep(c("0.5 s", "1 s", "1.5 s"), each=21))
write.csv(df,'IDM_com_result.csv')

plot <- ggplot(data = df, aes(x=x, y=val, colour=variable)) + 
  geom_line() + 
  #scale_fill_discrete(breaks=c("500 ms", "1000 ms", "2000 ms"))
  xlab("Stimulus Strength") + 
  ylab("P(Change of Mind)") +
  labs(colour = "Post-decision Period")
plot

#plot(vs, props_change, type = "l", xlab = "Drift Rate", ylab = "P(Change of Mind)")

#table(uncertain)
#table(change)








