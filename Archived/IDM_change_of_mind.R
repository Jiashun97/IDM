library(ggplot2)
library(dplyr)
## simulation
Cs <- seq(0,0.3,0.015) 
plot_result <- matrix(data = NA, nrow = 3, ncol = length(Cs))

h <- 0.5
durations <- seq(500,1500,500)

for (j in 1:3){
  print(j)
  duration <- durations[j]
  
  n_sim <- 1000 #######
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
                 variable = rep(c("0.5 s", "1 s", "1.5 s"), each=length(Cs)))

plot <- ggplot(data = df, aes(x=x, y=val, colour=variable)) + 
  geom_line() + 
  #scale_fill_discrete(breaks=c("500 ms", "1000 ms", "2000 ms"))
  ylim(0, 0.2) +
  xlab("Stimulus Strength") + 
  ylab("P(Change of Mind)") +
  labs(colour = "Post-decision Period")
plot

#plot(vs, props_change, type = "l", xlab = "Drift Rate", ylab = "P(Change of Mind)")

#table(uncertain)
#table(change)

write.csv(df,'IDM_com_result_h_0_5.csv')

## plot previous graph

df <- read.csv('IDM_com_result_h_0_6.csv')
plot <- ggplot(data = df, aes(x=x, y=val, colour=variable)) + 
  geom_line() + 
  #scale_fill_discrete(breaks=c("500 ms", "1000 ms", "2000 ms"))
  ylim(0, 0.2) +
  xlab("Stimulus Strength") + 
  ylab("P(Change of Mind)") +
  labs(colour = "Post-decision Period")
plot


### why IDM change of mind does not depend on the post-decisional time?


C <- 0.3
h <- 0.5
duration <- 500

n_sim <- 10000
result <- simul_IDM_com(C = C, h = h, duration = duration)
results_df <- data.frame(change = result$change, iter = result$iter, 
                        loc_init_y1 = result$loc_init[1],loc_init_y2 = result$loc_init[2],
                        loc_com_y1 = result$loc_com[1],loc_com_y2 = result$loc_com[2], 
                        iter_post = result$iter_post)
for (i in 1: n_sim){
  print(i)
  result <- simul_IDM_com(C = C, h = h, duration = duration)
  result_df <- data.frame(change = result$change, iter = result$iter, 
                          loc_init_y1 = result$loc_init[1],loc_init_y2 = result$loc_init[2],
                          loc_com_y1 = result$loc_com[1],loc_com_y2 = result$loc_com[2], 
                          iter_post = result$iter_post)
  results_df <- rbind(results_df, result_df)

}

#plot(results_df$loc_com_y1, results_df$loc_com_y2)



ggplot(data = filter(results_df, change == 1), aes(x=loc_com_y1, y=loc_com_y2)) + 
  geom_point() + 
  #scale_fill_discrete(breaks=c("500 ms", "1000 ms", "2000 ms"))
  xlim(0, 1) +
  ylim(0, 1) +
  xlab("Stimulus Strength") + 
  ylab("P(Change of Mind)")
  #labs(colour = "Post-decision Period")


ggplot(data = filter(results_df, change == 0), aes(x=loc_init_y1, y=loc_init_y2)) + 
  geom_point() + 
  #scale_fill_discrete(breaks=c("500 ms", "1000 ms", "2000 ms"))
  xlim(0, 1) +
  ylim(0, 1) +
  xlab("Stimulus Strength") + 
  ylab("P(Change of Mind)")
#labs(colour = "Post-decision Period")


ggplot(data = results_df, aes(x=loc_com_y1, y=loc_com_y2, color = change_str)) + 
  geom_point() + 
  #scale_fill_discrete(breaks=c("500 ms", "1000 ms", "2000 ms"))
  xlim(0, 1) +
  ylim(0, 1) +
  xlab("Stimulus Strength") + 
  ylab("P(Change of Mind)")
#labs(colour = "Post-decision Period")


results_df <- results_df %>%
  mutate(change_str = if_else(change == 1, "Change of Mind", "No Change")
  )


plot(x=results_df$loc_init_y1, y=results_df$loc_init_y2, xlim = c(0,1), ylim = c(0,1), 
     col = alpha("black", 0.4), cex = 0.1, xlab = "y1", ylab = "y2")
points(x=results_df$loc_com_y1, y=results_df$loc_com_y2, col = alpha("red", 0.4), cex = 0.5, pch=19) # solid circle
legend(0.65, 0.8, legend=c("Initial Decision", "Change of Mind"),
       col=c("black", "red"), lty=1:2, cex=0.8)


plot(x=results_df$loc_init_y1, y=results_df$loc_init_y2, xlim = c(0,1), ylim = c(0,1), 
     col = alpha("black", 0.4), cex = 0.1, xlab = "y1", ylab = "y2")
points(x=results_df$loc_init_y1, y=results_df$loc_init_y2, col = alpha("red", 0.4), cex = 0.1, pch=19) # solid circle
legend(0.65, 0.8, legend=c("Initial Decision", "Change of Mind"),
       col=c("black", "red"), lty=1:2, cex=0.8)

plot(x=filter(results_df, change == 1)$iter, y=filter(results_df, change == 1)$iter_post, #xlim = c(0,1), ylim = c(0,1), 
     xlab = "RT of Initial Decision", ylab = "RT of Change of Mind")
cor(x=filter(results_df, change == 1)$iter, y=filter(results_df, change == 1)$iter_post)


results_changeofmind <- filter(results_df, change == 1)

results_changeofmind$distance <- pmin(abs(results_changeofmind$loc_com_y1 - 0.5), abs(results_changeofmind$loc_com_y2 - 0.5))


plot(x=results_changeofmind$distance, y=results_changeofmind$iter_post, #xlim = c(0,1), ylim = c(0,1), 
     xlab = "Distance to Detectioin Box", ylab = "RT of Change of Mind")
cor(x=results_changeofmind$distance, y=results_changeofmind$iter_post)




