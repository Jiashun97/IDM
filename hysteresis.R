library("tidyverse")
library("ggplot2")
library("tidyr")


n_sim <- 100

Cs <- seq(-1,1,0.1)

states_all <-c()
for (i in 1:n_sim){
  print(i)
  init_loc <- c(0.3,0.3)
  states <- c()
  for(C in Cs){
    result_iter <- simul_IDM_non_stop(C = C, init_loc = init_loc, n_iter = 2000, Bns = 2500, Bs = 2000)
    init_loc <- result_iter$loc
    state <- result_iter$state
    states <- c(states, state)
  }
  states_all <- rbind(states_all, states)
}


df <- data.frame("C" = Cs, "state0" = 0, "state1" = 0, "state2" = 0)

for (i in 1:length(Cs)){
  table_iter <- as.data.frame(table(states_all[,i]))
  #print(table_iter)
  for (j in 1:length(table_iter$Var1)){
    state_iter <- paste("state", table_iter$Var1[j], sep = "") 
    freq_iter <- table_iter$Freq[j]
    df[state_iter][i,] <- freq_iter/n_sim
  }
}

df1 <- gather(df, key = "state", value = "p", state0, state1, state2)

ggplot(data=df1, aes(x=C, y=p, group=state,colour=state)) +
  xlab("Stimulus Strength (C)") + 
  ylab("Proportion") +
  geom_line() + 
  geom_point() +
  scale_colour_discrete(name  ="System State",
                        #breaks=c("Female", "Male"),
                        labels=c("State 0", "State 1", "State 2"))


##### backward
Cs <- seq(1,-1,-0.1)

states_all <-c()
for (i in 1:n_sim){
  print(i)
  init_loc <- c(0.3,0.3)
  states <- c()
  for(C in Cs){
    result_iter <- simul_IDM_non_stop(C = C, init_loc = init_loc, n_iter = 2000, Bns = 2500, Bs = 2000)
    init_loc <- result_iter$loc
    state <- result_iter$state
    states <- c(states, state)
  }
  states_all <- rbind(states_all, states)
}


df_backward <- data.frame("C" = Cs, "state0" = 0, "state1" = 0, "state2" = 0)

for (i in 1:length(Cs)){
  table_iter <- as.data.frame(table(states_all[,i]))
  #print(table_iter)
  for (j in 1:length(table_iter$Var1)){
    state_iter <- paste("state", table_iter$Var1[j], sep = "") 
    freq_iter <- table_iter$Freq[j]
    df_backward[state_iter][i,] <- freq_iter/n_sim
  }
}



df1_backward <- gather(df_backward, key = "state", value = "p", state0, state1, state2)

ggplot(data=df1_backward, aes(x=C, y=p, group=state,colour=state)) +
  xlab("Stimulus Strength (C)") + 
  ylab("Proportion") +
  geom_line() + 
  geom_point() +
  scale_colour_discrete(name  ="System State",
                        #breaks=c("Female", "Male"),
                        labels=c("State 0", "State 1", "State 2"))





