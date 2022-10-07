library("tidyverse")
library("ggplot2")
library("tidyr")
#v = 1
#init_loc = 0.4
#simul_DDM_non_stop(v = 0, init_loc = init_loc, n_iter = 2000)

n_sim <- 100

vs <- seq(-2,3,0.3)

states_all <-c()
for (i in 1:n_sim){
  print(i)
  init_loc <- 0.4
  states <- c()
  for(v in vs){
    result_iter <- simul_DDM_non_stop(v = v, init_loc = init_loc, n_iter = 2000)
    init_loc <- result_iter$loc
    state <- result_iter$state
    states <- c(states, state)
    #print(state)
  }
  states_all <- rbind(states_all, states)
}


df <- data.frame("v" = vs, "state0" = 0,"state1" = 0, "state2" = 0)

for (i in 1:length(vs)){
  table_iter <- as.data.frame(table(states_all[,i]))
  #print(table_iter)
  for (j in 1:length(table_iter$Var1)){
    state_iter <- paste("state", table_iter$Var1[j], sep = "") 
    freq_iter <- table_iter$Freq[j]
    df[state_iter][i,] <- freq_iter/n_sim
  }
}

df1 <- gather(df, key = "state", value = "p", state0, state1, state2)

ggplot(data=df1, aes(x=v, y=p, group=state,colour=state)) +
  xlab("Drift Rate (v)") + 
  ylab("Proportion") +
  geom_line() + 
  geom_point() +
  scale_colour_discrete(name  ="System State",
                        #breaks=c("Female", "Male"),
                        labels=c("state 0", "State 1", "State 2"))


##### backward
vs <- seq(2,-3,-0.3)

states_all <-c()
for (i in 1:n_sim){
  print(i)
  init_loc <- 0.4
  states <- c()
  for(v in vs){
    result_iter <- simul_DDM_non_stop(v = v, init_loc = init_loc, n_iter = 2000)
    init_loc <- result_iter$loc
    state <- result_iter$state
    states <- c(states, state)
    #print(state)
  }
  states_all <- rbind(states_all, states)
}


df <- data.frame("v" = vs, "state0" = 0,"state1" = 0, "state2" = 0)

for (i in 1:length(vs)){
  table_iter <- as.data.frame(table(states_all[,i]))
  #print(table_iter)
  for (j in 1:length(table_iter$Var1)){
    state_iter <- paste("state", table_iter$Var1[j], sep = "") 
    freq_iter <- table_iter$Freq[j]
    df[state_iter][i,] <- freq_iter/n_sim
  }
}

df1 <- gather(df, key = "state", value = "p", state0, state1, state2)

ggplot(data=df1, aes(x=v, y=p, group=state,colour=state)) +
  xlab("Drift Rate (v)") + 
  ylab("Proportion") +
  geom_line() + 
  geom_point() +
  scale_colour_discrete(name  ="System State",
                        #breaks=c("Female", "Male"),
                        labels=c("state 0", "State 1", "State 2"))



vs <- seq(2,-3,-0.3)

states_all <-c()
for (i in 1:n_sim){
  print(i)
  init_loc <- 0.4
  states <- c()
  for(v in vs){
    result_iter <- simul_DDM_non_stop(v = 0.3, init_loc = init_loc, n_iter = 1000)
    init_loc <- result_iter$loc
    state <- result_iter$state
    states <- c(states, state)
    #print(state)
  }
  states_all <- rbind(states_all, states)
}



