library(plotly) # contours
library(raster) # vetor field
library(rasterVis)

## analytically plot y1 nullcline
y1 <- seq(0.001, 0.999, 0.001)
y2 <- (W_pos*2*y1 + (B1-theta) - beta^(-1)*0.5*N*(log(y1/(1-y1))))/W_neg
plot(y1,y2, ylim = c(0,1))

## visualize y1 and y2 gradient seperately
B1 <- 3300
B2 <- 2700
resolution <- 0.005
dim_mat <- 1/resolution - 1
y1 <- seq(resolution,1-resolution, resolution)
y2 <- seq(resolution,1-resolution, resolution)

# gradient in y1 direction
dy1s <- matrix(NA, dim_mat,dim_mat)
for (i in 1:dim_mat){
  for (j in 1:dim_mat){
    #### here the index is reversed!!!!!!
    dy1s[j,i] <- sign(dFdy1(y1[i],y2[j], W_pos, W_neg, B1 = B1, theta, beta, N))#, W_pos, W_neg, B1, theta, beta, N)
  }
}
fig <- plot_ly(
  x = y1, 
  y = y2, 
  z = dy1s, 
  type = "contour",
  contours = list(start = min(dy1s), end = max(dy1s), size = 1000)
)
fig

# gradient in y2 direction
dy2s <- matrix(NA, dim_mat,dim_mat)
for (i in 1:dim_mat){
  for (j in 1:dim_mat){
    dy2s[i,j] <- sign(dFdy2(y1[i],y2[j], W_pos, W_neg, B2 = B2, theta, beta, N))#, W_pos, W_neg, B1, theta, beta, N)
  }
}
fig <- plot_ly(
  x = y1, 
  y = y2, 
  z = dy2s, 
  type = "contour",
  contours = list(start = min(dy1s), end = max(dy1s), size = 1000)
)
fig



##
df <- expand.grid(x=seq(0.005,0.995,0.005), y=seq(0.005,0.995,0.005))
colnames(df)
df$z <- with(df, free_energy(x, y, B1=B1, B2=B2))

r <- rasterFromXYZ(df)
projection(r) <- CRS("+init=EPSG:4326")
vectorplot(r, par.settings=RdBuTheme)

## use ggplot2
# Need to load grid for arrow() function
library(grid)
#df <- expand.grid(x=seq(0.01,0.99,0.01), y=seq(0.01,0.99,0.01))
df <- expand.grid(x=seq(0.05,0.95,0.05), y=seq(0.05,0.95,0.05))
#df$z <- with(df, free_energy(x, y, B1=B1, B2=B2))
df$vx <- dFdy1(df$x,df$y, W_pos, W_neg, B1, theta, beta, N)
df$vy <- dFdy2(df$x,df$y, W_pos, W_neg, B2, theta, beta, N)
df$speed <- sqrt(df$vx^2 + df$vy^2)
df$xend <- df$x-df$vx/df$speed/40
df$yend <- df$y-df$vy/df$speed/40
colnames(df)
# Make the plot with the subset, and use an arrowhead 0.1 cm long
ggplot(df) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.1, "cm")), size = 0.2) +
  xlim(0,1) + 
  ylim(0,1) +
  xlab("y1") +
  ylab("y2")

## visualize free energy surface
y1 <- seq(0.005,0.995,0.005)
y2 <- seq(0.005,0.995,0.005)
B1 <- 3500
B2 <- 3500
Fys <- matrix(NA, 199,199)
W_pos <- 50000
for (i in 1: 199){
  for (j in 1: 199){
    Fys[i,j] <- free_energy(y1[i],y2[j], B1=B1, B2=B2, W_pos = W_pos) #, W_neg, B1, theta, beta, N)
  }
}

fig <- plot_ly(
  x = y1, 
  y = y2, 
  z = Fys, 
  type = "contour",

  contours = list(start = min(Fys), end = max(Fys), size = 100)
)
fig <- fig %>% layout(title = paste("Self-excitation is", W_pos))
fig

#### sequential trajectory
W_pos <- 50000
result <- simul_continuous(W_pos = W_pos)
plot_trajectory(result, W_pos)

### simulate a series of decisions together
x <- simul_continuous_multiple(simulations = 2000, W_pos = 51000, C = 0.1)
x

### plot RT distribution
hist(x[x$Response==2,]$RT, prob = TRUE,xlab = "Reactioin Time", main = "Reaction Time Distribution")
hist(x[x$Response==1,]$RT, prob = TRUE, add = TRUE, col = rgb(1, 0, 0, 0.5))

lines(density(x[x$Response==2,]$RT), lwd = 2, col = 'black')
lines(density(x[x$Response==1,]$RT), lwd = 2, col = 'red')




### Conditional response function (CRF)
breaks <- seq(0.2,1.6,0.1)
h_c <- hist(x[x$Response==2,]$RT, prob = TRUE,
            xlab = "Reactioin Time", main = "Reaction Time Distribution",
            breaks=breaks)
h_e <- hist(x[x$Response==1,]$RT, prob = TRUE, add = TRUE, col = rgb(1, 0, 0, 0.5),
            breaks=breaks)

plot(h_c$mids,(h_c$counts/(h_c$counts + h_e$counts)), ylim=c(0,1), 
     xlab = "RT", ylab = "Proportion of Correct Responses", type = "l")




length(x[x$Response==1,]$RT)





