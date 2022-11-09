

DDM_traj <- simul_DDM(v = 0.5, a = 1, init_loc = 0.5)$traj
len <- length(DDM_traj)

par(mfrow = c(1, 1), mar = c(3,3,3,3))

plot(201:(len+200), DDM_traj, xlim = c(0,1500), type = "l", ylim = c(0,1.1), xlab = "", ylab = "") #,  xaxt='n', yaxt='n')
lines(1:2500, rep(1,2500), col = 'red', lwd = 2)
lines(1:2500, rep(0,2500), col = 'red', lwd = 2)
lines(seq(0,200,1), rep(0.5,201), lwd = 3)


#lines(c(0,500), c(0.5,0.75),type = "l", ylim = c(0,1.2), xlab = "Time (ms)", ylab = "Evidence", col = 'red')
arrows(201, 0.5, x1 = 500, y1 = 0.75, length = 0.1, col = c("blue", "blue", "blue"), lwd = 3)

mtext("Time (ms)", side = 1, line = 2)       
mtext("Evidence", side = 2, line = 2)       
