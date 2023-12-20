#   This file is part of the RTDist project
#   Copyright (c) 2014 Stijn Verdonck
#   Copyright (c) 2014 Kristof Meers
#
#   Verdonck, S., Meers, K., & Tuerlinckx, F. (in press). Efficient simulation
#       of diffusion-based choice RT models on CPU and GPU. Behavior Research
#       Methods. doi:10.3758/s13428-015-0569-0
# 
#   RTDist comes without any warranty of any kind. You are not allowed to
#   redistribute a copy of RTDist to others. If you want others to use RTDist,
#   refer them to http://ppw.kuleuven.be/okp/software/RTDist/. See the root
#   folder of this project for full license information in the LICENSE.txt file.
#
#   $Id: RTDistTestIndices.R 69 2015-03-30 16:35:24Z u0066818@kuleuven.be $

#test indices
IDMPars=c();
IDMPars$nDim=2;
IDMPars$N=matrix(1000,2,1);
IDMPars$W=matrix(c(52500,8400,8400,52500),2,2);
IDMPars$Theta=matrix(51450,2,1);
IDMPars$nStimuli=1;
IDMPars$B=matrix(2500,2,1);
IDMPars$h=matrix(0.4,2,2);
IDMPars$Ter=0.3;
IDMPars$spontaneousTime=1;
IDMPars$dynamics=0;
IDMPars$D=matrix(0.05,2,1);
IDMPars=checkIDMPars(IDMPars);
settings=c();
settings$seed=0;
settings=checkSettings(settings);

settings$nWalkers=100000;
settings$nGPUs=0;

#check indices B
IDMParsTest=IDMPars;
IDMParsTest$nStimuli=2;
#each column is a stimulus 
IDMParsTest$B=t(matrix(c(3000,3500,2000,1500),2,2));

out=IDMDist(IDMParsTest,settings);
distributions=out$distributions;
counts=c();
for(i in 1:4){
  counts[i]=sum(distributions[,i]);
}
(counts[1]-counts[2])/(counts[1]+counts[2])
(counts[3]-counts[4])/(counts[3]+counts[4])

#check indices h
IDMParsTest=IDMPars;
#each column is a box
IDMParsTest$h=t(matrix(c(0.4,0.4,1.0,0.4),2,2));
IDMParsTest$B=matrix(c(5000,5000),2,1);

out=IDMDist(IDMParsTest,settings);
distributions=out$distributions;
counts=c();
for(i in 1:2){
  counts[i]=sum(distributions[,i]);
}
(counts[1]-counts[2])/(counts[1]+counts[2])
