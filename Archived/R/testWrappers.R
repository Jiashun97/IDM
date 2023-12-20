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
#   $Id: testWrappers.R 69 2015-03-30 16:35:24Z u0066818@kuleuven.be $

testWrappers <- function(nGPUs=0,gpuIds=c(),loadPerGPU=c()){
  settings=c();
  settings$nGPUs=nGPUs;
  settings$gpuIds=gpuIds;
  settings$loadPerGPU=loadPerGPU;
  settings$nWalkers=1000000;
  settings$seed=1;
  settings=checkSettings(settings,1);
  
  #LDMPars
  LDMPars=c();
  LDMPars$nStimuli=4;
  LDMPars$a=0.08;
  LDMPars$v=matrix(c(0.4,0.25,0.1,0),1,LDMPars$nStimuli);
  LDMPars$Ter=0.3;
  LDMt=system.time(LDMout<-LDMDist(LDMPars,settings,1))[[3]];
  
  #LCAPars
  LCAPars=c();
  LCAPars$nDim=2;
  LCAPars$nStimuli=4;
  LCAPars$a=matrix(c(0.08,0.08),LCAPars$nDim,1);
  LCAPars$v=matrix(c(0,0.05,0.1,0.15,0.3,0.25,0.2,0.15),LCAPars$nDim,LCAPars$nStimuli,byrow=TRUE);
  LCAPars$Ter=0.3;
  LCAPars$Gamma=matrix(c(0.1,0.1,0.1,0.1),LCAPars$nDim,LCAPars$nDim);
  LCAt=system.time(LCAout<-LCADist(LCAPars,settings,1))[[3]];
  
  #IDMPars
  IDMPars=c();
  IDMPars$nDim=2;
  IDMPars$nStimuli=4;
  IDMPars$N=matrix(1000,IDMPars$nDim,1);
  IDMPars$Theta=matrix(51450,IDMPars$nDim,1);
  IDMPars$W=matrix(c(52500,8400,8400,52500),IDMPars$nDim,IDMPars$nDim);
  IDMPars$B=matrix(c(2000,2250,2500,2750,3500,3250,3000,2750),IDMPars$nDim,IDMPars$nStimuli,byrow=TRUE);
  IDMPars$h=matrix(0.4,IDMPars$nDim,IDMPars$nDim);
  IDMPars$dynamics=0;
  IDMPars$D=matrix(0.05,IDMPars$nDim,1);
  IDMPars$spontaneousTime=1;
  IDMPars$Ter=0.3;
  IDMt=system.time(IDMout<-IDMDist(IDMPars,settings,1))[[3]];
  
  LDMPdfs0=LDMout$distributions/settings$nWalkers;
  LCAPdfs0=LCAout$distributions/settings$nWalkers;
  IDMPdfs0=IDMout$distributions/settings$nWalkers;
  #LDMPdfs=LDMPdfs0;LCAPdfs=LCAPdfs0;IDMPdfs=IDMPdfs0;
  #save('LDMPdfs','LCAPdfs','IDMPdfs',file='testPdfs.RData');
  load('testPdfs.RData');

  par(mfrow=c(1,3));
  LDMLlh=sum(LDMPdfs0*log(pmax(LDMPdfs,0.000001)));
  if ((-24.90<LDMLlh)&&(LDMLlh<(-24.88))){
    print(paste('LDMDist passed test in',format(LDMt,nsmall=4),'sec'));
    plot(matrix(1:3000,3000,2*LDMPars$nStimuli),LDMPdfs,type='l');
    selection=(0:299)*10+1;
    points(matrix(selection,length(selection),2*LDMPars$nStimuli),LDMPdfs0[selection,],pch=20);
    title('LDMDist');
  }else{
    print('LDMDist did not pass test');
    allOk=FALSE;
  }
  
  LCALlh=sum(LCAPdfs0*log(pmax(LCAPdfs,0.000001)));
  if ((-26.25<LCALlh)&&(LCALlh<(-26.23))){
    print(paste('LCADist passed test in',format(LCAt,nsmall=4),'sec'));
    plot(matrix(1:3000,3000,2*LCAPars$nStimuli),LCAPdfs,type='l');
    selection=(0:299)*10+1;
    points(matrix(selection,length(selection),LCAPars$nDim*LCAPars$nStimuli),LCAPdfs0[selection,],pch=20);
    title('LCADist');
  }else{
    print('LCADist did not pass test');
    allOk=FALSE;
  }
  
  IDMLlh=sum(IDMPdfs0*log(pmax(IDMPdfs,0.000001)));
  if ((-22.05<IDMLlh)&&(IDMLlh<(-22.03))){
    print(paste('IDMDist passed test in',format(IDMt,nsmall=4),'sec'));
    plot(matrix(1:3000,3000,IDMPars$nDim*IDMPars$nStimuli),IDMPdfs,type='l');
    selection=(0:299)*10+1;
    points(matrix(selection,length(selection),IDMPars$nDim*IDMPars$nStimuli),IDMPdfs0[selection,],pch=20);
    title('IDMDist');
  }else{
    print('IDMDist did not pass test');
    allOk=FALSE;
  }
}
