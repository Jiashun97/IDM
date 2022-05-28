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
#   $Id: IDMWrapper.R 69 2015-03-30 16:35:24Z u0066818@kuleuven.be $

checkIDMPars <- function(pars,verbose=2){
  pars=as.list(pars);
  allOk=TRUE;  
  
  pars=checkParField(pars,'nDim','integer',c(),c(1,1),1,c(),verbose);allOk=pars$ok&&allOk;
  pars=checkParField(pars,'nStimuli','integer',c(),c(1,1),1,c(),verbose);allOk=pars$ok&&allOk;
  if (allOk){
    nDim=pars$nDim;
    nStimuli=pars$nStimuli;
    pars=checkParField(pars,'startPos','double',matrix(0.3,nDim,1),c(nDim,1),0.001,0.999,verbose);allOk=pars$ok&&allOk;
    pars=checkParField(pars,'beta','double',as.double(1/24),c(1,1),0,c(),verbose);allOk=pars$ok&&allOk;
    pars=checkParField(pars,'N','double',c(),c(nDim,1),0,c(),verbose);allOk=pars$ok&&allOk;
    pars=checkParField(pars,'W','double',c(),c(nDim,nDim),0,c(),verbose);allOk=pars$ok&&allOk;
    if (pars$ok){
      if (!min(pars$W==t(pars$W))){
        if (verbose>=1){
          print('ERROR: W not symmetric');
        }
        allOk=FALSE;
      }
    }
    pars=checkParField(pars,'Theta','double',c(),c(nDim,1),0,c(),verbose);allOk=pars$ok&&allOk;
    pars=checkParField(pars,'boxShape','integer',as.integer(0),c(1,1),0,2,verbose);allOk=pars$ok&&allOk;
    pars=checkParField(pars,'h','double',c(),c(nDim,nDim),0,c(),verbose);allOk=pars$ok&&allOk;
    pars=checkParField(pars,'B','double',c(),c(nDim,nStimuli),c(),c(),verbose);allOk=pars$ok&&allOk;
    pars=checkParField(pars,'Ter','double',as.double(0),c(1,1),0,c(),verbose);allOk=pars$ok&&allOk;
    if (allOk){
      pars=checkParField(pars,'sTer','double',as.double(0),c(1,1),0,2*pars$Ter,verbose);allOk=pars$ok&&allOk;
    }
    pars=checkParField(pars,'spontaneousTime','double',c(),c(1,1),0,c(),verbose);allOk=pars$ok&&allOk;
    pars=checkParField(pars,'dynamics','integer',c(),c(1,1),0,1,verbose);allOk=pars$ok&&allOk;
    if (pars$ok){
      if (pars$dynamics==0){
        pars=checkParField(pars,'D','double',c(),c(nDim,1),0,c(),verbose);allOk=pars$ok&&allOk;
      }else if (pars$dynamics==1){
        pars=checkParField(pars,'sigma','double',c(),c(nDim,1),0,c(),verbose);allOk=pars$ok&&allOk;
        pars=checkParField(pars,'deltat','double',c(),c(1,1),0,c(),verbose);allOk=pars$ok&&allOk;
      }
    }
    #has to be checked in relation to settings
    pars=checkParField(pars,'profile','double',as.double(1),c(),c(),c(),verbose);allOk=pars$ok&&allOk;
  }
  
  pars$ok=allOk;  
  return(pars);
}

IDMDist<-function(IDMPars,settings,verbose=2){
  #default output
  out=list();
  out$distributions=array(0,dim=c(0,0));
  out$ok=FALSE;
  out$extra=data.frame();
  out$used=list();
  
  #settings
  settings=checkSettings(settings,verbose);
  if (!(settings$ok)){
    return(out);
  }
  
  #IDMPars
  IDMPars=checkIDMPars(IDMPars,verbose);
  if (!(IDMPars$ok)){
    return(out);
  }
  
  #transform to general parameters
  pars=list();
  #if settings$seed=0, then set time-dependent seed
  if (settings$seed==0){  
    settings$seed=1+as.integer(2147483647*as.double(format(Sys.time(), format="%OS3"))/60);
  }
  set.seed(settings$seed-1);
  pars$nDim=IDMPars$nDim;
  pars$nStimuli=IDMPars$nStimuli;
  pars$N=IDMPars$N;
  pars$B=IDMPars$B;
  pars$eta=matrix(0,pars$nDim,pars$nStimuli);  
  pars$beta=IDMPars$beta;
  pars$W=IDMPars$W;
  pars$Theta=IDMPars$Theta;
  pars$spontaneousTime=IDMPars$spontaneousTime;
  pars$dynamics=IDMPars$dynamics;
  if(pars$dynamics==0){
    pars$D=IDMPars$D;
  }else if(pars$dynamics==1){
    pars$deltat=IDMPars$deltat;
    pars$sigma=IDMPars$sigma;
  }
  pars$startPos=matrix(rep(IDMPars$startPos,settings$nWalkers),pars$nDim,settings$nWalkers);
  pars$boxShape=IDMPars$boxShape;
  pars$h=IDMPars$h;
  if(IDMPars$sTer==0){
    pars$Ter=matrix(IDMPars$Ter,1,settings$nWalkers);
  }else{
    pars$Ter=IDMPars$Ter+IDMPars$sTer*(runif(settings$nWalkers)-0.5);
  }
  #extend profile if a constant
  if (length(IDMPars$profile)==1){
    pars$profile=matrix(IDMPars$profile,settings$nBins,1);
  }else{
    pars$profile=IDMPars$profile;
  }
  #check if profile is compatible with settings
  pars=checkParField(pars,'profile','double',c(),c(settings$nBins,1),c(),c(),verbose);
  if (!(pars$ok)){
    return(out);
  }
  
  if(pars$dynamics==0){
    out=RTDist(pars$nDim,pars$dynamics,pars$startPos,pars$eta,pars$nStimuli,pars$beta,pars$N,pars$W,pars$Theta,pars$B,pars$profile,pars$D,pars$h,pars$boxShape,pars$spontaneousTime,pars$Ter,settings$dt,settings$nBins,settings$binWidth,settings$seed,settings$nWalkers,settings$nGPUs,settings$gpuIds,settings$loadPerGPU);
  }else if(pars$dynamics==1){
    out=RTDist(pars$nDim,pars$dynamics,pars$startPos,pars$eta,pars$nStimuli,pars$beta,pars$N,pars$W,pars$Theta,pars$B,pars$profile,pars$sigma,pars$h,pars$boxShape,pars$spontaneousTime,pars$Ter,pars$deltat,settings$nBins,settings$binWidth,settings$seed,settings$nWalkers,settings$nGPUs,settings$gpuIds,settings$loadPerGPU);
  }
  
  #walkersAccountedFor check
  ok=as.logical(min(out$extra$walkersAccountedFor==settings$nWalkers));
  if (!ok){
    if (verbose>=1){
      print('WARNING: not all walkers accounted for');
    }
  }
  out$ok=ok;
  
  #add used
  out$used$IDMPars=IDMPars;
  out$used$pars=pars;
  out$used$settings=settings;
  return(out);
}
