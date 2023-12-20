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
#   $Id: LCAWrapper.R 69 2015-03-30 16:35:24Z u0066818@kuleuven.be $

checkLCAPars <- function(pars,verbose=2){
  pars=as.list(pars);
  allOk=TRUE;  
  
  pars=checkParField(pars,'nDim','integer',c(),c(1,1),1,c(),verbose);allOk=pars$ok&&allOk;
  pars=checkParField(pars,'nStimuli','integer',c(),c(1,1),1,c(),verbose);allOk=pars$ok&&allOk;
  if (allOk){
    nDim=pars$nDim;
    nStimuli=pars$nStimuli;
    pars=checkParField(pars,'c','double',matrix(0.1,nDim,1),c(nDim,1),0,c(),verbose);allOk=pars$ok&&allOk;
    pars=checkParField(pars,'startPos','double',matrix(0.0001,nDim,1),c(nDim,1),0.0001,0.9999,verbose);allOk=pars$ok&&allOk;
    pars=checkParField(pars,'Gamma','double',c(),c(nDim,nDim),c(),c(),verbose);allOk=pars$ok&&allOk;
    pars=checkParField(pars,'a','double',c(),c(nDim,1),0.0001,0.9999,verbose);allOk=pars$ok&&allOk;
    pars=checkParField(pars,'v','double',c(),c(nDim,nStimuli),c(),c(),verbose);allOk=pars$ok&&allOk;
    pars=checkParField(pars,'Ter','double',as.double(0),c(1,1),0,c(),verbose);allOk=pars$ok&&allOk;
    if (allOk){
      pars=checkParField(pars,'sTer','double',as.double(0),c(1,1),0,2*pars$Ter,verbose);allOk=pars$ok&&allOk;
    }
    #has to be checked in relation to settings
    pars=checkParField(pars,'profile','double',as.double(1),c(),c(),c(),verbose);allOk=pars$ok&&allOk;
  }
  
  pars$ok=as.integer(allOk);  
  return(pars);
}

LCADist<-function(LCAPars,settings,verbose=2){
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
  
  #LCAPars
  LCAPars=checkLCAPars(LCAPars,verbose);
  if (!(LCAPars$ok)){
    return(out);
  }
  
  #transform to general parameters
  pars=list();
  #if settings$seed=0, then set time-dependent seed
  if (settings$seed==0){  
    settings$seed=1+as.integer(2147483647*as.double(format(Sys.time(), format="%OS3"))/60);
  }
  set.seed(settings$seed-1);
  pars$boxShape=as.integer(0);
  pars$dynamics=as.integer(0);
  pars$beta=as.double(1);
  pars$spontaneousTime=as.double(0);
  
  pars$nDim=LCAPars$nDim;
  pars$N=matrix(0,pars$nDim,1);
  pars$Theta=matrix(0,pars$nDim,1);
  pars$nStimuli=LCAPars$nStimuli;
  pars$D=0.5*LCAPars$c^2;
  pars$B=LCAPars$v/matrix(rep(pars$D,pars$nStimuli),pars$nDim,pars$nStimuli);;
  pars$eta=matrix(0,pars$nDim,pars$nStimuli);  
  #convert Gamma to W (with transpose)
  pars$W=t(-LCAPars$Gamma/matrix(rep(pars$D,pars$nDim),pars$nDim,pars$nDim));
  diag(pars$W)=0.5*diag(LCAPars$Gamma)/pars$D;
  pars$startPos=matrix(rep(LCAPars$startPos,settings$nWalkers),pars$nDim,settings$nWalkers);
  pars$h=matrix(1,pars$nDim,pars$nDim);
  diag(pars$h)=1-LCAPars$a;

  if(LCAPars$sTer==0){
    pars$Ter=matrix(LCAPars$Ter,1,settings$nWalkers);
  }else{
    pars$Ter=LCAPars$Ter+LCAPars$sTer*(runif(settings$nWalkers)-0.5);
  }
  #extend profile if a constant
  if (length(LCAPars$profile)==1){
    pars$profile=matrix(LCAPars$profile,settings$nBins,1);
  }else{
    pars$profile=LCAPars$profile;
  }
  #check if profile is compatible with settings
  pars=checkParField(pars,'profile','double',c(),c(settings$nBins,1),c(),c(),verbose);
  if (!(pars$ok)){
    return(out);
  }
  
  out=RTDist(pars$nDim,pars$dynamics,pars$startPos,pars$eta,pars$nStimuli,pars$beta,pars$N,pars$W,pars$Theta,pars$B,pars$profile,pars$D,pars$h,pars$boxShape,pars$spontaneousTime,pars$Ter,settings$dt,settings$nBins,settings$binWidth,settings$seed,settings$nWalkers,settings$nGPUs,settings$gpuIds,settings$loadPerGPU);
  
  #walkersAccountedFor check
  ok=as.logical(min(out$extra$walkersAccountedFor==settings$nWalkers));
  if (!ok){
    if (verbose>=1){
      print('WARNING: not all walkers accounted for');
    }
  }
  out$ok=ok;
  
  #add used
  out$used$LCAPars=LCAPars;
  out$used$pars=pars;
  out$used$settings=settings;
  return(out);
}
