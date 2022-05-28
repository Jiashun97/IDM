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
#   $Id: LDMWrapper.R 89 2015-06-14 09:40:30Z u0066818@kuleuven.be $

checkLDMPars <- function(pars,verbose=2){
  pars=as.list(pars);
  allOk=TRUE;  
  
  pars=checkParField(pars,'nStimuli','integer',c(),c(1,1),1,c(),verbose);allOk=pars$ok&&allOk;
  if(allOk){
    nStimuli=pars$nStimuli;
    pars=checkParField(pars,'c','double',as.double(0.1),c(1,1),0,c(),verbose);allOk=pars$ok&&allOk;
    pars=checkParField(pars,'a','double',c(),c(1,1),0,0.99,verbose);allOk=pars$ok&&allOk;
    pars=checkParField(pars,'zr','double',as.double(0.5),c(1,1),0,1,verbose);allOk=pars$ok&&allOk;
    if (allOk){
      pars=checkParField(pars,'sz','double',as.double(0),c(1,1),0,2*pars$a*max(pars$zr,1-pars$zr),verbose);allOk=pars$ok&&allOk;
    }
    pars=checkParField(pars,'Ter','double',as.double(0),c(1,1),0,c(),verbose);allOk=pars$ok&&allOk;
    if (allOk){
      pars=checkParField(pars,'sTer','double',as.double(0),c(1,1),0,2*pars$Ter,verbose);allOk=pars$ok&&allOk;
    }
    pars=checkParField(pars,'gamma','double',as.double(0),c(1,1),c(),c(),verbose);allOk=pars$ok&&allOk;#for OU    
    pars=checkParField(pars,'v','double',c(),c(1,nStimuli),c(),c(),verbose);allOk=pars$ok&&allOk;
    pars=checkParField(pars,'eta','double',matrix(0,1,nStimuli),c(1,nStimuli),0,c(),verbose);allOk=pars$ok&&allOk;
    #has to be checked in relation to settings
    pars=checkParField(pars,'profile','double',as.double(1),c(),c(),c(),verbose);allOk=pars$ok&&allOk;
  }

  pars$ok=allOk;
  return(pars);
}

LDMDist<-function(LDMPars,settings,verbose=2){
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
  
  #LDMPars
  LDMPars=checkLDMPars(LDMPars,verbose);
  if (!(LDMPars$ok)){
    return(out);
  }
  
  #transform to general parameters
  pars=list();
  #if settings$seed=0, then set time-dependent seed
  if (settings$seed==0){  
    settings$seed=1+as.integer(2147483647*as.double(format(Sys.time(), format="%OS3"))/60);
  }
  set.seed(settings$seed-1);
  pars$boxShape=as.integer(3);
  pars$nDim=as.integer(1);
  pars$dynamics=as.integer(0);
  pars$N=as.double(0);
  pars$beta=as.double(1);
  pars$Theta=as.double(0);
  pars$spontaneousTime=as.double(0);
  
  pars$D=as.double(0.5*LDMPars$c^2);
  startPos=0.5+(LDMPars$zr-0.5)*LDMPars$a;
  if (LDMPars$sz==0){
    pars$startPos=startPos*matrix(1,1,settings$nWalkers);
  }else{
    pars$startPos=startPos+LDMPars$sz*(runif(settings$nWalkers)-0.5);    
  }
  pars$W=as.double(0.5*LDMPars$gamma)/pars$D;#OU
  pars$nStimuli=as.integer(length(LDMPars$v));
  pars$eta=as.double(LDMPars$eta)/pars$D;
  pars$h=matrix(c(0.5*(1-LDMPars$a),0.5*(1-LDMPars$a)),1,2);
  pars$B=as.double(LDMPars$v-LDMPars$gamma*pars$h[1])/pars$D;
  if(LDMPars$sTer==0){
    pars$Ter=matrix(LDMPars$Ter,1,settings$nWalkers);
  }else{
    pars$Ter=LDMPars$Ter+LDMPars$sTer*(runif(settings$nWalkers)-0.5);
  }
  #extend profile if a constant
  if (length(LDMPars$profile)==1){
    pars$profile=matrix(LDMPars$profile,settings$nBins,1);
  }else{
    pars$profile=LDMPars$profile;
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
  out$used$LDMPars=LDMPars;
  out$used$pars=pars;
  out$used$settings=settings;
  return(out);
}
