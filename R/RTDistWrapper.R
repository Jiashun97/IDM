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
#   $Id: RTDistWrapper.R 69 2015-03-30 16:35:24Z u0066818@kuleuven.be $

#load dynamic linked library
if (.Platform$OS.type == "windows") {
  dyn.load("RTDist.dll");
}else{
  if(substr(.Platform$pkgType,1,3)=="mac"){
    dyn.load("libiomp5.dylib");
    dyn.load("RTDist.bundle");
  }else{
    dyn.load("libiomp5.so");
    dyn.load("RTDist.so");
  }
}

#declare function 
RTDist <- function(nDim,dynamics,startPos,eta,nStimuli,beta,N,W,Theta,B,profile,D,h,boxShape,spontaneousTime,Ter,dt,nBins,binWidth,seed,nWalkers,nGPUs,gpuIds,loadPerGPU){
  extraLength=6;
  nBoxes=nDim;
  if (boxShape==3){
    if (nDim!=1){
      stop("boxShape=3 can only be used with nDim=1 (2 bound diffusion model)");
    }
    nBoxes=2;
  }  
  if (nGPUs>0){
    output <- .C("getFirstPassageTimeDistributionsGPUs",
                 distributions=integer(nStimuli*nBoxes*nBins),
                 extra=integer(nStimuli*extraLength), 
                 as.integer(nDim), 
                 as.integer(dynamics), 
                 as.single(startPos), 
                 as.single(eta),
                 as.integer(nStimuli),
                 as.single(beta), 
                 as.single(N),
                 as.single(W), 
                 as.single(Theta),
                 as.single(B), 
                 as.single(profile),
                 as.single(D),
                 as.single(h), 
                 as.integer(boxShape), 
                 as.single(spontaneousTime), 
                 as.single(Ter),
                 as.single(dt),
                 as.integer(nBins),
                 as.single(binWidth),
                 as.integer(seed),
                 as.integer(nWalkers),
                 as.integer(nGPUs),
                 as.integer(gpuIds),
                 as.integer(loadPerGPU)
    )
  }else{
    output <- .C("getFirstPassageTimeDistributionsCPUs",
                 distributions=integer(nStimuli*nBoxes*nBins),
                 extra=integer(nStimuli*extraLength), 
                 as.integer(nDim), 
                 as.integer(dynamics), 
                 as.single(startPos), 
                 as.single(eta),
                 as.integer(nStimuli),
                 as.single(beta), 
                 as.single(N),
                 as.single(W), 
                 as.single(Theta),
                 as.single(B), 
                 as.single(profile),
                 as.single(D),
                 as.single(h), 
                 as.integer(boxShape), 
                 as.single(spontaneousTime), 
                 as.single(Ter),
                 as.single(dt),
                 as.integer(nBins),
                 as.single(binWidth),
                 as.integer(seed),
                 as.integer(nWalkers))
  }
  out=list();
  out$distributions=matrix(output$distributions,nBins);  
  
  flood=array();
  multiArrival=array();
  binsExceeded=array();
  lemmingsSpontaneous=array();
  lemmingsStimulus=array();
  walkersAccountedFor=array();
  for (i in 1:nStimuli){    
    flood[i]=output$extra[extraLength*(i-1)+1];
    multiArrival[i]=output$extra[extraLength*(i-1)+2];
    binsExceeded[i]=output$extra[extraLength*(i-1)+3];
    lemmingsSpontaneous[i]=output$extra[extraLength*(i-1)+4];
    lemmingsStimulus[i]=output$extra[extraLength*(i-1)+5];
    walkersAccountedFor[i]=output$extra[extraLength*(i-1)+6];  
  }
  out$extra=data.frame(flood,multiArrival,binsExceeded,lemmingsSpontaneous,lemmingsStimulus,walkersAccountedFor);
  
  return(out);
}

checkParField<-function(pars,field,className,default,dimensions,minValue,maxValue,verbose){
  ok=TRUE;
  if (field %in% names(pars)){
    if (storage.mode(pars[[field]])!=className){
      storage.mode(pars[[field]])=className;
      if (verbose>=2){
        print(paste('INFO:',field, 'converted to', className));
      }
    }    
    if (length(dimensions)==2){
      if (length(dim(pars[[field]]))==0){
        fieldDim=c(1,1);
      }else{
        fieldDim=dim(pars[[field]]);
      }  
      if (max(fieldDim!=dimensions)){
        if (verbose>=1){
          print(paste('ERROR:',field, 'not correctly dimensioned: [',fieldDim[1],',',fieldDim[2],'] should be [',dimensions[1],',',dimensions[2],']'));
        }
        ok=FALSE;
      }
    }
  }else{
    if (length(default)==0){
      if (verbose>=1){
        print(paste('ERROR:',field, 'required'));
      }
      ok=FALSE;
    }else{
      if (verbose>=2){
        print(paste('INFO:','using default',field));
      }
      pars[[field]]=default;
    }
  }
  #check min/max
  if (ok){
    if (length(minValue)>0){
      if (max(pars[[field]]<minValue)){
        if (verbose>=1){
          print(paste('ERROR:',field, 'too low (maybe in combination with other parameter values)'));
        }
        ok=FALSE;
      }
    }
    if (length(maxValue)>0){
      if (max(pars[[field]]>maxValue)){
        if (verbose>=1){
          print(paste('ERROR:',field, 'too high (maybe in combination with other parameter values)'));
        }
        ok=FALSE;
      }
    }
  }
  
  pars$ok=ok;    
  return(pars);
}

checkSettings <- function(settings,verbose=2){
  settings=as.list(settings);
  allOk=TRUE;
  
  settings=checkParField(settings,'nWalkers','integer',as.integer(100000),c(1,1),1,c(),verbose);allOk=settings$ok&&allOk;
  settings=checkParField(settings,'seed','integer',as.integer(0),c(1,1),c(),c(),verbose);allOk=settings$ok&&allOk;
  settings=checkParField(settings,'dt','double',as.double(0.001),c(1,1),0,c(),verbose);allOk=settings$ok&&allOk;
  settings=checkParField(settings,'nBins','integer',as.integer(3000),c(1,1),1,c(),verbose);allOk=settings$ok&&allOk;
  settings=checkParField(settings,'binWidth','double',as.double(0.001),c(1,1),0,c(),verbose);allOk=settings$ok&&allOk;
  settings=checkParField(settings,'nGPUs','integer',as.integer(1),c(1,1),0,c(),verbose);allOk=settings$ok&&allOk;
  if (settings$ok&&(settings$nGPUs>0)){
    settings=checkParField(settings,'gpuIds','integer',matrix(0:(settings$nGPUs-1),1,settings$nGPUs),c(1,settings$nGPUs),0,c(),verbose);allOk=settings$ok&&allOk;
    settings=checkParField(settings,'loadPerGPU','double',matrix(rep(1,settings$nGPUs),1,settings$nGPUs),c(1,settings$nGPUs),0,c(),verbose);allOk=settings$ok&&allOk;
  }else{
    settings$gpuIds=matrix(0,0,0);
    settings$loadPerGPU=matrix(0,0,0);    
  }
  settings$ok=allOk;
  
  return(settings);
}
