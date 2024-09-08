PolyRBFCV<-function(data_path, mask=NULL,
                    local_bvals_path,local_bvecs_path,
                  order=orderls,Mb=Mbls, bdls, sig = 0.01){
  
  testind <- sample(1:length(bval[which(bval>min(bval))]), length(bval[which(bval>min(bval))])/4)
  
  bvaltr <- bval[-which(bval>min(bval))[testind]]
  bvalts <- bval[which(bval>min(bval))[testind]]
  ptr    <- p[,-which(bval>min(bval))[testind]]
  pts    <- p[,which(bval>min(bval))[testind]]
  
  if(!is.na(bdls)){
    out <- array(0, dim=c(length(Mb), length(order), length(bdls)))
    
    for(i in 1:length(Mb)){
      for(j in 1:length(order)){
        for(k in 1:length(bdls))
          fit <- PolyRBF(data[,,,-which(bval>min(bval))[testind]], bvaltr,ptr, bvalts,  pts,  order = order[j],Mb=Mb[i], sig = 0.01, bd=bdls[k])
        out[i,j,k] <- mean((data[,,,which(bval>min(bval))[testind]] - fit$DATApred[,,,1:length(bvalts)])^2)
      }
    }
    
    vec <- which( out==min(out,na.rm=T) , arr.ind = T )[1,]
    
    ret <- c(Mb[vec[1]], order[vec[2]], bdls[vec[3]])
  }
  
  if(is.na(bdls)){
    out <- matrix(0, length(Mb), length(order))
    
    for(i in 1:length(Mb)){
      for(j in 1:length(order)){
        for(k in 1:length(bdls))
          fit <- PolyRBF(data[,,,-which(bval>min(bval))[testind]], bvaltr,ptr, bvalts,  pts,  order = order[j],Mb=Mb[i], sig = 0.01, bd=NA)
        out[i,j] <- mean((data[,,,which(bval>min(bval))[testind]] - fit$DATApred[,,,1:length(bvalts)])^2)
      }
    }
    
    vec = which(out==min(out), arr.ind = T)
    ret <- c(Mb[vec[1]], order[vec[2]])
  }
  
  
  return(ret)
  ##############################################################################
}