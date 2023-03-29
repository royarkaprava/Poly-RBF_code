DWI_Processing<-function(data_path,
                         local_bvals_path,local_bvecs_path,
                         pred_bvals_path=NULL,pred_bvecs_path=NULL,
                         order,Mb=10, sig = 0.01){
  

  
  #load the cpp file
  Rcpp::sourceCpp("multitry1.cpp")
  
  ### Contributors : ###########################################################
  ### Writer: Arkaprava Roy
  ### Wrap-upper : Zhou Lan
  ### Other Writers: Zhengwu Zhang
  ##############################################################################
  
  ### data_path: the data of diffusion weighted signals
  ### local_bvals_path, local_bvecs_path: the b-values and b-vectors of the original diffusion weighted image
  ### pred_bvals_path, pred_bvecs_path: the b-values and b-vectors of the predicted diffusion weighted image; 
  ###                         if they are set as NULL, we use the local values instead.
  ### order: the polynomial order
  ### Mb=10, sig=0.01 tuning parameters

  
  #### Loading required R packages #############################################
  library(neurobase)
  ##############################################################################
  
  
  
  
  
  
  #### I. Preparing Steps: Reading the Input data ##############################
  ## 0. Checking Argument
  if(is.null(pred_bvals_path)){
    pred_bvals_path<-local_bvals_path
  }
  
  if(is.null(pred_bvecs_path)){
    pred_bvecs_path<-local_bvecs_path
  }
  
  ## 1. Read the data
  data=fast_readnii(data_path)
  
  ## 2. Read b-values and b-vectors
  bval <- read.table(local_bvals_path, quote="\"", comment.char="")
  p <- read.table(local_bvecs_path, quote="\"", comment.char="")
  bval <- as.numeric(unlist(bval))
  
  bvalHCP <- read.table(pred_bvals_path, quote="\"", comment.char="")
  pHCP <- read.table(pred_bvecs_path, quote="\"", comment.char="")
  bvalHCP <- as.numeric(unlist(bvalHCP))
  ##############################################################################
  
  
  ### II. Defining Important Parameters ########################################
  ind0HCP <- which(bvalHCP==min(bvalHCP))
  pts=pHCP; bvalts=bvalHCP; 
  
  if(length(ind0HCP)){
    bvalts <- bvalHCP[-ind0HCP]
    pts <- pHCP[, -ind0HCP]
  }
  
  bval <- as.numeric(unlist(bval))
  Mbval <- max(bval)
  bval <- bval/max(bval) 
  M <- length(bval)
  intercept = F
  ind0 <- which(bval==min(bval))
  dataf0 <- data[,,,ind0]
  data0 <- apply(dataf0, 1:3, mean)
  
  if(length(ind0)){
    bval <- bval[-ind0]
    p <- p[, -ind0]
  }
  ##############################################################################
  
  ### III. Modification of the MRI data ########################################
  data <- data[,,,-ind0]
  datar <- apply(data, 4, FUN=function(x){(x+1)/data0})
  if(sum(is.na(datar))) datar[which(is.na(datar))] <- 0
  datar <- array(datar, dim=dim(data))
  remove(data)
  data <- log(datar)
  if(sum(is.na(data))) data[which(is.na(data))] <- 0
  pts=pts; bvalts=bvalts/Mbval
  M <- length(bval)
  ##############################################################################
  
  ### IV. Processing ###########################################################
  print("Fitting starts now")
  bmat <- rep(1, M)#bval
  for(i in 1:order){
    bmat <- cbind(bmat, bval^i)
  }
  if(!intercept){
    bmat <- bmat[, -1]
  }
  
  bmat <- matrix(bmat, nrow=M)
  #Mb <- 20
  
  #M <- dim(data)[1]
  G <- ncol(bmat)
  Xmat <- matrix(0, M, G)#Y/bmat#Y/bmat#
  
  Mts <- length(bvalts)
  bmatts <- rep(1, Mts)#bvalts
  for(i in 1:order){
    bmatts <- cbind(bmatts, bvalts^i)
  }
  if(!intercept){
    bmatts <- bmatts[, -1]
  }
  
  bmatts <- matrix(bmatts, nrow=Mts)
  
  fibonacci_sphere<- function(samples=1000){
    points = matrix(0, 3, samples)
    phi = pi * (3. - sqrt(5.))  # golden angle in radians
    
    for(i in 1:samples){
      y = 1 - ((i-1) / (samples - 1)) * 2  # y goes from 1 to -1
      radius = sqrt(1 - y * y)  # radius at y
      
      theta = phi * (i-1)  # golden angle increment
      
      x = cos(theta) * radius
      z = sin(theta) * radius
      
      points[, i]=c(x, y, z)}
    
    return(points)
  }
  
  pmat <- fibonacci_sphere(Mb)#matrix(0, 3, Mb)
  
  ind <- combinat::combn(ncol(pmat), 2)
  f   <- lapply(1:ncol(ind), function(k){sum(pmat[, ind[1, k]]*pmat[, ind[2, k]])})
  indR <- ind[, which(abs(unlist(f)+1)==0)]
  
  if(length(indR)){
    indRt <- 1:(length(indR)/2)+1
    pmat  <- pmat[, -indRt]  
  }
  
  pmat <- cbind(pmat, -pmat)
  
  Mb <- ncol(pmat)
  
  Mbmat <- rbind(diag(Mb/2), diag(Mb/2))
  
  dis <- as.matrix(dist(t(cbind(p, pmat))))
  dis <- dis[(1:M),(1:Mb+M)]
  
  bd <- sqrt(2)*mean(dis[dis>0])
  
  basismat <- (dis<=(3*bd/sqrt(2)))*exp(-dis^2/bd^2)
  
  basismat <- basismat %*% Mbmat
  
  dis <- as.matrix(dist(t(cbind(pts, pmat))))
  dis <- dis[(1:Mts),(1:Mb+Mts)]
  
  basismatts <- (dis<=(3*bd/sqrt(2)))*exp(-dis^2/bd^2)
  
  basismatts <- basismatts %*% Mbmat
  
  process <- function(Y){
    Mb <- ncol(basismat)
    
    beta <- matrix(0, Mb, G)
    
    Xmat <- basismat %*% beta
    
    muts <- array(0, dim=c(dim(Y)[1:2], Mts))
    mu <- array(0, dim=c(dim(Y)[1:2], M))
    
    upc(beta, Xmat, muts, mu, basismat, basismatts, Y, bmat, bmatts, sig)
    
    return(muts)
  }
  #Fit the model
  out <- apply(data, 1, process)
  
  out2 <- array(t(out), dim=c(dim(data)[1:3],Mts))
  
  dimo <- dim(out2)[1:3]
  out2 <- apply(out2, 4, FUN=function(x){exp(x)*data0})
  out2 <- array(out2, dim=c(dimo,Mts))
  
  bval1 <- read.table(local_bvals_path, quote="\"", comment.char="")
  p1 <- read.table(local_bvecs_path, quote="\"", comment.char="")
  
  bval1 <- as.numeric(unlist(bval1))
  bval2 <- bval1/max(bval1)
  
  #ind0 <- which(bval2==min(bval2))
  
  #out2 <- array(0, dim=c(dim(data)[1:3],Mts+length(ind0)))
  ##############################################################################
  
  
  
  #### Return Object ###########################################################
  outf <- array(0, dim=c(dimo, dim(out2)[4]+length(ind0)))
  outf[,,,1:dim(out2)[4]] <- out2#array(t(out), dim=c(dim(data)[1:3],Mts))
  outf[,,,1:length(ind0)+dim(out2)[4]] <- dataf0
  
  pf    <- cbind(pts, p1[, ind0])
  bvalf <- c(bvalts*Mbval, bval1[ind0])
  
  if(sum(is.nan(outf))) outf[which(is.nan(outf))] <- 0
  
  
  #predicted signal, bvecs and bvals in this order
  return(list(DATApred=outf, bvecs=pf, bvals=bvalf))
  ##############################################################################
}