source("poly_rbf_rev.R")
source("poly_rbf_rev_cv.R")

############################Cross-validation based fitting#########
if(is.null(pred_bvals_path)){
  pred_bvals_path<-local_bvals_path
}

if(is.null(pred_bvecs_path)){
  pred_bvecs_path<-local_bvecs_path
}

# 1. Read the data
data=fast_readnii(data_path)

if(!is.null(mask)){
  data <- array(apply(data, 4, FUN=function(x){x*mask}), dim(data))
}

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
#pts=pHCP; bvalts=bvalHCP; 

if(length(ind0HCP)){
  bvalts <- bvalHCP[-ind0HCP]
  pts <- pHCP[, -ind0HCP]
}

fit <- PolyRBFCV(data, bval, p,  order = c(2,3,4,5,6),Mb=c(5,10,15), bdls, sig = 0.01)

out <- PolyRBF(data, bval,p, bvalts,  pts,  order = fit[2],Mb=fit[1], sig = 0.01, bd=NA)

