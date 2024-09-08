##Installation of neurobase which is needed
install.packages("devtools")
devtools::install_github("muschellij2/neurobase")

#This is the reduced cfin_multib data to be able to push in Github
datapath <- "data/cfin_multib/datareduced.nii.gz"

local_bvals_path = "data/cfin_multib/bval.bval"

local_bvecs_path <- "data/cfin_multib/bvec.bvec"

pred_bvals_path <- "HCP_protocol/bvals"

pred_bvecs_path <- "HCP_protocol/bvecs"

#Creating a mask
mask <- array(0, dim=c(31,  41,  19)) #Dimension of the data in datapath
mask[10:20, 10:30, 5:12] <- 1

#Predicted on HCP grid
K=4
N=20
out <- PolyRBFgivenN_K(data_path=datapath, mask=mask,
                      local_bvals_path,local_bvecs_path,
                      pred_bvals_path=pred_bvals_path,pred_bvecs_path=pred_bvecs_path,
                      order=K,Mb=N, sig = 0.01)

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


#Fitted signal on its own grid of bvals and bvecs
#Here mask=NULL, which will prompt the method to be applied on the whole brain
out <- PolyRBF(data_path=datapath, mask=NULL,
                      local_bvals_path,local_bvecs_path,
                      pred_bvals_path=NULL,pred_bvecs_path=NULL,
                      order=K,Mb=N, sig = 0.01)


data=fast_readnii(datapath)
plot(density(data-out$DATApred))