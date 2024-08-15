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
out <- PolyRBF(data_path=datapath, mask=mask,
                      local_bvals_path,local_bvecs_path,
                      pred_bvals_path=pred_bvals_path,pred_bvecs_path=pred_bvecs_path,
                      order=K,Mb=N, sig = 0.01)

#Fitted signal on its own grid of bvals and bvecs
#Here mask=NULL, which will prompt the method to be applied on the whole brain
out <- PolyRBF(data_path=datapath, mask=NULL,
                      local_bvals_path,local_bvecs_path,
                      pred_bvals_path=NULL,pred_bvecs_path=NULL,
                      order=K,Mb=N, sig = 0.01)


data=fast_readnii(datapath)
plot(density(data-out$DATApred))