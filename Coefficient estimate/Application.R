##Installation of neurobase which is needed
install.packages("devtools")
devtools::install_github("muschellij2/neurobase")

#This is the reduced cfin_multib data to be able to push in Github
datapath <- "data/cfin_multib/datareduced.nii.gz"

local_bvals_path = "data/cfin_multib/bval.bval"

local_bvecs_path <- "data/cfin_multib/bvec.bvec"

#Creating a mask
mask <- array(0, dim=c(31,  41,  19)) #Dimension of the data in datapath
mask[10:20, 10:30, 5:12] <- 1


#Estimated Coef for each voxel
K=4
N=10
out <- PolyRBFcoef(data_path=datapath, mask=mask,
               local_bvals_path,local_bvecs_path,
               order=K,Mb=N, sig = 0.01)
