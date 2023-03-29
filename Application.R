##Installation of neurobase which is needed
install.packages("devtools")
devtools::install_github("muschellij2/neurobase")


datapath <- "data/cfin_multib/datareduced.nii.gz"

local_bvals_path = "data/cfin_multib/bval.bval"

local_bvecs_path <- "data/cfin_multib/bvec.bvec"

pred_bvals<- "HCP_protocol/bvals"

pred_bvecs="HCP_protocol/bvecs"

#Predicted on HCP grid
out <- DWI_Processing(data_path=datapath,
                      local_bvals_path,local_bvecs_path,
                      pred_bvals_path=pred_bvals,pred_bvecs_path=pred_bvecs,
                      order=4,Mb=10, sig = 0.01)

#Fitted signal on its own grid of bvals and bvecs
out <- DWI_Processing(data_path=datapath,
                      local_bvals_path,local_bvecs_path,
                      pred_bvals_path=NULL,pred_bvecs_path=NULL,
                      order=4,Mb=10, sig = 0.01)