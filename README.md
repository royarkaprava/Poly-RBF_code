# Poly-RBF_code

Application.R has an example implementation using the shared data.
HCP protocol specifications are provided for Harmonization purposes in HCP_protocol folder

DWI_Processing.R has the main implementation function 'PolyRBFgivenN_K' which depends on the Rcpp code multitry1.cpp

Codes associated with cross-validation and an example application are also provided as:

poly_rbf_rev.R and poly_rbf_rev_cv.R are the new codes to implement the methods with cross-validation.

Application_with_CV.R provides an example implementation of the method
