#include <RcppArmadillo.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include<omp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;

// [[Rcpp::export]]

void upc(cube &beta, mat &Xmat, const mat &basismat, const cube Y, const mat &bmat, const double sig, const int Mb) {
  int itr = 0;
  int G = Xmat.n_cols;
  
  int n1 = Y.n_rows;
  int n2 = Y.n_cols;
  
  uvec ind = regspace<uvec>(1,  1,  G-1);
  colvec pvec;
  colvec diff;
  int M = Y.n_slices;
  
  mat crX1 = zeros<mat>(Mb*G, Mb*G);
  
  mat gen1 = zeros<mat>(G*Mb, M);
  
  colvec d    = sig*ones<vec>(Mb*G);//(1-pvec)/nu0(g) + pvec/nu1(g);
  mat X0 = zeros<mat>(M, Mb*G);
  for(int g = 0; g < G; g++){
    uvec indg = regspace<uvec>(g*Mb,  1,  (g+1)*Mb-1);
    
    X0.cols(indg) = diagmat(bmat.col(g)) * basismat;
    
  }
  
  crX1 = X0.t()*X0;
  gen1 = inv(crX1+diagmat(d))*X0.t();
  
  
  int i;
  int j;
  
  
  for(i = 0;i <n1; i++){
    for(j=0;j<n2;j++){
      colvec Yij = Y.subcube(span(i), span(j), span());
      if(mean(abs(Yij))>0){
        colvec gen = gen1 * Yij;
        
        beta.subcube(span(i), span(j), span())  = gen;
      }
    }
  }
}
