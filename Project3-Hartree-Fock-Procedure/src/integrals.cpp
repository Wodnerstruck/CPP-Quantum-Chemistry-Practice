#include "integrals.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <armadillo>
using namespace arma;

integrals::integrals( std::string filename){
    //get nbasis
    if (filename == "./data/h2o/STO-3G")
        dim = 7;
        else if (filename == "./data/h2o/DZ")
            dim = 14;
            else
                dim = 26;
    
    cout<<"The number of basis functions is : "<< dim << endl;

    //build matrices
    S = mat(dim,dim,arma::fill::zeros);
    T = S;
    V = S;
    ERI.resize(((dim*dim+dim)/2)*(((dim*dim+dim)/2)+1)/2,0);
     
    std::ifstream is ;
    is.open(filename + ("/enuc.dat"));
    is >> E_NN;
    cout<<"E_NN stored"<<endl;
    cout<< E_NN << endl;
    is.close();
    is.clear();
    is.open(filename + ("/s.dat"));
     cout << is.is_open()<<endl;
    for (int i = 0; i < (dim*dim+dim)/2; i++)
    {
        int j,k;
        is >> j >> k;
        is>>S(j-1,k-1);
        S(k-1,j-1) = S(j-1,k-1);
    }
    cout<< "S Matrix stored"<<endl;
    is.close();
    is.clear();
    is.open(filename + ("/t.dat") );
    for (int i = 0; i < (dim*dim+dim)/2; i++)
    {
        int j,k;
        is>>j>>k>>T(j-1,k-1);
        T(k-1,j-1) = T(j-1,k-1 );
    }
    is.close();
    is.clear();
    is.open(filename + ("/v.dat") );
    for (int i = 0; i < (dim*dim+dim)/2; i++)
    {
        int j,k;
        is>>j>>k>>V(j-1,k-1);
        V(k-1,j-1) = V(j-1,k-1);
    }
    is.close();
    is.clear();
    HC = T + V;
    is.open(filename + ("/eri.dat")) ;
    while (is.peek()!= EOF)
    {
        int i,j,k,l;
        is>>i>>j>>k>>l;
        i-=1;
        j-=1;
        k-=1;
        l-=1;
        int ij = INDEX(i,j);
        int kl = INDEX(k,l);
        int ijkl = INDEX(ij,kl);
        is>>ERI[ijkl];

    }
    is.close();
    //build orthogonalization matrix X
    mat eig_vec;
    vec eig_val;
    arma::eig_sym(eig_val,eig_vec,S);//diagnalzation of S matrix: s=L^T * S * L
    eig_val = arma::pow(eig_val,-0.5);//s^(-1/2)
    X = eig_vec * arma::diagmat(eig_val) * eig_vec.t();//build X
    
}
integrals::~integrals(){

}