#include "SCF.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <armadillo>
using namespace arma;
//initial guess H_core
SCF::SCF(integrals Int)
{
    F = Int.HC;
    F_prime = Int.X.t() * F * Int.X; //F‘
    arma::mat C_pri; //C'
    arma::eig_sym(eps_0,C_pri,F_prime);//diag
    C = Int.X * C_pri;//C Matrix
    P_2 = C.cols(0,4) * C.cols(0,4).t();
    cout<<"Initial transformed Fock Matrix : "<<endl;
    F_prime.print();
    cout<<"Initial Coefficients :"<<endl;
    C.print();
    cout<<"Initial Density Matrix :"<<endl;
    P_2.print();
    E_2 = 0.0;
  
    for (int i = 0; i < (Int.dim); i++)
    {
        for (int j = 0; j < (Int.dim); j++)
        {
           
            E_2 += P_2(i,j) * (Int.HC(i,j) + F(i,j));
        }
        
    }
    cout << "Initial Electronic Energy : " << std::fixed 
        << std::setprecision(12) << E_2 << " (Hartree)" << endl;
    cout << "Initial Total Energy : " << std::fixed 
        << std::setprecision(12) << E_2 + (Int.E_NN) << " (Hartree)" <<endl;
    cout << "\n\n\n"; 
}

SCF::~SCF()
{
}
//iterate process
void SCF::iterat(integrals Int){
    int cycles = 0;
    while (true)
    {   
        
        E_1 = E_2;
        P_1 = P_2;
        //build new Fock Matrix
        F = Int.HC;
        for(int i=0; i < Int.dim; i++)
            for(int j=0; j < Int.dim; j++) {
                for(int k=0; k < Int.dim; k++)
                    for(int l=0; l < Int.dim; l++) {
                        int ij = INDEX(i,j);
                        int kl = INDEX(k,l);
                        int ijkl = INDEX(ij,kl);
                        int ik = INDEX(i,k);
                        int jl = INDEX(j,l);
                        int ikjl = INDEX(ik,jl);
  
                        F(i,j) += P_1(k,l) * (2.0 * Int.ERI[ijkl] - Int.ERI[ikjl]);
            }
        }
        
        //build new Density Matrix
        F_prime = Int.X.t() * F * Int.X;
        arma::mat C_pri;//C'
        arma::eig_sym(eps_0,C_pri,F_prime);//diag
        C = Int.X * C_pri;//New C
        P_2 = C.cols(0,4) * C.cols(0,4).t();//New Density Matrix
        //New Electronic Energy
        E_2 = 0.0;
        for (int i = 0; i < (Int.dim); i++)
    {
        for (int j = 0; j < (Int.dim); j++)
        {
           
            E_2 += P_2(i,j) * (Int.HC(i,j) + F(i,j));
        }
        
    }
    // Cycle number
    cycles++; 
    //IF Converge
    double delta_E = abs(E_2 - E_1); //Energy difference
    arma::mat Delta_P = P_2 - P_1; //DM difference
    double RMS_D = sqrt(arma::accu(arma::pow(Delta_P,2)));// DM Root Mean Square
    if ((delta_E < 1e-12) && (RMS_D < 1e-12)) //Convergence threshold
    {
        cout<< "SCF Done :   E = " << std::fixed  << std::setprecision(12) 
            << E_2 + (Int.E_NN) << "    A.U.  after " << cycles << " Cycles" << endl;
        break;
    }
        else
            continue;
    

    }
    
}