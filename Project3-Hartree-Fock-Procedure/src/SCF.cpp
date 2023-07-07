#include "SCF.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <armadillo>
using namespace arma;
//initial guess H_core
SCF::SCF(const integrals &Int)
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
void SCF::iterat( const integrals &Int,std::string M_ ){
    if(M_ == "MP2" || M_ == "mp2")
        Is_mp2 = 1;
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
    //MP2
    if (Is_mp2 == 1){
        cout << "MP2 Calculation : " << endl;
        // vector构造n*n*n*n的四维数组
        arma::field<arma::mat> FD(Int.dim,Int.dim);
        for (int i = 0; i < Int.dim; i++)
        {
            for (int j = 0; j < Int.dim; j++)
            {
                FD(i,j) = arma::zeros(Int.dim,Int.dim);
            }   
        }
        arma::field<arma::mat>MO_1 = FD;
        for(int i=0; i < Int.dim; i++)
            for(int j=0; j < Int.dim; j++) {
                for(int k=0; k < Int.dim; k++)
                    for(int l=0; l < Int.dim; l++) {
                        int ij = INDEX(i,j);
                        int kl = INDEX(k,l);
                        int ijkl = INDEX(ij,kl);
                        MO_1(i,j)(k,l) += Int.ERI[ijkl];
                    } 
                } 
        // 4th index  
        arma::field<arma::mat>MO_2 = FD;
        for (int i = 0; i < Int.dim; i++)
        {
            for (int j = 0; j < Int.dim; j++)
            {
                for (int k = 0; k < Int.dim; k++)
                {
                    for (int l = 0; l < Int.dim; l++)
                    {
                        for (int m =0; m < Int.dim; m++)
                            MO_2(i,j)(k,l) += MO_1(i,j)(k,m) * C(m,l);
                    }
                    
                }
            }
            
        }
        //3rd index
        arma::field<arma::mat>MO_3 = FD;
        for (int i = 0; i < Int.dim; i++)
        {
            for (int j = 0; j < Int.dim; j++)
            {
                for (int k = 0; k < Int.dim; k++)
                {
                    for (int l = 0; l < Int.dim; l++)
                    {
                        for (int m =0; m < Int.dim; m++)
                            MO_3(i,j)(k,l) += MO_2(i,j)(m,l) * C(m,k);
                    }
                    
                }
            }   
        }
        //2nd index
        arma::field<arma::mat>MO_4 = FD;
        for (int i = 0; i < Int.dim; i++)
        {
            for (int j = 0; j < Int.dim; j++)
            {
                for (int k = 0; k < Int.dim; k++)
                {
                    for (int l = 0; l < Int.dim; l++)
                    {
                        for (int m =0; m < Int.dim; m++)
                            MO_4(i,j)(k,l) += MO_3(i,m)(k,l) * C(m,j);
                    }
                    
                }
            }   
        }
        //1st index
        arma::field<arma::mat>ERI_MO = FD;
        for (int i = 0; i < Int.dim; i++)
        {
            for (int j = 0; j < Int.dim; j++)
            {
                for (int k = 0; k < Int.dim; k++)
                {
                    for (int l = 0; l < Int.dim; l++)
                    {
                        for (int m =0; m < Int.dim; m++)
                            ERI_MO(i,j)(k,l) += MO_4(m,j)(k,l) * C(m,i);
                    }
                    
                }
            }   
        }
       
        //MP2 Energy
        int ndocc = 5;
        for (int i = 0; i < ndocc; i++)
            for (int a = ndocc; a < Int.dim; a++)
                for (int j = 0; j < ndocc; j++)
                    for (int b = ndocc; b < Int.dim; b++)
                        E_MP2 += ERI_MO(i,a)(j,b) * (2.0 * ERI_MO(i,a)(j,b) - ERI_MO(i,b)(j,a)) / (eps_0(i) + eps_0(j) - eps_0(a) - eps_0(b));
        cout << "MP2 Energy : " << std::fixed << std::setprecision(12) << E_MP2 << " " << endl;
        cout << "Total Energy : " << std::fixed << std::setprecision(12) << E_2 + (Int.E_NN) + E_MP2 << " " << endl;


    }        
        
        
        
}
    
