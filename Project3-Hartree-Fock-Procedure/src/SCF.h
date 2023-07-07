#ifndef SCF_H
#define SCF_H

#include "integrals.h"
#include <string>
class SCF
{

public:
    arma::mat F_prime; //transformed Fock Matrix
    arma::mat F; //Fock matrix
    arma::mat C; //Coeff Matrix
    arma::mat P_1; //Density matrix former
    arma::mat P_2;//Density matrix last
    arma::vec eps_0; //orbital energy
    
    double E_1 = 0; //energy of former step;
    double E_2;//energy of last step
    bool Is_mp2 = 0;
    double E_MP2 = 0.0;
    SCF(const integrals &Int);
    ~SCF();
    void iterat( const integrals &Int,std::string M_ = "NMP2");

};







#endif