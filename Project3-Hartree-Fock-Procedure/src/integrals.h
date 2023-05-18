#ifndef INTEGRALS_H
#define INTEGRALS_H



#include <armadillo>

#include "index.h"
#include <vector>

using namespace arma;

class integrals
{
public:
    int dim = 0;//nabasis
    double E_NN;
    arma::mat S;
    arma::mat T;
    arma::mat V;
    arma::mat HC;//Hcore
    std::vector<double> ERI;
    arma::mat X;
    integrals( std::string filename);
    ~integrals();
};


#endif
