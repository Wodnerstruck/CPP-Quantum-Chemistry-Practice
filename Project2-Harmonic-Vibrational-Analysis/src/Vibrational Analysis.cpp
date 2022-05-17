#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdio>

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
    mat A;
    A.load("./H2O_HS.txt");
    mat hess = A.t();
    hess.reshape(9,9);
    hess.print();

    mat m = {{16,16,16,1,1,1,1,1,1}};
    mat mass = sqrt(kron(m.t(),m));
    
    mat mass_weighted_hess = hess / mass;
    vec eig_val;
    mat eig_vec;
    eig_sym(eig_val,eig_vec,mass_weighted_hess);
    eig_val.print();


    return 0;
}
