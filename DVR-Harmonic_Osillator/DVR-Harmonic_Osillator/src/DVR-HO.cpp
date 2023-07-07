#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
    double delta_x = 0.01;
    double L = 5.00;
    int dim = int(2 * L / delta_x) + 1 ;

    cout << "range: [" << -L << "," << L << "]" << endl;
    cout << "dx: " << delta_x << endl;

    arma::mat T = arma::mat(dim,dim,arma::fill::zeros);
    T.diag() += 0.5 * (2 / pow(delta_x,2));
    T.diag(1) -= 0.5 * (1 / pow(delta_x,2));
    T.diag(-1) -= 0.5 * (1 / pow(delta_x,2));
    
    arma::mat V = arma::mat(dim,dim,arma::fill::zeros);
    for(int i = 0;i < dim;i++){
        double x_i = (-1 * L) + (i * delta_x);
        V(i,i) = 0.5 * pow(x_i,2);
    }
    
    arma::mat H = T + V;
    arma::vec E;
    arma::mat Phi;
    arma::eig_sym(E,Phi,H);
    cout << "Eigen values:  ";
    cout << E(0) << "   " << E(1) << endl ;
    cout << "dim: " << dim << endl;
}