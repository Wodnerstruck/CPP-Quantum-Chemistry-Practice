#include "molecule.h"
#include "masses.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

int main()
{
    Molecule mol("./data/acetaldehyde.dat", 0);
    cout << "Number of Atom: " << mol.natom << endl;
    cout << "Input Cartesian coordinates\n";
    mol.print_geom();
    // bond length
    cout << "Interatomic distances\n";
    // Bond length
    for (int i = 0; i < mol.natom; i++)
    {
        for (int j = 0; j < i; j++)
        {
            printf("%d %d %8.5f\n", i, j, mol.bond(i, j));
        }
    }
    // Bond angles
    cout << "\nBond angles\n";
    for (int i = 0; i < mol.natom; i++)
    {
        for (int j = 0; j < i; j++)
        {
            for (int k = 0; k < j; k++)
            {
                if (mol.bond(i, j) < 4.0 && mol.bond(j, k) < 4.0) // long distance ingnored
                {
                    printf("%d-%d-%d %10.6f\n", i, j, k, mol.angle(i, j, k) * (180 / acos(-1.0)));
                }
            }
        }
    }
    // Out of plane angles
    cout << "\nOut of plane angles\n";
    for (int i = 0; i < mol.natom; i++)
    {
        for (int j = 0; j < mol.natom; j++)
        {
            for (int k = 0; k < mol.natom; k++)
            {
                for (int l = 0; l < j; l++) // j k l in the plane ,k is the center (connected i j l),l<j avoid double counting
                {
                    if (i != j && i != k && i != l && j != k && j != l && k != l && mol.bond(i, k) < 4.0 && mol.bond(k, j) < 4.0 && mol.bond(k, l) < 4.0)
                    {
                        printf("%d-%d-%d-%d %10.6f\n", i, j, k, l, mol.opp(i, j, k, l) * (180 / acos(-1.0)));
                    }
                }
            }
        }
    }

    // Dihedral angle
    cout << "\nDiheral angle:\n";
    for (int i = 0; i < mol.natom; i++)
    {
        for (int j = 0; j < i; j++)
        {
            for (int k = 0; k < j; k++)
            {
                for (int l = 0; l < k; l++)
                {
                    if (mol.bond(i, j) < 4.0 && mol.bond(j, k) < 4.0 && mol.bond(k, l) < 4.0)
                    {
                        printf("%d-%d-%d-%d %10.6f\n", i, j, k, l, mol.dihedral(i, j, k, l) * (180 / acos(-1.0)));
                    }
                }
            }
        }
    }

    // Center of mass
    double M = 0.000000;
    double m_x = 0.000000;
    double m_y = 0.000000;
    double m_z = 0.000000;
    double mi;
    for (int i = 0; i < mol.natom; i++)
    {

        mi = masses[mol.zvals[i]];
        M += mi;
        m_x += mi * mol.geom[i][0];
        m_y += mi * mol.geom[i][1];
        m_z += mi * mol.geom[i][2];
    }
    double cm_x = m_x / M;
    double cm_y = m_y / M;
    double cm_z = m_z / M;

    printf("\nCenter of Mass: %20.12f %20.12f %20.12f\n", cm_x, cm_y, cm_z);
    // Translate molecular to center of mass
    mol.translate(-cm_x, -cm_y, -cm_z);

    // Principal Moments of Inertia
    mat A(3,3,arma::fill::zeros);
    for (int i = 0; i < mol.natom; i++)
    {
        mi = masses[mol.zvals[i]];
        A(0,0) += mi * (pow(mol.geom[i][1],2) + pow(mol.geom[i][2],2));
        A(1,1) += mi * (pow(mol.geom[i][0],2) + pow(mol.geom[i][2],2));
        A(2,2) += mi * (pow(mol.geom[i][0],2) + pow(mol.geom[i][1],2));
        A(0,1) -= mi * mol.geom[i][0] * mol.geom[i][1];
        A(0,2) -= mi * mol.geom[i][0] * mol.geom[i][2];
        A(1,2) -= mi * mol.geom[i][1] * mol.geom[i][2];


    }
    A(1,0) = A(0,1);
    A(2,0) = A(0,2);
    A(2,1) = A(1,2);
    
    A.print();
    vec eig_val;
    mat eig_vec;
    eig_sym(eig_val,eig_vec,A);

    cout<<"\nPrincipal Moments of Inertia:\n";
    for (int i = 0; i < eig_val.size(); i++)
    {
        cout<<eig_val(i)<<" ";
    }
    
    return 0;
}
