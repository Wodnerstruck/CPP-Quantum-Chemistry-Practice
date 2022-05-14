#include"molecule.h"
#include<iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>

using namespace std;


//function end
int main()
{
    Molecule mol("./data/acetaldehyde.dat",0);
    cout<<"Number of Atom: "<<mol.natom<<endl;
    cout<<"Input Cartesian coordinates\n";
    mol.print_geom();
    cout<<"Interatomic distances\n";
    for(int i=0;i<mol.natom;i++){
        for(int j=0;j<i;j++){
            printf("%d %d %8.5f\n",i , j, mol.bond(i,j));
        }
    }

    cout<<"\nBond angles\n";
    for (int i = 0; i < mol.natom; i++)
    {
        for (int j = 0; j < i; j++)
        {
            for (int k = 0; k < j; k++)
            {
                if(mol.bond(i,j) < 4.0 && mol.bond(j,k) < 4.0)//long distance ingnored
                {
                    printf("%d-%d-%d %10.6f\n",i,j,k,mol.angle(i,j,k));
                }
            }
            
        }
        
    }
    


}
