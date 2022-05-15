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
    //bond length
    cout<<"Interatomic distances\n";
    for(int i=0;i<mol.natom;i++){
        for(int j=0;j<i;j++){
            printf("%d %d %8.5f\n",i , j, mol.bond(i,j));
        }
    }
    //bond angles
    cout<<"\nBond angles\n";
    for (int i = 0; i < mol.natom; i++)
    {
        for (int j = 0; j < i; j++)
        {
            for (int k = 0; k < j; k++)
            {
                if(mol.bond(i,j) < 4.0 && mol.bond(j,k) < 4.0)//long distance ingnored
                {
                    printf("%d-%d-%d %10.6f\n",i,j,k,mol.angle(i,j,k) * (180 / acos(-1.0)));
                }
            }
            
        }
        
    }
    //out of plane angles
    cout<<"\nOut of plane angles\n";
    for (int i = 0; i < mol.natom; i++)
    {
        for(int j = 0;j < mol.natom;j++ )
        {
            for (int k = 0; k < mol.natom; k++)
            {
                for (int l = 0; l < j; l++)
                {
                    if(i!=j&&i!=k&&i!=l&&j!=k&&j!=l&&k!=l&&mol.bond(i,k)<4.0&&mol.bond(k,j)<4.0&&mol.bond(k,l)<4.0)
                    {
                        printf("%d-%d-%d-%d %10.6f\n",i,j,k,l,mol.opp(i,j,k,l) * (180/acos(-1.0)));
                    }
                }
                
            }
            
        }
    }
    
    
    


    
return 0;

}
