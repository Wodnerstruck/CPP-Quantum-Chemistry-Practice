#include "molecule.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <cassert>

Molecule::Molecule(const char *filename, int q)
{
  charge = q;
// open filename
  std::ifstream is(filename);
  if (!is) {
      printf("%s not found\n",filename);
      abort();
  }

  // read the number of atoms from filename
  is >> natom;

  // allocate space
  zvals = new int[natom];
  geom = new double* [natom];
  for(int i=0; i < natom; i++)
    geom[i] = new double[3];
  //store coordinates
  for(unsigned int i=0; i < natom; i++)
    is >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];

  is.close();
}

Molecule::~Molecule()
{
  delete[] zvals;
  for(int i=0; i < natom; i++)
    delete[] geom[i];
  delete[] geom;
}

void Molecule::print_geom ()
{
    for(int i=0;i<natom;i++){
        printf("%d %20.12f %20.12f %20.12f\n", zvals[i],geom[i][0],geom[i][1],geom[i][2]);
    }

}
double Molecule::bond(int a1,int a2)
{
     return (sqrt((geom[a1][0]-geom[a2][0])*(geom[a1][0]-geom[a2][0])
     +(geom[a1][1]-geom[a2][1])*(geom[a1][1]-geom[a2][1])
     +(geom[a1][2]-geom[a2][2])*(geom[a1][2] - geom[a2][2])));
}
double Molecule::unit(int a,int b,int coor)
{
  return -(geom[a][coor]- geom[b][coor])/bond(a,b);
}

double Molecule::angle(int a,int b, int c)
{
  return (acos(unit(b,a,0)*unit(b,c,0)+unit(b,a,1)*unit(b,c,1)
  +unit(b,a,2)*unit(b,c,2)))*(180/acos(-1.0));

}


