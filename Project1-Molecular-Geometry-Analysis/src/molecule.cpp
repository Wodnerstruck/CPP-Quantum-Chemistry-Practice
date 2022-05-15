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
     return (sqrt((geom[a1][0]-geom[a2][0]) * (geom[a1][0]-geom[a2][0])
     +(geom[a1][1]-geom[a2][1]) * (geom[a1][1]-geom[a2][1])
     +(geom[a1][2]-geom[a2][2]) * (geom[a1][2] - geom[a2][2])));
}
double Molecule::unit(int a,int b,int coor)
{
  return -(geom[a][coor]- geom[b][coor])/bond(a,b);
}

double Molecule::angle(int a,int b, int c)
{
  return (acos(unit(b,a,0) * unit(b,c,0)+unit(b,a,1) * unit(b,c,1)
  +unit(b,a,2) * unit(b,c,2))) ;
  //Plz trans radian to degree of angle in the main function,don't trans in this 
  //function,otherwise you would get wrong values when function opp call function angle
}

double Molecule::opp(int a,int b,int c,int d)
{
  //cross product of cb(kj) & cd(kl)
  double e_bcd_x = unit(c,b,1) * unit(c,d,2)-unit(c,b,2)*unit(c,d,1);
 
  double e_bcd_y = unit(c,b,2) * unit(c,d,0)-unit(c,b,0)*unit(c,d,2);
 
  double e_bcd_z = unit(c,b,0) * unit(c,d,1)-unit(c,b,1)*unit(c,d,0);
 
//dot product with ca(ki) ,then divided by sin(theta[bcd(jkl)])

  double e_xx = e_bcd_x * unit(c,a,0);
  double e_yy = e_bcd_y * unit(c,a,1);
  double e_zz = e_bcd_z * unit(c,a,2);
  //calc theta
  double theta = (e_xx+e_yy+e_zz) / sin(angle(b,c,d));
  // the value may out of (-1,1)
  if (theta<-1.0) theta = asin(-1.0);
    else if(theta>1.0) theta = asin(1.0);
      else theta = asin(theta);

        return theta;

}


