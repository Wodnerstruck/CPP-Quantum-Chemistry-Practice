#ifndef MOLECULE_H
#define MOLECULE_H


#include<string>
using namespace std;
 
class Molecule
{
  public:
    int natom;
    int charge;
    int *zvals;
    double **geom;
   // string point_group;
 
    void print_geom();
    //void rotate(double phi);
   // void translate(double x, double y, double z);
    double bond(int a1, int a2);
    double unit(int a,int b,int coor);//coor : 0 -x 1 -y 2 -z
    double angle(int a1, int a2, int a3);
    //double torsion(int a1, int a2, int a3, int a4);
 
    Molecule(const char *filename,int q );
    ~Molecule();
};

#endif