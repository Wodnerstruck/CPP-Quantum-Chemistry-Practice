#include "SCF.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cmath>

using namespace std;

int main(int argc, char** argv)
{
    integrals Int("./data/h2o/DZP");
    SCF H2O(Int);
    H2O.iterat(Int,"MP2");
    return 0;

}
