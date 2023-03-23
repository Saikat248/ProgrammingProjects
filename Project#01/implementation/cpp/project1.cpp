#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;


int main(int argc, char *argv[])
{
    ifstream input("geom.dat");

    int natom;

    input >> natom;

    int *zval = new int[natom];
    double *x = new double[natom];
    double *y = new double[natom];
    double *z = new double[natom];

    for (int i = 0; i < natom; i++){
        input >> zval[i] >> x[i] >> y[i] >> z[i];
    }
    input.close();

    cout << "No of atoms: " << natom << endl;
    cout << "Input Cartesian coordinates:\n";
    for(int i=0; i < natom; i++)
        printf("%d %20.12f %20.12f %20.12f\n", (int) zval[i], x[i], y[i], z[i]);

   delete[] zval;
   delete[] x;  delete[] y;  delete[] z;

    return 0;
}