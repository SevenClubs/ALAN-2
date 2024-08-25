#include <iostream>
#include <cmath>
#include "ALAN_Lab2.h"

using namespace std;


// Computes infinity matrix norm
float infiniteNorm (Matrix m)
{
    float max = 0;
    for (unsigned i = 0; i < m.rows; ++i) 
    {
        float sum = 0;
        for (unsigned j = 0; j < m.columns; ++j)
            sum += abs(m.values[i][j]);
        if (sum > max)
            max = sum;
    }
    return max;
}



int main() {
    Matrix A1 = genMatrix1();
    cout << "Matrice A_1:\n"; 
    printMatrix(A1);
    cout << "\nNorma infinito matrice A_1 = " << infiniteNorm(A1) << "\n\n";

    Matrix A2 = genMatrix2();
    cout << "Matrice A_2:\n";
    printMatrix(A2);
    cout << "\nNorma infinito matrice A_2 = " << infiniteNorm(A2) << "\n\n";

    Matrix A3 = genMatrix3();
    cout << "Matrice A_3:\n";
    printMatrix(A3);
    cout << "\nNorma infinito matrice A_3 = " << infiniteNorm(A3) << "\n\n";

    Matrix A4 = genMatrix4();
    cout << "Matrice A_4:\n";
    printMatrix(A4);
    cout << "\nNorma infinito matrice A_4 = " << infiniteNorm(A4) << "\n\n";

    return 0;
}
