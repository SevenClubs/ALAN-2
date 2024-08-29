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
    cout << "\nMatrice A1:\n"; 
    printMatrix(A1);
    
    Matrix A2 = genMatrix2();
    cout << "\nMatrice A2:\n";
    printMatrix(A2);
    
    Matrix A3 = genMatrix3();
    cout << "\nMatrice A3:\n";
    printMatrix(A3);
    
    Matrix A4 = genMatrix4();
    cout << "\nMatrice A4:\n";
    printMatrix(A4);

    cout << "\nNorma infinito matrice A1 = " << infiniteNorm(A1) << "\n";
    cout << "\nNorma infinito matrice A2 = " << infiniteNorm(A2) << "\n"; 
    cout << "\nNorma infinito matrice A3 = " << infiniteNorm(A3) << "\n"; 
    cout << "\nNorma infinito matrice A4 = " << infiniteNorm(A4) << "\n\n";

    return 0;
}
