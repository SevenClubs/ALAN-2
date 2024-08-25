#include <iostream>
#include <cmath>

struct Matrix 
{
    unsigned rows;
    unsigned columns;
    float** values;
};

Matrix createMatrix(unsigned r, unsigned c, float* a);
void deleteMatrix(Matrix& m);
Matrix subMatrix(const Matrix& A, unsigned r, unsigned c);
Matrix genMatrix1();
Matrix genMatrix2();
Matrix genMatrix3();
Matrix genMatrix4();
void printMatrix (Matrix& m);
