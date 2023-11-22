#include <iostream>
#include <cmath>

using namespace std;

// Funzione che calcola la norma infinito di una matrice

float veryAbs(float a)
{
    if (a>=0 ) return a;
    return a - 2*a;
}

struct Matrix 
{
    unsigned rows;
    unsigned columns;
    float** values;
};

Matrix createMatrix(unsigned r, unsigned c, float* a)
{
    Matrix m;
    m.rows = r;
    m.columns = c;
    m.values = new float* [m.rows];
    for (unsigned i = 0; i < r; ++i) 
        m.values[i] = new float [m.columns];
    for (unsigned i = 0; i < m.rows; ++i) 
        for (unsigned j = 0; j < m.columns; ++j)
            m.values[i][j] = a[i*m.columns + j];
    return m;
}

void printMatrix (Matrix m)
{
    for (unsigned i = 0; i < m.rows; ++i) 
    {
        for (unsigned j = 0; j < m.columns; ++j)
            cout << m.values[i][j] << "  ";
        cout << "\n";
    }
}

float infiniteNorm (Matrix m)
{
    float max = 0;
    for (unsigned i = 0; i < m.rows; ++i) 
    {
        float sum = 0;
        for (unsigned j = 0; j < m.columns; ++j)
            sum += veryAbs(m.values[i][j]);
        if (sum > max)
            max = sum;
    }
    return max;
}

float binCoeff(int a, int b)
{
    int x, y;
    if (a == b || b == 0)
        return 1;
    if (a / 2 >= b)
    {
        x = a - b + 1;
        y = b;
    }
    else
    {
        x = b + 1;
        y = a - b;
    }
    float res = x++;
    while (x < a || y > 1)
    {
        if (x <= a)
            res *= x++;
        if (y > 1)
            res /= y--;
    }
    return res;
}

int main() {
    const int nRows = 4;
    const int nCols = 4;
    float matr[nRows * nCols] = {3, 1, -1, 0, 0, 7, -3, 0, 0, -3, 9, -2, 0, 0, 4, -10};
    Matrix m = createMatrix(nRows, nCols, matr);
    printMatrix(m);
    cout << "\n" << infiniteNorm(m) << "\n\n";

    float matr2[nRows * nCols] = {2, 4, -2, 0, 1, 3, 0, 1, 3, -1, 1, 2, 0, -1, 2, 1};
    Matrix m2 = createMatrix(nRows, nCols, matr2);
    printMatrix(m2);
    cout << "\n" << infiniteNorm(m2) << "\n\n";

    const int M = 7;
    const int N = 7;
    float matr3[M * N];
    for (unsigned i = 1; i <= M; ++i)
    {
        for (unsigned j = 1; j <= N; ++j)
        {
            matr3[(i - 1) * M + (j - 1)] = binCoeff(i + j - 2, i - 1);
        }
    }

    Matrix m3 = createMatrix(M, N, matr3);
    printMatrix(m3);
    cout << "\n" << infiniteNorm(m3) << "\n\n";

    return 0;
}
