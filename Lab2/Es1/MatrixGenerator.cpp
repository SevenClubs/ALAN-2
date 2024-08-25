#include "ALAN_Lab2.h"

using namespace std;


// Array a contains the matrix values
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



// Delete matrix for avoiding memory leaks
void deleteMatrix(Matrix& m)
{
    for (unsigned i = 0; i < m.rows; ++i)
        delete [] m.values[i];
    delete [] m.values;
    m.values = nullptr;
    m.rows = 0;
    m.columns = 0;
}



// Return a copy of matrix A without r-th row and c-th column
// r in 1, ..., A.rows, c in 1, ..., A.columns
Matrix subMatrix(const Matrix& A, unsigned r, unsigned c)
{
    if (r < 1 || r > A.rows || c < 1 || c > A.columns) 
    {
        string ERROR = "Invalid arguments in subMatrix";
        throw ERROR;
    }
    Matrix m;
    m.rows = A.rows - 1;
    m.columns = A.columns - 1;
    m.values = new float* [m.rows];
    for (unsigned i = 0; i < m.rows; ++i) 
        m.values[i] = new float [m.columns];
    for (unsigned i = 0; i < A.rows; ++i)
    {
        if (i < r - 1)
        {
            for (unsigned j = 0; j < A.columns; ++j)
                if (j < c - 1)
                    m.values[i][j] = A.values[i][j];
                else if (j >= c)
                    m.values[i][j - 1] = A.values[i][j];
        }
        else if (i >= r)
        {
            for (unsigned j = 0; j < A.columns; ++j)
                if (j < c - 1)
                    m.values[i - 1][j] = A.values[i][j];
                else if (j >= c)
                    m.values[i - 1][j - 1] = A.values[i][j];
        }
    }
    return m;
}



void printMatrix (Matrix& m)
{
    for (unsigned i = 0; i < m.rows; ++i) 
    {
        for (unsigned j = 0; j < m.columns; ++j)
            cout << m.values[i][j] << "  ";
        cout << "\n";
    }
}



// Compute binomial coefficient 
float binCoeff(int a, int b)
{
    int x, y;
    if (a == b || b == 0)
        return 1;
    
    // a! can be a very big number, so it is used a different algorithm to avoid that numerator and denominator are too different in size, avoiding excessive approximation
    if (a / 2 >= b)
    {
        // a! / (b!*(a - b)!) == a*(a - 1)*...*(a - b + 1) / b!
        x = a - b + 1;
        y = b;
    }
    else
    {
        // a! / (b! (a - b)!) == 1*2*...*(b) / (a - b)!
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



// Create matrix A1
Matrix genMatrix1() 
{
    const unsigned n = 4;
    float matr[n * n] = {3, 1, -1, 0, 0, 7, -3, 0, 0, -3, 9, -2, 0, 0, 4, -10};
    return createMatrix(n, n, matr);
}



// Create matrix A2
Matrix genMatrix2()
{
    const unsigned n = 4;
    float matr2[n * n] = {2, 4, -2, 0, 1, 3, 0, 1, 3, -1, 1, 2, 0, -1, 2, 1};
    return createMatrix(n, n, matr2);
}



// Create matrix A3
Matrix genMatrix3()
{
    const unsigned n = 10;
    float matr3[n * n];
    for (unsigned i = 1; i <= n; ++i)
    {
        for (unsigned j = 1; j <= n; ++j)
        {
            matr3[(i - 1) * n + (j - 1)] = binCoeff(i + j - 2, i - 1);
        }
    }
    return createMatrix(n, n, matr3);
}



// Create matrix A4
Matrix genMatrix4()
{
    const unsigned d0 = 2;
    const unsigned d1 = 2;
    const unsigned n = 10 * (d1 + 1) + d0;
    float matr4[n * n];
    for (unsigned i = 1; i <= n; ++i)
    {
        for (unsigned j = 1; j <= n; ++j)
        {
            if (i == j)
                matr4[(i - 1) * n + (j - 1)] = 2;
            else if (i == j + 1 || j == i + 1)
                matr4[(i - 1) * n + (j - 1)] = 1;
            else
                matr4[(i - 1) * n + (j - 1)] = 0;
        }
    }
    return createMatrix(n, n, matr4);
}
