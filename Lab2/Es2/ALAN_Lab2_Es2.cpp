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



// Returns true if float value is approximately 0
bool isZero(float x)
{
    return (abs(x) < pow(10, -6));
}



// Generic matrix product
Matrix matrixProduct(const Matrix& A, const Matrix& B, int& status)
{
    if (A.columns != B.rows)
    {
        status = -1; // data not valid
        return {0, 0, nullptr};
    }
    float res[A.rows * B.columns];
    for (unsigned i = 0; i < A.rows; ++i)
        for (unsigned j = 0; j < B.columns; ++j)
            res[i * B.columns + j] = 0;
    for (unsigned i = 0; i < A.rows; ++i)
    {
        for (unsigned j = 0; j < B.columns; ++j)
        {
            for (unsigned k = 0; k < A.columns; ++k)
            {
                res[i * B.columns + j] += A.values[i][k] * B.values[k][j];
            }
        }
    }
    return createMatrix(A.rows, B.columns, res);
}



// Knowing A and x, returns b as Ax = b
Matrix calcConstantMatrix(const Matrix& M) 
{
    float x[M.rows];
    for (unsigned i = 0; i < M.rows; ++i)
        x[i] = 1;
    Matrix x_vec = createMatrix(M.rows, 1, x);
    int status;
    Matrix res = matrixProduct(M, x_vec, status);
    deleteMatrix(x_vec);
    return res;
}



// Gauss method without partial pivoting
void GaussMethod(const Matrix& A, const Matrix& b, Matrix& x, int& status)
{
    if (b.columns != 1 || A.columns != b.rows) 
    {
        status = -1; // data not valid
        return;
    }

    float Ab_vec[(A.rows + 1) * A.columns];
    for (unsigned i = 0; i < A.rows; ++i)
        for (unsigned j = 0; j <= A.columns; ++j)
        {
            if (j != A.columns)
                Ab_vec[i * (A.columns + 1) + j] = A.values[i][j];
            else
                Ab_vec[i * (A.columns + 1) + j] = b.values[i][0];
        }

    Matrix Ab = createMatrix(A.rows, A.columns + 1, Ab_vec);

    // printMatrix(Ab);

    unsigned n_pivot = Ab.rows; 
    for (unsigned i = 0; i < Ab.rows && i < Ab.columns; ++i) 
    {
        // Check if the current pivot is zero
        bool notNullColumn = true;
        if (isZero(Ab.values[i][i]))
        {
            notNullColumn = false;
            for (unsigned j = i + 1; j < Ab.rows && !notNullColumn; ++j)
            {
                if(!isZero(Ab.values[j][i]))
                {
                    for (unsigned k = 0; k < Ab.columns; ++k)
                    {
                        float temp = Ab.values[i][k];
                        Ab.values[i][k] = Ab.values[j][k];
                        Ab.values[j][k] = temp;
                    }
                    notNullColumn = true;
                }
            }
        }
        if (notNullColumn)
        {
            // Gauss method
            for (unsigned j = i + 1; j < Ab.rows; ++j)
            {
                float m = -Ab.values[j][i] / Ab.values[i][i];
                for (unsigned k = i; k < Ab.columns; ++k)
                {
                Ab.values[j][k] += m * Ab.values[i][k];
                }
            }
        // printMatrix(Ab);
        }
        // If there is a column of zeroes, there is no pivot for that column
        else
        {
            --n_pivot;
        }
    }
    //cout << "\n" << "N_pivot = " << n_pivot << "\n";
    cout << "\n";
    // printMatrix(Ab);

    if (n_pivot > Ab.rows)
    {
        status = -2; // no existing solution
    }
    else if (n_pivot == Ab.rows)
    {
        status = 0; // exists one vector solution
        float x_vec[b.rows];
        for (unsigned i = n_pivot; i > 0; --i)
        {
            float temp = Ab.values[i - 1][Ab.columns - 1];
            for (unsigned j = i + 1; j < Ab.columns; ++j)
            {
                temp -= x_vec[j - 1] * Ab.values[i - 1][j - 1];
                // cout << "temp[" << i << "][" << j << "] = " << temp << "\n";
            }

            x_vec[i - 1] = temp / Ab.values[i - 1][i - 1];
        }
        x = createMatrix(b.rows, 1, x_vec);
    }
    else
    {
        if (isZero(Ab.values[n_pivot - 1][Ab.columns - 2]))
            status = -2; // no existing solution
        else
            status = Ab.rows - n_pivot; // exist infinity^(status) solutions
    }
    deleteMatrix(Ab);
}



// Gauss method with partial pivoting
void GaussMethodPart(const Matrix& A, const Matrix& b, Matrix& x, int& status)
{
    if (b.columns != 1 || A.columns != b.rows) 
    {
        status = -1; // data not valid
        return;
    }

    float Ab_vec[(A.rows + 1) * A.columns];
    for (unsigned i = 0; i < A.rows; ++i)
        for (unsigned j = 0; j <= A.columns; ++j)
        {
            if (j != A.columns)
                Ab_vec[i * (A.columns + 1) + j] = A.values[i][j];
            else
                Ab_vec[i * (A.columns + 1) + j] = b.values[i][0];
        }

    Matrix Ab = createMatrix(A.rows, A.columns + 1, Ab_vec);

    // printMatrix(Ab);

    unsigned n_pivot = Ab.rows; 
    for (unsigned i = 0; i < Ab.rows && i < Ab.columns; ++i) 
    {
        // Check if there is a column of zeroes
        bool notNullColumn = true;
        if (isZero(Ab.values[i][i]))
        {
            notNullColumn = false;
            for (unsigned j = i + 1; j < Ab.rows && !notNullColumn; ++j)
            {
                if(!isZero(Ab.values[j][i]))
                {
                    for (unsigned k = 0; k < Ab.columns; ++k)
                    {
                        float temp = Ab.values[i][k];
                        Ab.values[i][k] = Ab.values[j][k];
                        Ab.values[j][k] = temp;
                    }
                    notNullColumn = true;
                    // printMatrix(Ab);
                    // cout << "\n";
                }
            }
        }
        if (notNullColumn)
        {
            // Partial pivoting
            unsigned pivot = i;
            for (unsigned j = i + 1; j < Ab.rows; ++j) 
            {
                if (abs(Ab.values[j][i]) > Ab.values[pivot][i])
                    pivot = j;
            }
            if (pivot != i) 
            {
                for (unsigned j = i; j < Ab.columns; ++j) 
                {
                    float temp = Ab.values[i][j];
                    Ab.values[i][j] = Ab.values[pivot][j];
                    Ab.values[pivot][j] = temp;
                }
                // printMatrix(Ab);
                // cout << "\n";
            }

            // Gauss method
            for (unsigned j = i + 1; j < Ab.rows; ++j)
            {
                float m = -Ab.values[j][i] / Ab.values[i][i];
                for (unsigned k = i; k < Ab.columns; ++k)
                {
                Ab.values[j][k] += m * Ab.values[i][k];
                // printMatrix(Ab);
                // cout << "\n";
                }
            }
        // printMatrix(Ab);
        }
        else
        {
            --n_pivot;
        }
    }
    cout << "\n";
    // printMatrix(Ab);

    if (n_pivot > Ab.rows)
    {
        status = -2; // no existing solution
    }
    else if (n_pivot == Ab.rows)
    {
        status = 0; // exists one vector solution
        float x_vec[b.rows];
        for (unsigned i = n_pivot; i > 0; --i)
        {
            float temp = Ab.values[i - 1][Ab.columns - 1];
            for (unsigned j = i + 1; j < Ab.columns; ++j)
            {
                temp -= x_vec[j - 1] * Ab.values[i - 1][j - 1];
                // cout << "temp[" << i << "][" << j << "] = " << temp << "\n";
            }

            x_vec[i - 1] = temp / Ab.values[i - 1][i - 1];
        }
        x = createMatrix(b.rows, 1, x_vec);
    }
    else
    {
        if (isZero(Ab.values[n_pivot - 1][Ab.columns - 2]))
            status = -2; // no existing solution
        else
            status = Ab.rows - n_pivot; // exist infinity^(status) solutions
    }
    deleteMatrix(Ab);
}



// Knowing A and b, returns x as Ax = b (no partial pivoting)
Matrix calcVariableMatrix(const Matrix& A, const Matrix& b, int& status)
{
    Matrix x;
    GaussMethod(A, b, x, status);
    return x;
}



// Knowing A and b, returns x as Ax = b (with partial pivoting)
Matrix calcVariableMatrixPart(const Matrix& A, const Matrix& b, int& status)
{
    Matrix x;
    GaussMethodPart(A, b, x, status);
    return x;
}


int main() {
    int status;

    Matrix A1 = genMatrix1();
    // printMatrix(A1);
    Matrix b1 = calcConstantMatrix(A1);
    cout << "\nVector b1:\n";
    printMatrix(b1);
    Matrix x1 = calcVariableMatrix(A1, b1, status);
    if (status == 0)
    {
        cout << "\nSolution vector (no partial pivoting): \n";
        printMatrix(x1);
    }
    else
    {
        cout << "\nWrong result\n";
        return 1;
    }
    Matrix x1_bis = calcVariableMatrixPart(A1, b1, status);
    if (status == 0)
    {
        cout << "\nSolution vector (with partial pivoting): \n";
        printMatrix(x1_bis);
    }
    else
    {        
        cout << "\nWrong result\n";
        return 1;
    }
    cout << "\nError on solution without partial pivoting = " << infiniteNorm(x1) - 1
        << "\nError on solution with partial pivoting = " << infiniteNorm(x1_bis) - 1 << "\n";



    Matrix A2 = genMatrix2();
    // printMatrix(A2);
    Matrix b2 = calcConstantMatrix(A2);
    cout << "\nVector b2:\n";
    printMatrix(b2);
    Matrix x2 = calcVariableMatrix(A2, b2, status);
    if (status == 0)
    {
        cout << "\nSolution vector (no partial pivoting): \n";
        printMatrix(x2);
    }
    else
    {
        cout << "\nWrong result\n";
        return 1;
    }
    Matrix x2_bis = calcVariableMatrixPart(A2, b2, status);
    if (status == 0)
    {
        cout << "\nSolution vector (with partial pivoting): \n";
        printMatrix(x2_bis);
    }
    else
    {
        cout << "\nWrong result\n";
        return 1;
    }
    cout << "\nError on solution without partial pivoting = " << infiniteNorm(x2) - 1
        << "\nError on solution with partial pivoting = " << infiniteNorm(x2_bis) - 1 << "\n";



    Matrix A3 = genMatrix3();
    // printMatrix(A3);
    Matrix b3 = calcConstantMatrix(A3);
    cout << "\nVector b3:\n";
    printMatrix(b3);
    Matrix x3 = calcVariableMatrix(A3, b3, status);
    if (status == 0)
    {
        cout << "\nSolution vector (no partial pivoting): \n";
        printMatrix(x3);
    }
    else
    {
        cout << "\nWrong result\n";
        return 1;
    }
    Matrix x3_bis = calcVariableMatrixPart(A3, b3, status);
    if (status == 0)
    {
        cout << "\nSolution vector (with partial pivoting): \n";
        printMatrix(x3_bis);
    }
    else
    {
        cout << "\nWrong result\n";
        return 1;
    }
    cout << "\nError on solution without partial pivoting = " << infiniteNorm(x3) - 1
        << "\nError on solution with partial pivoting = " << infiniteNorm(x3_bis) - 1 << "\n";


    Matrix A4 = genMatrix4();
    // printMatrix(A4);
    Matrix b4 = calcConstantMatrix(A4);
    cout << "\nVector b4:\n";
    printMatrix(b4);
    Matrix x4 = calcVariableMatrix(A4, b4, status);
    if (status == 0)
    {
        cout << "\nSolution vector (no partial pivoting): \n";
        printMatrix(x4);
    }
    else
    {
        cout << "\nWrong result\n";
        return 1;
    }
    Matrix x4_bis = calcVariableMatrixPart(A4, b4, status);
    if (status == 0)
    {
        cout << "\nSolution vector (with partial pivoting): \n";
        printMatrix(x4_bis);
    }
    else
    {
        cout << "\nWrong result\n";
        return 1;
    }
    cout << "\nError on solution without partial pivoting = " << infiniteNorm(x4) - 1
        << "\nError on solution with partial pivoting = " << infiniteNorm(x4_bis) - 1 << "\n";

    cout << "\n";
    return 0;
}
