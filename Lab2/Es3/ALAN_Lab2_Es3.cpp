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
    return (abs(x) < pow(10, -10));
}



// Returns the determinant of a matrix
float determinant(const Matrix& A, int& status)
{
    if (A.columns != A.rows)
    {
        status = -1; // non-square matrix given
        return 0;
    }
    if (A.columns == 0)
    {
        status = 1; // degenerate case
        return 0;
    }
    if (A.columns == 1) 
    {
        status = 0;
        return A.values[0][0];
    }
    float res = 0;
    for (unsigned i = 1; i <= A.rows; ++i)
    {
        if (!isZero(A.values[i - 1][0]))
        {
            Matrix temp = subMatrix(A, i, 1);
            // cout << "\nCurrent submatrix rows = " << temp.rows << "\nColumns = " << temp.columns << "\nValues:\n";
            // printMatrix(temp);
            res += A.values[i - 1][0] * determinant(temp, status);
            deleteMatrix(temp);
        }
    }
    // cout << "\nDeterminant = " << res << "\n";
    status = 0;
    return res;
}



// Returns inverse matrix of the given matrix
Matrix inverseMatrix(const Matrix& A, int& status)
{
    if (A.columns != A.rows)
    {
        status = -1; // non-square matrix given
        return {0, 0, nullptr};
    }
    if (A.columns == 0)
    {
        status = 1; // degenerate case
        return {0, 0, nullptr};
    }
    if (isZero(determinant(A, status)))
    {
        status = -2; // non-invertible matrix
        return {0, 0, nullptr};
    }
    float val[A.rows * A.columns];
    for (unsigned i = 1; i <= A.rows; ++i)
    {
        for (unsigned j = 1; j <= A.columns; ++j)
        {
            Matrix temp = subMatrix(A, i, j);
            val[(j - 1) * A.rows + (i - 1)] = determinant(temp, status);
            deleteMatrix(temp);
            if (status != 0) // could not compute determinant
                return {0, 0, nullptr};
        }
    }
    float det_A = determinant(A, status);
    if (status != 0) // could not compute determinant
        return {0, 0, nullptr};

    for (unsigned i = 0; i < A.rows * A.columns; ++i)
        val[i] /= det_A;
    status = 0;
    return createMatrix(A.rows, A.columns, val);
}



// Knowing A and x, returns b as Ax = b
Matrix calcConstantMatrix(const Matrix& M) 
{
    float res[M.rows];
    for (unsigned i = 0; i < M.rows; ++i)
        res[i] = 0;
    for (unsigned i = 1; i <= M.rows; ++i)
    {
        for (unsigned j = 1; j <= M.columns; ++j)
        {
            res[i - 1] += M.values[i - 1][j - 1];
        }
    }
    return createMatrix(M.rows, 1, res);
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
        // Check if there is a column of zeroes
        bool notNullColumn = true;
        // Check if the current pivot is zero
        if (isZero(Ab.values[i][i]))
        {
            // If so, look for another pivot in the same column
            for (unsigned j = i + 1; j < Ab.rows && notNullColumn; ++j)
            {
                if(!isZero(Ab.values[j][i]))
                {
                    notNullColumn = false;
                    for (unsigned k = 0; k < Ab.columns; ++k)
                    {
                        float temp = Ab.values[i][k];
                        Ab.values[i][k] = Ab.values[j][k];
                        Ab.values[j][k] = temp;
                    }
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
        return;
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



// Computes pert_b = b + delta_b where delta_b = infiniteNorm(b) * (0.01, -0.01, ...)
Matrix calcPerturbation(const Matrix& b)
{
    float norm_b = infiniteNorm(b);
    float pert_b_val[b.rows];
    for (unsigned i = 0; i < b.rows; ++i) 
    {
        if (i % 2 == 0)
            pert_b_val[i] = b.values[i][0] - norm_b * 0.01;
        else
            pert_b_val[i] = b.values[i][0] + norm_b * 0.01;
    }
    return createMatrix(b.rows, 1, pert_b_val);
}

int main() {
    int status;
    
    Matrix A1 = genMatrix1();
    Matrix A1_inv = inverseMatrix(A1, status);
    if (status == 0)
    {
        cout << "\nInverse matrix of A1:\n";
        printMatrix(A1_inv);
        cout << "\ndet(A1) = " << determinant(A1, status) << ", det(A1_inv) = " << determinant(A1_inv, status) << "\n";
        cout << "det(A1) * det(A1_inv) = " << determinant(A1, status) * determinant(A1_inv, status) << "\n";
        cout << "\nMaximum condition number of A1 = " << infiniteNorm(A1) * infiniteNorm(A1_inv) << "\n";
    }
    else
    {
        cout << "\nCould not compute inverse of A1\nStatus = " << status << "\n";
        return 1;
    }
    Matrix b1 = calcConstantMatrix(A1);
    cout << "\nVector b1:\n";
    printMatrix(b1);
    cout << "\n";
    Matrix pert_b1 = calcPerturbation(b1);
    cout << "\nVector pert_b1:\n";
    printMatrix(pert_b1);
    cout << "\n";
    float in_pert = abs(infiniteNorm(b1) - infiniteNorm(pert_b1)) / infiniteNorm(b1);
    Matrix x1 = calcVariableMatrix(A1, b1, status);
    if (status == 0)
    {
        cout << "\nVector x1: \n";
        printMatrix(x1);
    }
    else
    {
        cout << "\nWrong result\n";
        return 1;
    }
    Matrix pert_x1 = calcVariableMatrix(A1, pert_b1, status);
    if (status == 0)
    {
        cout << "\nVector pert_x1: \n";
        printMatrix(pert_x1);
    }
    else
    {
        cout << "\nWrong result\n";
        return 1;
    }
    float out_pert = abs(infiniteNorm(pert_x1) - infiniteNorm(x1)) / infiniteNorm(x1);
    cout << "\nInput perturbance = " << in_pert << "\nOutput perturbance = " << out_pert;
    cout << "\nEmpirical condition number = " << out_pert / in_pert << "\n";



    Matrix A2 = genMatrix2();
    Matrix A2_inv = inverseMatrix(A2, status);
    if (status == 0)
    {
        cout << "\nInverse matrix of A2:\n";
        printMatrix(A2_inv);
        cout << "\ndet(A2) = " << determinant(A2, status) << ", det(A2_inv) = " << determinant(A2_inv, status) << "\n";
        cout << "det(A2) * det(A2_inv) = " << determinant(A2, status) * determinant(A2_inv, status) << "\n";
        cout << "\nMaximum condition number of A2 = " << infiniteNorm(A2) * infiniteNorm(A2_inv) << "\n";
    }
    else
    {
        cout << "\nCould not compute inverse of A2\nStatus = " << status << "\n";
        return 1;
    }
    Matrix b2 = calcConstantMatrix(A2);
    cout << "\nVector b2:\n";
    printMatrix(b2);
    cout << "\n";
    Matrix pert_b2 = calcPerturbation(b2);
    cout << "\nVector pert_b2:\n";
    printMatrix(pert_b2);
    cout << "\n";
    in_pert = abs(infiniteNorm(b2) - infiniteNorm(pert_b2)) / infiniteNorm(b2);
    Matrix x2 = calcVariableMatrix(A2, b2, status);
    if (status == 0)
    {
        cout << "\nVector x2: \n";
        printMatrix(x2);
    }
    else
    {
        cout << "\nWrong result\n";
        return 1;
    }
    Matrix pert_x2 = calcVariableMatrix(A2, pert_b2, status);
    if (status == 0)
    {
        cout << "\nVector pert_x2: \n";
        printMatrix(pert_x2);
    }
    else
    {
        cout << "\nWrong result\n";
        return 1;
    }
    out_pert = abs(infiniteNorm(pert_x2) - infiniteNorm(x2)) / infiniteNorm(x2);
    cout << "\nInput perturbance = " << in_pert << "\nOutput perturbance = " << out_pert;
    cout << "\nEmpirical condition number = " << out_pert / in_pert << "\n";
    


    Matrix A3 = genMatrix3();
    Matrix A3_inv = inverseMatrix(A3, status);
    if (status == 0)
    {
        cout << "\nInverse matrix of A3:\n";
        printMatrix(A3_inv);
        cout << "\ndet(A3) = " << determinant(A3, status) << ", det(A3_inv) = " << determinant(A3_inv, status) << "\n";
        cout << "det(A3) * det(A3_inv) = " << determinant(A3, status) * determinant(A3_inv, status) << "\n";
        cout << "\nMaximum condition number of A3 = " << infiniteNorm(A3) * infiniteNorm(A3_inv) << "\n";
    }
    else
    {
        cout << "\nCould not compute inverse of A3\nStatus = " << status << "\n";
        return 1;
    }
    Matrix b3 = calcConstantMatrix(A3);
    cout << "\nVector b3:\n";
    printMatrix(b3);
    cout << "\n";
    Matrix pert_b3 = calcPerturbation(b3);
    cout << "\nVector pert_b3:\n";
    printMatrix(pert_b3);
    cout << "\n";
    in_pert = abs(infiniteNorm(b3) - infiniteNorm(pert_b3)) / infiniteNorm(b3);
    Matrix x3 = calcVariableMatrix(A3, b3, status);
    if (status == 0)
    {
        cout << "\nVector x3: \n";
        printMatrix(x3);
    }
    else
    {
        cout << "\nWrong result\n";
        return 1;
    }
    Matrix pert_x3 = calcVariableMatrix(A3, pert_b3, status);
    if (status == 0)
    {
        cout << "\nVector pert_x3: \n";
        printMatrix(pert_x3);
    }
    else
    {
        cout << "\nWrong result\n";
        return 1;
    }
    out_pert = abs(infiniteNorm(pert_x3) - infiniteNorm(x3)) / infiniteNorm(x3);
    cout << "\nInput perturbance = " << in_pert << "\nOutput perturbance = " << out_pert;
    cout << "\nEmpirical condition number = " << out_pert / in_pert << "\n";
    


    Matrix A4 = genMatrix4();
    Matrix A4_inv = inverseMatrix(A4, status);
    if (status == 0)
    {
        cout << "\nInverse matrix of A4:\n";
        printMatrix(A4_inv);
        cout << "\ndet(A4) = " << determinant(A4, status) << ", det(A4_inv) = " << determinant(A4_inv, status) << "\n";
        cout << "det(A4) * det(A4_inv) = " << determinant(A4, status) * determinant(A4_inv, status) << "\n";
        cout << "\nMaximum condition number of A4 = " << infiniteNorm(A4) * infiniteNorm(A4_inv) << "\n";
    }
    else
    {
        cout << "\nCould not compute inverse of A4\nStatus = " << status << "\n";
        return 1;
    }
    Matrix b4 = calcConstantMatrix(A4);
    cout << "\nVector b4:\n";
    printMatrix(b4);
    cout << "\n";
    Matrix pert_b4 = calcPerturbation(b4);
    cout << "\nVector pert_b4:\n";
    printMatrix(pert_b4);
    cout << "\n";
    in_pert = abs(infiniteNorm(b4) - infiniteNorm(pert_b4)) / infiniteNorm(b4);
    Matrix x4 = calcVariableMatrix(A4, b4, status);
    if (status == 0)
    {
        cout << "\nVector x4: \n";
        printMatrix(x4);
    }
    else
    {
        cout << "\nWrong result\n";
        return 1;
    }
    Matrix pert_x4 = calcVariableMatrix(A4, pert_b4, status);
    if (status == 0)
    {
        cout << "\nVector pert_x4: \n";
        printMatrix(pert_x4);
    }
    else
    {
        cout << "\nWrong result\n";
        return 1;
    }
    out_pert = abs(infiniteNorm(pert_x4) - infiniteNorm(x4)) / infiniteNorm(x4);
    cout << "\nInput perturbance = " << in_pert << "\nOutput perturbance = " << out_pert;
    cout << "\nEmpirical condition number = " << out_pert / in_pert << "\n";
    

    return 0;
}
