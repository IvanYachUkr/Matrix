#include <iostream>
#include <vector>
#include <ostream>
#include <iomanip>
#include <cmath>


class Matrix {
public:
    std::vector<std::vector<double> > matrix;
    int rows_number;
    int columns_number;

    Matrix(int rows , int cols) {
        rows_number = rows;
        columns_number = cols;
        for (int i = 0; i < rows; i++) {
            std::vector<double> column;

            for (int j = 0; j < cols; j++) {
                column.push_back(0.0);
            }
            matrix.push_back(column);
        }
    }

    Matrix operator+(const Matrix &matrix1) {
        Matrix M(rows_number, columns_number);

        for (int i = 0; i < rows_number; i++) {
            for (int j = 0; j < columns_number; j++) {
                M.matrix[i][j] = matrix1.matrix[i][j] + matrix[i][j];
            }
        }
        return M;
    };

    Matrix operator+(const double num) {
        Matrix M(rows_number, columns_number);

        for (int i = 0; i < rows_number; i++) {
            for (int j = 0; j < columns_number; j++) {
                M.matrix[i][j] = num * (i == j) + matrix[i][j];
            }
        }
        return M;
    };

    Matrix operator-(const Matrix &matrix1) {
        Matrix M(rows_number, columns_number);

        for (int i = 0; i < rows_number; i++) {
            for (int j = 0; j < columns_number; j++) {
                M.matrix[i][j] = matrix[i][j] - matrix1.matrix[i][j];
            }
        }
        return M;
    };

    Matrix operator-(const double num) {
        Matrix M(rows_number, columns_number);

        for (int i = 0; i < rows_number; i++) {
            for (int j = 0; j < columns_number; j++) {
                M.matrix[i][j] = matrix[i][j] - num * (i == j);
            }
        }
        return M;
    };

    Matrix operator*(const Matrix &matrix1) {
        Matrix M(rows_number, matrix1.columns_number);
        for (int i = 0; i < rows_number; i++) {
            for (int j = 0; j < matrix1.columns_number; j++) {

                for (int m = 0; m < columns_number; m++) {
                    M.matrix[i][j] += matrix[i][m] * matrix1.matrix[m][j];

                }

            }
        }
        return M;
    };

    Matrix operator*(const double num) {
        Matrix M(rows_number, columns_number);

        for (int i = 0; i < rows_number; i++) {
            for (int j = 0; j < columns_number; j++) {
                M.matrix[i][j] = num * matrix[i][j];
            }
        }
        return M;
    };

    Matrix operator+() {
        //transpose operation
        Matrix M(columns_number, rows_number);
        for (int i = 0; i < columns_number; i++) {
            for (int j = 0; j < rows_number; j++) {
                M.matrix[i][j] = matrix[j][i];
            }
        }
        return M;
    }

    friend std::ostream &operator<<(std::ostream &out, const Matrix &m);

    static double dot_product(Matrix vector1, const Matrix& vector2) {

        return ((+vector1) * vector2).matrix[0][0];
    }

    static double norm(Matrix vector1, Matrix vector2) {
        double res = 0;
        for (int i = 0; i < vector1.rows_number; ++i) {
            double curr_res = vector1.matrix[i][0] - vector2.matrix[i][0];


            if (fabs(curr_res) > res) {
                res = fabs(curr_res);
            }
        }

        return res;

    }



    std::vector<double>& operator[](int row)
    {

        return matrix[row];
    }

    static Matrix Strassen_2x2(Matrix A, Matrix B){
        /*
         Strassen method for multiplying matrices 2x2
         */
        double p1 = (A[0][0] + A[1][1])*(B[0][0]+B[1][1]);
        double p2 = (A[1][0] + A[1][1])*(B[0][0]);
        double p3 = (A[0][0])*(B[0][1]-B[1][1]);
        double p4 = (A[1][1])*(B[1][0]-B[0][0]);
        double p5 = (A[0][0] + A[0][1])*(B[1][1]);
        double p6 = (A[1][0] - A[0][0])*(B[0][0]+B[0][1]);
        double p7 = (A[0][1] - A[1][1])*(B[1][0]+B[1][1]);

        double c_1_1 = p1+p4-p5+p7;
        double c_1_2 = p3+p5;
        double c_2_1 = p2+p4;
        double c_2_2 = p1+p3-p2+p6;

        Matrix C(2,2);
        C[0][0] = c_1_1;
        C[0][1] = c_1_2;
        C[1][0] = c_2_1;
        C[1][1] = c_2_2;
        return C;


    }
    static Matrix fill_submatrix(Matrix A, int rows, int row_num_of_square = 0, int column_num_of_square = 0){
        /*
        method to fill submatrices
         by default copies whole matrix
        A - full matrix
        to copy square submatrix we use row_num_of_square and column_num_of_square to define it
            (A 1_1    A 1_2)
        A =
            (A 2_1    A 2_2)
         */
        Matrix B(rows, rows);

        int begin_row[] = {0, 0, A.rows_number / 2};
        int end_row[] = {A.rows_number, A.rows_number / 2, A.rows_number};

        int begin_column[] = {0, 0, A.columns_number/2};
        int end_column[] = {A.columns_number,A.columns_number/2, A.columns_number};

        int first_index_row = begin_row[row_num_of_square];
        int last_idx_row = end_row[row_num_of_square];

        int first_index_column = begin_column[column_num_of_square];
        int last_idx_column = end_column[column_num_of_square];

        for (int i = first_index_row; i < last_idx_row; ++i) {
            for (int j = first_index_column; j < last_idx_column; ++j) {
                B[i - first_index_row][j - first_index_column]= A[i][j];
            }
        }
        return B;
    }

    static Matrix Strassen(Matrix A, Matrix B){
        /*
         Strassen method for multiplying matrices A[m x n] * B [n x k]
         */

        //if both matrices are 2x2 apply strassen method for 2x2 matrices
        if (A.rows_number == A.columns_number
        and
        B.rows_number == B.columns_number
        and
        A.rows_number == 2
        and
        B.rows_number == 2){
            return Strassen_2x2(A, B);
        }


        Matrix C(A.rows_number, B.columns_number);


        //make both matrices square (2^n x 2^n)

        int max_rows_columns = std::max(std::max(A.rows_number,A.columns_number), B.columns_number);
        double a = log2(max_rows_columns);
        int rows_num_w_padding = pow(2, floor(a) + ((a - floor(a)) > 0) );
        int rows_num_submatrix = rows_num_w_padding / 2;

        Matrix A_w_padding = fill_submatrix(A, rows_num_w_padding);
        Matrix B_w_padding = fill_submatrix(B, rows_num_w_padding);


        Matrix A_1_1 = fill_submatrix(A_w_padding, rows_num_submatrix, 1, 1);
        Matrix A_1_2 = fill_submatrix(A_w_padding, rows_num_submatrix, 1, 2);
        Matrix A_2_1 = fill_submatrix(A_w_padding, rows_num_submatrix, 2, 1);
        Matrix A_2_2 = fill_submatrix(A_w_padding, rows_num_submatrix, 2, 2);



        Matrix B_1_1 = fill_submatrix(B_w_padding, rows_num_submatrix, 1, 1);
        Matrix B_1_2 = fill_submatrix(B_w_padding, rows_num_submatrix, 1, 2);
        Matrix B_2_1 = fill_submatrix(B_w_padding, rows_num_submatrix, 2, 1);
        Matrix B_2_2 = fill_submatrix(B_w_padding, rows_num_submatrix, 2, 2);

        Matrix p1 = Strassen((A_1_1 + A_2_2),(B_1_1+B_2_2));
        Matrix p2 = Strassen((A_2_1 + A_2_2),(B_1_1));
        Matrix p3 = Strassen((A_1_1),(B_1_2-B_2_2));
        Matrix p4 = Strassen((A_2_2),(B_2_1-B_1_1));
        Matrix p5 = Strassen((A_1_1 + A_1_2),(B_2_2));
        Matrix p6 = Strassen((A_2_1 - A_1_1),(B_1_1+B_1_2));
        Matrix p7 = Strassen((A_1_2 - A_2_2),(B_2_1+B_2_2));
        Matrix c_1_1 = p1+p4-p5+p7;
        Matrix c_1_2 = p3+p5;
        Matrix c_2_1 = p2+p4;
        Matrix c_2_2 = p1+p3-p2+p6;

        Matrix A_mult_B(rows_num_w_padding, rows_num_w_padding);
        for (int i = 0; i < rows_num_submatrix; ++i) {
            for (int j = 0; j < rows_num_submatrix; ++j) {
                A_mult_B[i][j] = c_1_1[i][j];
            }
        }
        for (int i = 0; i < rows_num_submatrix; ++i) {
            for (int j = rows_num_submatrix; j < rows_num_w_padding; ++j) {
                A_mult_B[i][j] = c_1_2[i][j - rows_num_submatrix];
            }
        }
        for (int i = rows_num_submatrix; i < rows_num_w_padding; ++i) {
            for (int j = 0; j < rows_num_submatrix; ++j) {
                A_mult_B[i][j] = c_2_1[i - rows_num_submatrix][j];
            }
        }
        for (int i = rows_num_submatrix; i < rows_num_w_padding; ++i) {
            for (int j = rows_num_submatrix; j < rows_num_w_padding; ++j) {
                A_mult_B[i][j] = c_2_2[i-rows_num_submatrix][j-rows_num_submatrix];
            }
        }

        for (int i = 0; i < C.rows_number; ++i) {
            for (int j = 0; j < C.columns_number; ++j) {
                C[i][j] = A_mult_B[i][j];
            }
        }

        return C;
    }
};



std::ostream& operator<<(std::ostream& out, const Matrix& m)
{
    if (out) {
        out << std::fixed << std::setprecision(5);
        for (int r = 0; r < m.rows_number; ++r) {
            for (int c = 0; c < m.columns_number; ++c) {
                out << (c > 0 ? " " : "") << std::setw(4);
                out << m.matrix[r][c];
            }
            out << std::endl;
        }
    }
    return out;
}

int main() {

    //
    //first simple check
    //

    Matrix X1(12,17);
    Matrix X2(17,11);
    for (int i = 0; i < X1.rows_number; ++i) {
        for (int j = 0; j < X1.columns_number; ++j) {
            X1[i][j] = i*i + 2*j + 1;
        }
    }
    for (int i = 0; i < X2.rows_number; ++i) {
        for (int j = 0; j < X2.columns_number; ++j) {
            X2[i][j] = i + j + 2;
        }
    }
    //X1*X2 using Strassen multiplication method
    Matrix Strassen_mult_0 = Matrix::Strassen(X1, X2);
    //A*B using regular multiplication method
    Matrix basic_mult_0 = X1 * X2;


    std::cout << X1 << std::endl;
    std::cout << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl;
    std::cout << X2 << std::endl;
    std::cout << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl;
    std::cout << Strassen_mult_0 << std::endl;
    std::cout << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl;
    std::cout << X1*X2 << std::endl;

    //compare results
    //if a[i][j] != b[i][j] print false
    for (int i = 0; i < Strassen_mult_0.rows_number; ++i) {
        for (int j = 0; j < Strassen_mult_0.columns_number; ++j) {
            if (Strassen_mult_0[i][j] != basic_mult_0[i][j])
                std::cout<<"false"<<std::endl;
        }
    }



    //
    // Now check with bigger matrices
    //


    //create matrices
    Matrix A(154,71);
    Matrix B(71,117);
    for (int i = 0; i < A.rows_number; ++i) {
        for (int j = 0; j < A.columns_number; ++j) {
            A[i][j] = i + j + 1;
        }
    }
    for (int i = 0; i < B.rows_number; ++i) {
        for (int j = 0; j < B.columns_number; ++j) {
            B[i][j] = i + j + 2;
        }
    }
    //A*B using Strassen multiplication method
    Matrix Strassen_mult = Matrix::Strassen(A, B);
    //A*B using regular multiplication method
    Matrix basic_mult = A * B;

    //in online compiler on my computer it takes approximately 27s to print all these matrices
    std::cout << A << std::endl;
    std::cout << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl;
    std::cout << B << std::endl;
    std::cout << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl;
    std::cout << Strassen_mult << std::endl;
    std::cout << "''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''" << std::endl;
    std::cout << A*B << std::endl;

    //compare results
    //if a[i][j] != b[i][j] print false
    for (int i = 0; i < Strassen_mult.rows_number; ++i) {
        for (int j = 0; j < Strassen_mult.columns_number; ++j) {
            if (Strassen_mult[i][j] != basic_mult[i][j])
                std::cout<<"false"<<std::endl;
        }
    }


    //check for very big matrices without printing
    //in online compiler on my computer it takes approximately 35s to calculate products
    // and check if they are identical
    Matrix C(15723,15726);
    Matrix D(15726,15729);
    for (int i = 0; i < C.rows_number; ++i) {
        for (int j = 0; j < C.columns_number; ++j) {
            C[i][j] = i + j + 1;
        }
    }
    for (int i = 0; i < D.rows_number; ++i) {
        for (int j = 0; j < D.columns_number; ++j) {
            D[i][j] = i + j + 2;
        }
    }

    Matrix Strassen_mult_big_check = Matrix::Strassen(A, B);

    Matrix basic_mult_big_check = A * B;

    for (int i = 0; i < Strassen_mult.rows_number; ++i) {
        for (int j = 0; j < Strassen_mult.columns_number; ++j) {
            if (Strassen_mult[i][j] != basic_mult[i][j])
                std::cout<<"false"<<std::endl;
        }
    }

    return 0;
}
