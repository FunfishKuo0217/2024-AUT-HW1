make#include "algebra.h"
#include <stdexcept>
#include <random>
#include <vector>

namespace algebra {

    template<typename T>
    MATRIX<T> create_matrix(std::size_t rows, std::size_t columns, std::optional<MatrixType> type, std::optional<T> lowerBound, std::optional<T> upperBound) {
        // MATRIX<T> matrix(rows, vector<T>(columns, 0));
        // invalid dimension
        if(rows <= 0 || columns <= 0){throw std::invalid_argument("Invalid dimension.");}
        MATRIX<T> matrix(rows, std::vector<T>(columns, 0));

        if(type == MatrixType::Ones){
            for(size_t i = 0; i < rows; i++){
                for(size_t j = 0; j < columns; j++){
                    matrix[i][j] = 1;
                }
            }
        }
        else if(type == MatrixType::Identity){
            if(rows != columns){throw std::invalid_argument("Identity matrix must be square (rows must equal columns).");}
            for(size_t i = 0; i < rows; i++){
              for(size_t j = 0; j < columns; j++){
                if(i == j) matrix[i][j] = 1.0;
              }
            }
        }
        else if(type == MatrixType::Random){
            if(!lowerBound || !upperBound){throw std::invalid_argument("Lowerbound and upperbound must be provided.");}
            if(lowerBound >= upperBound){throw std::invalid_argument("Lowerbound should be smaller.");}
            for(size_t i = 0; i < rows; i++){
                for(size_t j = 0; j < columns; j++){
                    matrix[i][j] = randm(lowerBound, upperBound);
                }
            }
        }

        return matrix;
    }

    template<typename T>
    MATRIX<T> sum_sub(const MATRIX<T>& matrixA, const MATRIX<T>& matrixB, std::optional<std::string> operation){
        size_t rowsA = matrixA.size(); 
        size_t rowsB = matrixB.size();
        if(rowsA == 0 || rowsB == 0){return {};}
        size_t columnsA = matrixA[0].size(); 
        size_t columnsB = matrixB[0].size();
        if(rowsA != rowsB || columnsA != columnsB){throw std::invalid_argument("MatA and MatB are mismatched!");}  
        
        MATRIX<T> matrix(rowsA, std::vector<T>(columnsA, 0));
        if (operation == "sub"){
            for(size_t i = 0; i < rowsA; i++){
                for(size_t j = 0; j < columnsA; j++){
                    matrix[i][j] = matrixA[i][j] - matrixB[i][j];  
                }
            }
        }
        else{
            for(size_t i = 0; i < rowsA; i++){
                for(size_t j = 0; j < columnsA; j++){
                    matrix[i][j] = matrixA[i][j] + matrixB[i][j];  
                }
            }
        }
        return matrix;
    }
    
    
    template<typename T>
    MATRIX<T> multiply(const MATRIX<T>& matrix, const T scalar) {
      // if(matrix.size() == 0){throw std::invalid_argument("Invalid matrix");}
      MATRIX<T> result = matrix;
      for(auto& row : result){
          for(auto& element : row){
              element *= scalar;
          }
      }
      return result;
    }

    template<typename T>
    MATRIX<T> multiply(const MATRIX<T>& matrixA, const MATRIX<T>& matrixB){
        size_t rowA = matrixA.size();
        size_t rowB = matrixB.size();
        if(rowA == 0 || rowB == 0){throw std::invalid_argument("Invalid matrix");}
        
        size_t colA = matrixA[0].size();
        size_t colB = matrixB[0].size();

        // check mismatched
        if(colA != rowB){throw std::invalid_argument("Mismatched matrix!");}
        
        size_t result_row = rowA;
        size_t result_col = colB;
        MATRIX<T> result(result_row, std::vector<T>(result_col, 0));
        
        for(size_t r = 0; r < result_row; r++){
            for(size_t c = 0; c < result_col; c++){
                // caculation
                for(size_t k = 0; k < rowB; k++){
                    result[r][c] += matrixA[r][k]*matrixB[k][c];
                }
            }
        }
        return result;
    }
    template<typename T>
    MATRIX<T> hadamard_product(const MATRIX<T>& matrixA, const MATRIX<T>& matrixB){
        if(matrixA.size() == 0 || matrixB.size() == 0){return {};}
        if(matrixA.size() != matrixB.size() || matrixA[0].size() != matrixB[0].size()){throw std::invalid_argument("Mismatched matrix!");}

        size_t result_row = matrixA.size();
        size_t result_col = matrixA[0].size();
        MATRIX<T> result = matrixA;
        for(size_t r = 0; r < result_row; r++){
            for(size_t c = 0; c < result_col; c++){
                result[r][c] *= matrixB[r][c];
            }
        }
        return result;
    }
    
    template<typename T>
    MATRIX<T> transpose(const MATRIX<T>& matrix){
        size_t row = matrix.size();
        if(row == 0){return {};}
        size_t col = matrix[0].size();
        MATRIX<T> result(col, std::vector<T>(row, 0));
        for(size_t r = 0; r < col; r++){
            for(size_t c = 0; c < row; c++){
                result[r][c] = matrix[c][r];
            }
        }
        return result;
    }

    template<typename T>
    T trace(const MATRIX<T>& matrix){
        //check square matrix
        if(matrix.size() == 0 || matrix.size() != matrix[0].size()){throw std::invalid_argument("Invalid matrix");}
        T trace = 0;
        for(size_t i = 0; i < matrix.size(); i++){
            trace += matrix[i][i];
        }
        return trace;
    }

    // helper function
    double** MatrixConstruct(int n, int m) {
        double** A = new double*[n];
        if (!A)
            return NULL;
        for (int i = 0; i < n; ++i) {
            A[i] = new double[m];
            if (!A[i])
                return NULL;
        }
        return A;
    }

    void DeleteMatrix(double** A, int n) {
        for (int i = 0; i < n; ++i)
            delete[] A[i];
        delete[] A;
        A = NULL;
    }


    double Det(double** A, int n){
        if(n == 1){return A[0][0];}
        if(n == 2){return A[0][0]*A[1][1] - A[0][1]*A[1][0];}
        double det = 0;
        double** M = MatrixConstruct(n-1, n-1);
        
        for(int a = 0; a < n; a++){
            for(int i = 1, minor_i = 0; i < n; ++i, ++minor_i){
                for(int j = 0, minor_j = 0; j < n; ++j){
                    if(j == a){continue;}
                    M[minor_i][minor_j++] = A[i][j];
                }
            }

            // add up term
            if(a%2 == 0){det += A[0][a] * Det(M, n-1);}
            else{det -= A[0][a] * Det(M, n-1);}
        }

        DeleteMatrix(M, n-1);
        return det;
    }

    template<typename T>
    double determinant(const MATRIX<T>& matrix){
        int sz = matrix.size();
        if(sz == 0 || (size_t) sz != matrix[0].size()){throw std::invalid_argument("Invalid matrix");}
        double** A = MatrixConstruct(sz, sz);
        for(int i = 0; i < sz; i++){
            for(int j = 0; j < sz; j++){
                A[i][j] = (double) matrix[i][j];
            }
        }
        return Det(A, sz);
    }

    // instantialization
    template MATRIX<double> create_matrix<double>(unsigned long, unsigned long, 
                                                  std::optional<MatrixType>, 
                                                  std::optional<double>, 
                                                  std::optional<double>);
    template MATRIX<float> create_matrix<float>(unsigned long, unsigned long, 
                                                  std::optional<MatrixType>, 
                                                  std::optional<float>, 
                                                  std::optional<float>);
    template MATRIX<int> create_matrix<int>(unsigned long, unsigned long, 
                                                  std::optional<MatrixType>, 
                                                  std::optional<int>, 
                                                  std::optional<int>);

    template MATRIX<double> sum_sub<double>(const MATRIX<double>& matrixA, 
                                                              const MATRIX<double>& matrixB, 
                                                              std::optional<std::string> operation);
    template MATRIX<float> sum_sub<float>(const MATRIX<float>& matrixA, 
                                                              const MATRIX<float>& matrixB, 
                                                              std::optional<std::string> operation);
    template MATRIX<int> sum_sub<int>(const MATRIX<int>& matrixA, 
                                                              const MATRIX<int>& matrixB, 
                                                              std::optional<std::string> operation);

    template MATRIX<double> multiply<double>(const MATRIX<double>& matrix, const double scalar);
    template MATRIX<float> multiply<float>(const MATRIX<float>& matrix, const float scalar);
    template MATRIX<int> multiply<int>(const MATRIX<int>& matrix, const int scalar);
    
    template MATRIX<double> multiply<double>(const MATRIX<double>& matrixA, const MATRIX<double>& matrixB);
    template MATRIX<float> multiply<float>(const MATRIX<float>& matrixA, const MATRIX<float>& matrixB);
    template MATRIX<int> multiply<int>(const MATRIX<int>& matrixA, const MATRIX<int>& matrixB);  

    template MATRIX<double>  hadamard_product<double>(const MATRIX<double>& matrixA, const MATRIX<double>& matrixB);
    template MATRIX<float> hadamard_product<float>(const MATRIX<float>& matrixA, const MATRIX<float>& matrixB);
    template MATRIX<int> hadamard_product<int>(const MATRIX<int>& matrixA, const MATRIX<int>& matrixB);

    template MATRIX<double> transpose<double>(const MATRIX<double>& matrix);
    template MATRIX<float> transpose<float>(const MATRIX<float>& matrix);
    template MATRIX<int> transpose<int>(const MATRIX<int>& matrix);
    
    template double trace<double>(const MATRIX<double>& matrix);
    template float trace<float>(const MATRIX<float>& matrix);
    template int trace<int>(const MATRIX<int>& matrix);

    template double determinant<double>(const MATRIX<double>& matrix);
    template double determinant<float>(const MATRIX<float>& matrix);
    template double determinant<int>(const MATRIX<int>& matrix);
}