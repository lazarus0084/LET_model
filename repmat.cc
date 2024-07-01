#include "iitpave2.h"

MatrixXd repmat(const MatrixXd& mat, int rows, int cols) {
    MatrixXd result(mat.rows() * rows, mat.cols() * cols);
    for(int i = 0; i < rows; ++i) {
        for(int j = 0; j < cols; ++j) {
            result.block(i * mat.rows(), j * mat.cols(), mat.rows(), mat.cols()) = mat;
        }
    }
    return result;
}

