#include "iitpave2.h"
MatrixXd repmat(int mat, int rows, int cols) {
    return MatrixXd::Constant(rows, cols, static_cast<double>(mat));
}

// Function overload for double
MatrixXd repmat(double mat, int rows, int cols) {
    return MatrixXd::Constant(rows, cols, mat);
}

// Function overload for VectorXd
MatrixXd repmat(const VectorXd& mat, int rows, int cols) {
    int size = mat.size();
    MatrixXd result(size * rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.block(i * size, j, size, 1) = mat;
        }
    }

    return result;
}


// Function overload for RowVectorXd
MatrixXd repmat(const RowVectorXd& mat, int rows, int cols) {
    int size = mat.size();
    MatrixXd result(rows, cols * size);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.block(i, j * size, 1, size) = mat;  // Use block operation to assign the row vector
        }
    }

    return result;
}

// Function overload for MatrixXd
MatrixXd repmat(const MatrixXd& mat, int rows, int cols) {
    int numRows = mat.rows();
    int numCols = mat.cols();
    MatrixXd result(numRows * rows, numCols * cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.block(i * numRows, j * numCols, numRows, numCols) = mat;
        }
    }

    return result;
}

// Function to generate a sequence j:k
MatrixXd LinSeq(int j, int k) {
    int m = k - j + 1;
    MatrixXd result(m, 1);
    for (int idx = 0; idx < m; ++idx) {
        result(idx, 0) = j + idx;
    }
    return result;
}

// Function to generate a sequence j:i:k
MatrixXd LinSeq(int j, int i, int k) {
    int m = static_cast<int>((k - j) / i) + 1;
    MatrixXd result(m, 1);
    for (int idx = 0; idx < m; ++idx) {
        result(idx, 0) = j + idx * i;
    }
    return result;
}
MatrixXd getSubMatrix(const MatrixXd& Mat, const VectorXd& indices) {
    int numRows = indices.size();
    int numCols = Mat.cols();
    MatrixXd submatrix(numRows, numCols);
    submatrix.setZero();

    VectorXd adjustedIndices = indices.array() - 1; // Adjust indices by subtracting 1

    // Extract rows based on adjustedIndices
    for (int i = 0; i < numRows; ++i) {
        int rowIdx = static_cast<int>(adjustedIndices(i));
        if (rowIdx >= 0 && rowIdx < Mat.rows()) {
            submatrix.row(i) = Mat.row(rowIdx);
        } else {
            cerr << "Row index " << rowIdx << " out of bounds." << endl;
        }
    }

    return submatrix;
}
MatrixXd GetColumn(const MatrixXd& Xp, const std::string& row_index, int column_index) {
    int rows = Xp.rows();
    int cols = Xp.cols();
    MatrixXd result;

    if (row_index == ":") {
        // Get all rows for a specific column
        result.resize(rows, 1);
        result = Xp.col(column_index - 1);  // Adjust column index by -1
    } else {
        // Get a specific row and column
        int row = std::stoi(row_index) - 1; // Adjust row index by -1
        result.resize(1, 1);
        result(0, 0) = Xp(row, column_index - 1);  // Adjust column index by -1
    }

    return result;
}
// Overload for concatenating MatrixXd with MatrixXd
MatrixXd JoinAllRows(const MatrixXd& mat1, const MatrixXd& mat2) {
    // Check if matrices have the same number of rows
    assert(mat1.rows() == mat2.rows());

    // Create a new matrix for concatenation
    MatrixXd concatenated(mat1.rows(), mat1.cols() + mat2.cols());

    // Copy mat1 into the left part of concatenated matrix
    concatenated.block(0, 0, mat1.rows(), mat1.cols()) = mat1;

    // Copy mat2 into the right part of concatenated matrix
    concatenated.block(0, mat1.cols(), mat2.rows(), mat2.cols()) = mat2;

    return concatenated;
}

// Overload for concatenating VectorXd with MatrixXd
MatrixXd JoinAllRows(const VectorXd& vec1, const MatrixXd& mat2) {
    // Check if vector and matrix have the same number of rows
    assert(vec1.size() == mat2.rows());

    // Create a new matrix for concatenation
    MatrixXd concatenated(vec1.size(), 1 + mat2.cols());

    // Copy vec1 into the first column of concatenated matrix
    concatenated.col(0) = vec1;

    // Copy mat2 into the remaining columns of concatenated matrix
    concatenated.block(0, 1, mat2.rows(), mat2.cols()) = mat2;

    return concatenated;
}
// Function to sort each column of the matrix
void sortColumns(MatrixXd& matrix) {
    int rows = matrix.rows();
    int cols = matrix.cols();

    for (int col = 0; col < cols; ++col) {
        // Extract the column into a temporary vector
        vector<double> column;
        for (int row = 0; row < rows; ++row) {
            column.push_back(matrix(row, col)); // Access using () for Eigen MatrixXd
        }

        // Sort the column vector
        sort(column.begin(), column.end());

        // Place sorted values back into the matrix column
        for (int row = 0; row < rows; ++row) {
            matrix(row, col) = column[row]; // Assign back using () for Eigen MatrixXd
        }
    }
}
VectorXd ExtractVecElements(const VectorXd& rho, const VectorXd& rho_non0) {
    VectorXd ExtractedElements(rho_non0.size());
    for (int i = 0; i < rho_non0.size(); ++i) {
        int index = static_cast<int>(rho_non0(i)) - 1; // Adjust to 0-based indexing
        ExtractedElements(i) = rho(index);
    }
    return ExtractedElements;
}
void saveit(const Eigen::MatrixXd& mat) {
    std::ofstream file("mat.txt");
    if (file.is_open()) {
        file << mat;
        file.close();
        std::cout << "Matrix saved to mat.txt" << std::endl;
    } else {
        std::cerr << "Unable to open file for writing" << std::endl;
    }
}

MatrixXd linearIndexing(const MatrixXd &A, const MatrixXd &B) {
    int rowsA = A.rows();
    int colsA = A.cols();
    int sizeB = B.size();

    MatrixXd C(B.rows(), B.cols());

    for (int i = 0; i < B.rows(); ++i) {
        for (int j = 0; j < B.cols(); ++j) {
            int index = static_cast<int>(B(i, j)) - 1;  // Converting to 0-based index
            if (index < 0 || index >= rowsA * colsA) {
                cerr << "Error: Index " << B(i, j) << " is out of bounds." << endl;
                exit(EXIT_FAILURE);
            }
            int row = index % rowsA;
            int col = index / rowsA;
            C(i, j) = A(row, col);
        }
    }

    return C;
}

