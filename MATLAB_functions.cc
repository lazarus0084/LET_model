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
// Function overload for RowVectorXi
MatrixXi repmat(const RowVectorXi& mat, int rows, int cols) {
    int size = mat.size();
    MatrixXi result(rows, cols * size);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.block(i, j * size, 1, size) = mat.cast<int>();  // Use block operation to assign the row vector
        }
    }

    return result;
}
// Function overload for VectorXi
MatrixXi repmat(const VectorXi& mat, int rows, int cols) {
    int size = mat.size();
    MatrixXi result(size * rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.block(i * size, j, size, 1) = mat;
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
MatrixXi repmat(const MatrixXi& mat, int rows, int cols) {
    int numRows = mat.rows();
    int numCols = mat.cols();
    MatrixXi result(numRows * rows, numCols * cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.block(i * numRows, j * numCols, numRows, numCols) = mat;
        }
    }

    return result;
}
// Function overload for SparseVector<double>
Eigen::SparseMatrix<double> repmat(const Eigen::SparseVector<double>& mat, int rows, int cols) {
    int size = mat.size();
    Eigen::SparseMatrix<double> result(size * rows, cols);

    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(mat.nonZeros() * rows * cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            for (Eigen::SparseVector<double>::InnerIterator it(mat); it; ++it) {
                tripletList.emplace_back(i * size + it.index(), j, it.value());
            }
        }
    }

    result.setFromTriplets(tripletList.begin(), tripletList.end());
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
MatrixXd LinSeq(int start, int step, int end) {
    int m = static_cast<int>((end - start) / step) + 1;
    MatrixXd result(m, 1);
    for (int idx = 0; idx < m; ++idx) {
        result(idx, 0) = start + idx * step ;
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
// function overload for saveit
void saveit(const Eigen::MatrixXi& mat) {
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
// Function to create an empty sparse matrix with reserved space
Eigen::SparseMatrix<double> spalloc(int m, int n, int nzmax) {
Eigen::SparseMatrix<double> s(m, n);
    s.reserve(nzmax);
    return s;
}


// Function for horizontal concatenation of two sparse matrices
SparseMatrix<double> H_SparseConcatenate(const SparseMatrix<double>& A, const SparseMatrix<double>& B) {
    if (A.rows() != B.rows()) {
        throw std::invalid_argument("Matrices must have the same number of rows for horizontal concatenation");
    }

    SparseMatrix<double> C(A.rows(), A.cols() + B.cols());

    // Insert elements from A
    for (int k = 0; k < A.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            C.insert(it.row(), it.col()) = it.value();
        }
    }

    // Insert elements from B
    for (int k = 0; k < B.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(B, k); it; ++it) {
            C.insert(it.row(), it.col() + A.cols()) = it.value();
        }
    }

    return C;
}

// Function for vertical concatenation of two sparse matrices
SparseMatrix<double> V_SparseConcatenate(const SparseMatrix<double>& A, const SparseMatrix<double>& B) {
    if (A.cols() != B.cols()) {
        throw std::invalid_argument("Matrices must have the same number of columns for vertical concatenation");
    }

    SparseMatrix<double> C(A.rows() + B.rows(), A.cols());

    // Insert elements from A
    for (int k = 0; k < A.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            C.insert(it.row(), it.col()) = it.value();
        }
    }

    // Insert elements from B
    for (int k = 0; k < B.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(B, k); it; ++it) {
            C.insert(it.row() + A.rows(), it.col()) = it.value();
        }
    }

    return C;
}
// Function to convert linear indices to row and column indices
MatrixXi ind2sub(int rows, int columns, const VectorXi& linear_indices) {
    // Size of the input vector
    int num_indices = linear_indices.size();

    // Vectors to store row and column indices
    VectorXi row_indices(num_indices);
    VectorXi col_indices(num_indices);

    // Compute row and column indices
    for (int i = 0; i < num_indices; ++i) {
        int idx = linear_indices(i);
        row_indices(i) = (idx - 1) % rows + 1; // MATLAB 1-based indexing for rows
        col_indices(i) = (idx - 1) / rows + 1; // MATLAB 1-based indexing for columns
    }

    // Create a matrix to hold row and column indices
    MatrixXi indices_matrix(num_indices, 2);
    indices_matrix.col(0) = row_indices;
    indices_matrix.col(1) = col_indices;

    return indices_matrix;
}

SparseMatrix<double> extractSparseSubMat(const SparseMatrix<double>& BC0, 
                                         const VectorXi& row_indices, 
                                         const VectorXi& column_indices) {
    // Determine the size of the submatrix
    int rows = row_indices.size();
    int cols = column_indices.size();
    
    // Create a sparse matrix to hold the submatrix
    SparseMatrix<double> submatrix(rows, cols);

    // Create a triplet list to construct the submatrix
    std::vector<Triplet<double>> tripletList;
    
    // Loop over the given row and column indices to extract the submatrix
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // Get the value from the original matrix
            double value = BC0.coeff(row_indices[i], column_indices[j]);
            if (value != 0.0) {
                // Add the non-zero value to the triplet list
                tripletList.emplace_back(i, j, value);
            }
        }
    }

    // Construct the sparse submatrix from the triplet list
    submatrix.setFromTriplets(tripletList.begin(), tripletList.end());
    
    return submatrix;
}
// Function to perform linear interpolation without using a class
RowVectorXd linearInterp(const RowVectorXd& x, const RowVectorXd& y, const RowVectorXd& xi) {
    // Ensure x and y are the same size
    if (x.size() != y.size()) {
        std::cerr << "Error: Size mismatch between x and y!" << std::endl;
        return RowVectorXd::Zero(xi.size());
    }

    // Sort the (x, y) pairs based on x
    std::vector<std::pair<double, double>> xyPairs(x.size());
    for (int i = 0; i < x.size(); ++i) {
        xyPairs[i] = {x(i), y(i)};
    }

    // Sort by x
    std::sort(xyPairs.begin(), xyPairs.end(), [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
        return a.first < b.first;
    });

    // Create sorted x and y vectors
    RowVectorXd xSorted(x.size()), ySorted(x.size());
    for (int i = 0; i < x.size(); ++i) {
        xSorted(i) = xyPairs[i].first;
        ySorted(i) = xyPairs[i].second;
    }

    // Preallocate the output
    RowVectorXd yi(xi.size());

    // Loop through each xi and perform interpolation
    for (int i = 0; i < xi.size(); ++i) {
        double xiValue = xi(i);

        // Find the interval in xSorted that contains xi
        auto it = std::lower_bound(xSorted.data(), xSorted.data() + xSorted.size(), xiValue);

        // Handle out-of-bounds cases
        if (it == xSorted.data()) {
            yi(i) = ySorted(0);  // Use first value if xi is below x range
        } else if (it == xSorted.data() + xSorted.size()) {
            yi(i) = ySorted(ySorted.size() - 1);  // Use last value if xi is above x range
        } else {
            // We're in the range: interpolate between the surrounding points
            int idx = it - xSorted.data();
            if (xSorted(idx) == xiValue) {
                yi(i) = ySorted(idx);  // xi matches an x value exactly
            } else {
                // Linear interpolation formula
                double xa = xSorted(idx - 1), xb = xSorted(idx);
                double ya = ySorted(idx - 1), yb = ySorted(idx);
                double t = (xiValue - xa) / (xb - xa);
                yi(i) = ya + t * (yb - ya);
            }
        }
    }

    return yi;
}
// Define the besselj function that takes an order `n` and a matrix `m1`
MatrixXd besselj(int n, const MatrixXd& m1) {
    MatrixXd result(m1.rows(), m1.cols());

    // Apply the Bessel function of order n to each element of the matrix
    for (int i = 0; i < m1.rows(); ++i) {
        for (int j = 0; j < m1.cols(); ++j) {
            if (n == 0) {
                result(i, j) = j0(m1(i, j));  // j0 for Bessel function of order 0
            } else if (n == 1) {
                result(i, j) = j1(m1(i, j));  // j1 for Bessel function of order 1
            } else {
                result(i, j) = jn(n, m1(i, j));  // jn for higher orders
            }
        }
    }

    return result;
}
// Custom cumsum function for Eigen vectors
VectorXd cumsum(const VectorXd& v) {
    VectorXd result(v.size());
    if (v.size() == 0) return result;

    result(0) = v(0);
    for (int i = 1; i < v.size(); ++i) {
        result(i) = result(i - 1) + v(i);
    }
    return result;
}