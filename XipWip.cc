#include "iitpave2.h"

MatrixXd numint_coeff(int N, int n, MatrixXd Xp, MatrixXd X1, double H,float a){

 // Number of load points (xl) and deformation points (xd)
    int xl = X1.rows(); // Number of load points (rows in Xl)
    int xd = Xp.rows(); // Number of deformation points (rows in Xd)
  
    MatrixXd mat(1,1);
    mat << a/H;

    // Repeat the matrixAlpha
    MatrixXd Alpha = repmat(mat, xd, 1);
     mat << 3.33; // converting double to Xd
    X1    = repmat(X1, xd, 1);

    // Output the result
    cout << Alpha << endl << "printed Alpha" << endl;
    cout << X1 << endl << "printed X1" << endl;

    return MatrixXd::Identity(1, 1);
}

