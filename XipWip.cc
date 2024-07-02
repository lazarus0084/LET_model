#include "iitpave2.h"

MatrixXd numint_coeff(int N, int n, MatrixXd Xp, MatrixXd X1, double H,float a){
   myfunctions my;
 // Number of load points (xl) and deformation points (xd)
    int x1 = X1.rows(); // Number of load points (rows in Xl)
    int xd = Xp.rows(); // Number of deformation points (rows in Xd)

    MatrixXd Alpha = my.repmat((a/H), xd, 1);
    X1    = my.repmat(X1, x1, 1);
    Xp    = my.repmat(Xp, xd, 1);  //Xd is reorganized below to correspond to the
    MatrixXd dx = my.repmat(my.LinSeq(1,xd),x1,1);
    MatrixXd dx1 = (dx + my.repmat(my.LinSeq(0,xd,xd*(x1-1)),xd,1)).transpose();
    dx = dx1;
    cout << "Matrix Xp:\n" << Xp << endl;
    cout << "Matrix dx:\n" << dx << endl;


    MatrixXd Xpp = my.getSubMatrix(Xp, dx.transpose());
    Xp = Xpp;

    // Display the extracted matrix
    cout << "Extracted Matrix:\n" << Xp << endl;

    

    
    return MatrixXd::Identity(1, 1);
}

