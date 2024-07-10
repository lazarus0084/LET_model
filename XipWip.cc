#include "iitpave2.h"

MatrixXd numint_coeff(int N, int n, MatrixXd Xp, MatrixXd Xl, double H,float a){
   
 // Number of load points (xl) and deformation points (xd)
    int xl = Xl.rows(); // Number of load points (rows in Xl)
    int xp = Xp.rows(); // Number of deformation points (rows in Xd)
 
    MatrixXd alpha = repmat((a/H), xp, 1);
    Xl    = repmat(Xl, xp, 1);
    Xp    = repmat(Xp, xl, 1);  //Xd is reorganized below to correspond to the
    
   RowVectorXd temp_arg_rowvec = LinSeq(1,xp).transpose();
   MatrixXd dx = repmat(temp_arg_rowvec,xl,1);
 
   MatrixXd dx1 = dx + repmat(LinSeq(0,xp,xp*(xl-1)),xp,1).transpose();
   dx = dx1;
   MatrixXd Xpp = getSubMatrix(Xp, dx.transpose()); // using transpose() for row selection takes a vector input (column)
   Xp = Xpp;
   MatrixXd Arr = GetColumn(Xp, ":", 1) - GetColumn(Xl, ":", 1);
   VectorXd r = ((Xp.col(0) - Xl.col(0)).array().square() + (Xp.col(1) - Xl.col(1)).array().square()).sqrt();
   //Introduce the parameter rho (notice: since r is a vector, rho is also  a vector)
   VectorXd rho = r/H; 
   MatrixXd B0r = besselroots(0, N, 1);
   MatrixXd B1r = besselroots(1, N, 1);
   double n3 = rho.size();       //Number of radii
   double n4 = alpha.size();     //Number of loads: NB: n3 = n4, due to repmat(..)
                                 //commands above

    // Find non-zero radii/rho
    vector<int> rho_Non0;
    for (int i = 0; i < rho.size(); ++i) {
        if (rho(i) != 0) {
            rho_Non0.push_back(i);
        }
    }
    MatrixXd mat1;
    mat1.resize(n3, N); mat1.setZero();
    MatrixXd mat2 = repmat(B1r,n4,1).array()/repmat(alpha,1,N).array();
    MatrixXd roots_z = JoinAllRows(mat1,mat2);
   
    
    VectorXd arg1 = rho(rho_Non0);
    VectorXi columnIndices = VectorXi::LinSpaced(N, 0, N-1);
    //rho_Non0 is the rowIndices!
    if (rho_Non0.size() != 0) {
      roots_z(rho_Non0, columnIndices) = (repmat(B0r,rho_Non0.size(),1)).array()/(repmat(arg1,1,N)).array();
    }
    // Sort each column (which is now each row after transposition)
    roots_z.transposeInPlace();
   
    for (int i = 0; i < roots_z.rows(); ++i) {
        sort(roots_z.row(i).data(), roots_z.row(i).data() + roots_z.cols());
    }
    // Transpose back to original orientation
    roots_z.transposeInPlace();







         // Print rows and columns of roots_z
    cout << "roots_z dimensions: (" << roots_z.rows() << ", " << roots_z.cols() << ")" << endl;
    cout << "Element of roots_z after sorting at position (7,599): " << roots_z(7,599) << endl;
    
    return MatrixXd::Identity(1, 1);
   
}

