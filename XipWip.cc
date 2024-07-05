#include "iitpave2.h"

MatrixXd numint_coeff(int N, int n, MatrixXd Xp, MatrixXd Xl, double H,float a){
   myfunctions my;
 // Number of load points (xl) and deformation points (xd)
    int xl = Xl.rows(); // Number of load points (rows in Xl)
    int xp = Xp.rows(); // Number of deformation points (rows in Xd)
 
    MatrixXd Alpha = my.repmat((a/H), xp, 1);
    Xl    = my.repmat(Xl, xp, 1);
    Xp    = my.repmat(Xp, xl, 1);  //Xd is reorganized below to correspond to the
    
   RowVectorXd temp_arg_rowvec = my.LinSeq(1,xp).transpose();
   MatrixXd dx = my.repmat(temp_arg_rowvec,xl,1);
 
   MatrixXd dx1 = dx + my.repmat(my.LinSeq(0,xp,xp*(xl-1)),xp,1).transpose();
   dx = dx1;
   MatrixXd Xpp = my.getSubMatrix(Xp, dx.transpose()); // using transpose() for row selection takes a vector input (column)
   Xp = Xpp;
   MatrixXd Arr = my.GetColumn(Xp, ":", 1) - my.GetColumn(Xl, ":", 1);
   VectorXd r = ((Xp.col(0) - Xl.col(0)).array().square() + (Xp.col(1) - Xl.col(1)).array().square()).sqrt();
   //Introduce the parameter rho (notice: since r is a vector, rho is also  a vector)
   MatrixXd rho = r/H;
      
   MatrixXd B0r = besselroots(0, N, 1);
   MatrixXd Blr = besselroots(1, N, 1);
   


   #if 0
   // Example usage of besselroots function
    int o = 0;
    int N = 10;
    int k = 1;
 
    cout << "Bessel roots: " << roots << endl;
    
 
    cout << "Matrix rho:\n" << rho << endl;
  

   #endif
   return MatrixXd::Identity(1, 1);

}

