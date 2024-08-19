#include "iitpave2.h"

void init_LET(Pave& iitpave) {
    int N = iitpave.N;   // Number of Bessel roots in integration
    int n = iitpave.n;  // Number of integration points in between each Bessel root
    float a = 150.8;      // MM
    MatrixXd Xp = iitpave.Xp; // Output coordinates !!Xd in alva file
    MatrixXd Xl = iitpave.Xl; // Transpose Load position coordinates
   //MatrixXd a = iitpave.a.transpose();   // Transpose Load radii
    VectorXd E  = iitpave.E;   // Transpose Layer Young's moduli
    VectorXd nu = iitpave.nu; // Transpose Layer Poisson's ratios
    //RowVectorXd zi = iitpave.zi.transpose(); // // Define H as the last element of zi
    VectorXd zi(2);
    zi << 260,760;
    cout << "Forced populating zi LET.CC line 15." << endl ;
    double H = zi(zi.size() - 2); //Thickness of last layer above subgrade
    double Laml = zi(0)/ H ; //Relative height of top layer
    iitpave.Laml = Laml ; iitpave.H = H ;
   
 MatrixXd XipWip = numint_coeff(N, n, Xp, Xl, H , a) ;  //*filename XipWip sIntegration points and weights
//    // Display size of XipWip
// cout << "Size of XipWip: " << XipWip.rows() << " rows x " << XipWip.cols() << " columns." << endl;
// cout << "Element at (4, 6): " << XipWip(9, 42) << endl;
MatrixXd ABCD = arb_func(E.size(),zi,E,nu,iitpave) ;
cout << "Size of ABCD : " << ABCD.rows() << " X " << ABCD.cols() << endl ;
LET_response(iitpave);
}

    




