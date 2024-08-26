#include "iitpave2.h"

void init_LET(Pave& iitpave) {
    int N = iitpave.N;   // Number of Bessel roots in integration
    int n = iitpave.n; 
    double a = 150.8;     // mm
    iitpave.a = a;
    MatrixXd Xp = iitpave.Xp; // Evaluation coordinates data
    MatrixXd Xl = iitpave.Xl; // Load application position coordinates
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
 iitpave.XipWip = XipWip;
// Display size of XipWip
// cout << "Size of XipWip: " << XipWip.rows() << " rows x " << XipWip.cols() << " columns." << endl;
// cout << "Element at (4, 6): " << XipWip(9, 42) << endl;
MatrixXd ABCD = arb_func(E.size(),zi,E,nu,iitpave) ;
iitpave.ABCD = ABCD;
LET_response(iitpave);

}

    




