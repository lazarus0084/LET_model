#include "iitpave2.h"

void init_LET(Pave& iitpave) {
    int N = iitpave.N;   // Number of Bessel roots in integration
    int n = iitpave.n;  // Number of integration points in between each Bessel root
    float a = 310;      // MM
    MatrixXd Xp = iitpave.Xp; // Output coordinates !!Xd in alva file
    MatrixXd X1 = iitpave.X1; // Transpose Load position coordinates
   //MatrixXd a = iitpave.a.transpose();   // Transpose Load radii
    MatrixXd E = iitpave.E;   // Transpose Layer Young's moduli
    MatrixXd nu = iitpave.nu; // Transpose Layer Poisson's ratios
    MatrixXd zi = iitpave.zi; // // Define H as the last element of zi
   
    
   double H = zi(zi.size() - 2); //Thickness of last layer above subgrade
   double Lam1 = zi(0)/H ; //Relative height of top layer
   iitpave.Lam1 = Lam1 ; iitpave.H = H ;
   
   MatrixXd XipWip = numint_coeff(N, n, Xp, X1, H , a);  //*filename XipWip sIntegration points and weights
   
        


}

    




