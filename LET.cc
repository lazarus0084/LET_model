#include "iitpave2.h"

void init_LET(Pave& iitpave) {
    int N = iitpave.N;   // Number of Bessel roots in integration
    int n = iitpave.n;  // Number of integration points in between each Bessel root
    float a = 150.8;      // MM
    MatrixXd Xp = iitpave.Xp; // Output coordinates !!Xd in alva file
    MatrixXd Xl = iitpave.Xl; // Transpose Load position coordinates
   //MatrixXd a = iitpave.a.transpose();   // Transpose Load radii
    MatrixXd E  = iitpave.E.transpose();   // Transpose Layer Young's moduli
    MatrixXd nu = iitpave.nu.transpose(); // Transpose Layer Poisson's ratios
    MatrixXd zi = iitpave.zi.transpose(); // // Define H as the last element of zi

    
   double H = zi(zi.size() - 2); //Thickness of last layer above subgrade
   double Laml = zi(0)/H ; //Relative height of top layer
   iitpave.Laml = Laml ; iitpave.H = H ;
   
   MatrixXd XipWip = numint_coeff(N, n, Xp, Xl, H , a);  //*filename XipWip sIntegration points and weights
// Print matrices with clear labels and separation
#if 0
cout << "Xp:" << endl << Xp << endl << endl;
cout << "Printed Xp" << endl << "-----------------------------" << endl;

cout << "Xl:" << endl << Xl << endl << endl;
cout << "Printed X1" << endl << "-----------------------------" << endl;

cout << "N:" << endl << N << endl << endl;
cout << "Printed N" << endl << "-----------------------------" << endl;

cout << "n:" << endl << n << endl << endl;
cout << "Printed n" << endl << "-----------------------------" << endl;

cout << "a:" << endl << a << endl << endl;
cout << "Printed a" << endl << "-----------------------------" << endl;

cout << "nu:" << endl << nu << endl << endl;
cout << "Printed nu" << endl << "-----------------------------" << endl;

cout << "zi:" << endl << zi << endl << endl;
cout << "Printed zi" << endl << "-----------------------------" << endl;

cout << "E:" << endl << E << endl << endl;
cout << "Printed E" << endl << "-----------------------------" << endl;
cout << "Laml:" << endl << Laml << endl << endl;
cout << "Printed Laml" << endl << "-----------------------------" << endl;
#endif 


}

    




