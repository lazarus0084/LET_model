#include "iitpave2.h"

void arb_func(int n, const VectorXd& zi_temp, const VectorXd& E, const VectorXd& nu, const Pave& iitpave) {
    
//------------------------------------------------------------------------------
// DESCRIPTION:
// This function evaluates the coefficients of integration Ai, Bi, Ci and Di
// of each layer. These are unitless functions that embody the layered 
// system properties and connectivity.

// INPUT PARAMETERS
// n    : Number of layers (including half space layer)
// zi   : Distance from surface to the bottom of each layer
// E    : Layer Young's moduli
// nu   : Layer Poissons ratio's
// bond : Interface bonding type
//------------------------------------------------------------------------------

// To speed up the computational time, the number of matrix inversions
// is limited to 96, corresponding to 96 predetermined values of the
// integration variable m in the range of 0 to 100,000, as follows:

    // Print the initial zi vector
 
    VectorXd zi(2);
    zi << 260,760;
   

RowVectorXd m(96);
m << 1e-10, 0.05, 0.10, 0.20, 0.40, 0.60, 0.80, 1.00, 1.20, 1.40, 1.60,
     1.80, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.20, 3.40, 3.60, 3.80,
     4.00, 4.20, 4.40, 4.60, 4.80, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50,
     8.00, 8.50, 9.00, 9.50, 10.00, 11.00, 12.00, 13.00, 14.00, 15.00,
     16.00, 17.00, 18.00, 19.00, 20.00, 25.00, 30.00, 35.00, 40.00,
     45.00, 50.00, 55.00, 60.00, 65.00, 70.00, 75.00, 80.00, 85.00,
     90.00, 95.00, 100.00, 110.00, 120.00, 130.00, 140.00, 150.00,
     160.00, 170.00, 180.00, 190.00, 200.00, 210.00, 220.00, 230.00,
     240.00, 250.00, 260.00, 270.00, 280.00, 290.00, 300.00, 350.00,
     400.00, 450.00, 500.00, 600.00, 700.00, 800.00, 900.00, 1000.00,
     2000.00, 10000.00, 100000.00;
//Number of integration points
int Lm   = m.size();               // Lm = Length of m
RowVectorXd ONES(Lm); ONES.setOnes();

// Evaluate the parameters ABCD
double H = zi(zi.size()-1);       // Depth to the upper boundary of the 
                                  //inifnite (bottom) layer from the top 
                                // surface 
VectorXd lam(zi.size() + 2);      // Ensure `lam_temp` is sized appropriately
lam << 0, zi, 1e20 * zi.maxCoeff(); // Include zi and the max coefficient scaled by 1e20
lam = lam/H;

//VectorXd lam = lam_temp.transpose(); // Transpose to get a column vector                             // the number of layer interfaces (including
                                // the surface of the top layer). We have 
                                // added a infinite value (1e20*max(zi)) in 
                                //the last place. But it will not influence 
                                // on the outcome, since it will only be 
                                // used in the last A and C coefficients 
                                // equations which are remmoved from the 
                                // system of equations before solved, since 
                                // these A and C values are zero. These two 
                                // equations are however organized, as it is 
                                // easier to do in the loop instead of 
                                // having a special case at the end.

//lam = [lam0 lam1 lam2 ... lam_n-1 0] = [z0 z1 z2 ... z_n-1 0]/H
//I.e. lami = lam(i+1)


if (zi.size() > E.size()){
    cout << "The length of zi is too high or the length of E is to short";
}
//--------------------------------------------------------------------------
// Setup the boundary and continuity conditions
//--------------------------------------------------------------------------
// Define the i and j-coordinates
RowVectorXd j = RowVectorXd::LinSpaced(n-1,0,n-2);
RowVectorXd i = j.array() + 1;


// Define the equations
MatrixXd Ai = VectorXd::Ones(n-1)*ONES;
//Equations based on the vertical stresses
// sigmaz = [                                                    ... 
//  1*ones(n-1,1)*ONES                                           ...   Ai
//  exp(-(lam(i)-lam(i-1))*m)                                    ...  % Bi
// -(1-2*nu(j)*ONES-lam(i)*m)                                    ...  % Ci
//  (1-2*nu(j)*ONES+lam(i)*m).*exp(-(lam(i)-lam(i-1))*m)         ...  % Di
// -exp(-(lam(i+1)-lam(i))*m)                                    ...  % A{i+1}
// -1*ones(n-1,1)*ONES                                           ...  % B{i+1}
//  (1-2*nu(j+1)*ONES-lam(i)*m).*exp(-(lam(i+1)-lam(i))*m)       ...  % C{i+1}
// -(1-2*nu(j+1)*ONES+lam(i)*m)                                  ...  % D{i+1}
// ];
cout << Ai;

 }




