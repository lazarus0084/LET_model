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
RowVectorXi j = RowVectorXi::LinSpaced(n-1,0,n-2);
RowVectorXi i = j.array() + 1;


// Define the equations
MatrixXd Ai = VectorXd::Ones(n-1)*ONES;
MatrixXd Bi_temp = -(lam(i.array()) - lam(j.array())); Bi_temp = Bi_temp*m;
MatrixXd Bi = Bi_temp.array().exp();
MatrixXd Ci_a =  2*nu(j.array())*ONES; Ci_a = 1 - Ci_a.array() ;Ci_a = Ci_a - lam(i.array())*m;
MatrixXd Ci = - Ci_a;             //-(1-2*nu(j)*ONES-lam(i)*m)  
MatrixXd Di_a =  2*nu(j.array())*ONES; Di_a = 1 - Di_a.array() ;Di_a = Di_a + lam(i.array())*m; //(1-2*nu(j)*ONES+lam(i)*m).*exp(-(lam(i)-lam(i-1))*m)         ...  % Di
MatrixXd Di = Di_a.array() * Bi.array();
//for i+1 & j+1 coefficients
i = i.array() + 1;
j = j.array() + 1;
MatrixXd Aiplus1 = VectorXd::Ones(n-1)*ONES;
MatrixXd Biplus1_temp = -(lam(i.array()) - lam(j.array())); Biplus1_temp = Biplus1_temp*m;
MatrixXd Biplus1 = Biplus1_temp.array().exp();
MatrixXd Ciplus1_a =  2*nu(j.array())*ONES; Ciplus1_a = 1 - Ciplus1_a.array() ;Ciplus1_a = Ciplus1_a - lam(i.array())*m;
MatrixXd Ciplus1 = - Ciplus1_a;             //-(1-2*nu(j)*ONES-lam(i)*m)  
MatrixXd Diplus1_a =  2*nu(j.array())*ONES; Diplus1_a = 1 - Diplus1_a.array() ;Diplus1_a = Diplus1_a + lam(i.array())*m; //(1-2*nu(j)*ONES+lam(i)*m).*exp(-(lam(i)-lam(i-1))*m)         ...  % Di
MatrixXd Diplus1 = Diplus1_a.array() * Biplus1.array();

//Equations based on the vertical stresses
// Concatenate matrices horizontally
MatrixXd sigmaz(Ai.rows(), Ai.cols() + Bi.cols() + Ci.cols() + Di.cols() + Aiplus1.cols() + Biplus1.cols() + Ciplus1.cols() + Diplus1.cols());
sigmaz << Ai, Bi, Ci, Di, Aiplus1, Biplus1, Ciplus1, Diplus1;





i = i.array() - 1;
j = j.array() - 1;
// Define the equations based on shear stresses
MatrixXd Ai_tau = VectorXd::Ones(n-1)*ONES;
MatrixXd Bi_tau_temp = -(lam(i.array()) - lam(j.array())); Bi_tau_temp = Bi_tau_temp*m;
MatrixXd Bi_tau = Bi_tau_temp.array().exp();
MatrixXd Ci_tau_a =  2*nu(j.array())*ONES;
MatrixXd Ci_tau = Ci_tau_a + lam(i.array())*m; //-(1-2*nu(j)*ONES-lam(i)*m)  
MatrixXd Di_tau_a =  2*nu(j.array())*ONES; Di_tau_a = Di_tau_a - lam(i.array())*m; //(1-2*nu(j)*ONES+lam(i)*m).*exp(-(lam(i)-lam(i-1))*m)         ...  % Di
MatrixXd Di_tau = Di_tau_a.array() * Bi_tau.array(); 

i = i.array() + 1;
j = j.array() + 1;       

MatrixXd Ai_tauplus1 = VectorXd::Ones(n-1)*ONES;

MatrixXd Bi_tauplus1_temp = -(lam(i.array()) - lam(j.array())); Bi_tauplus1_temp = Bi_tauplus1_temp*m;

MatrixXd Bi_tauplus1 = Bi_tauplus1_temp.array().exp();

MatrixXd Ci_tauplus1_a =  2*nu(j.array())*ONES;

MatrixXd Ci_tauplus1 = Ci_tauplus1_a + lam(i.array())*m; //-(1-2*nu(j)*ONES-lam(i)*m)  

MatrixXd Di_tauplus1_a =  2*nu(j.array())*ONES; Di_tauplus1_a = Di_tauplus1_a - lam(i.array())*m; //(1-2*nu(j)*ONES+lam(i)*m).*exp(-(lam(i)-lam(i-1))*m)    
    
MatrixXd Di_tauplus1 = Di_tauplus1_a.array() * Bi_tauplus1.array(); 


// Combine the matrices
MatrixXd tauz(Ai_tau.rows(), Ai_tau.cols() + Bi_tau.cols() + Ci_tau.cols() + Di_tau.cols() + Ai_tauplus1.cols() + Bi_tauplus1.cols() + Ci_tauplus1.cols() + Di_tauplus1.cols());
tauz << Ai_tau, Bi_tau, Ci_tau, Di_tau, Ai_tauplus1, Bi_tauplus1, Ci_tauplus1, Di_tauplus1;


    // Print the sizes of Ai, Bi, Ci, and Di
    std::cout << "Size of Ai: " << Ai.rows() << " x " << Ai.cols() << std::endl;
    std::cout << "Size of Bi: " << Bi.rows() << " x " << Bi.cols() << std::endl;
    std::cout << "Size of Ci: " << Ci.rows() << " x " << Ci.cols() << std::endl;
    std::cout << "Size of Di: " << Di.rows() << " x " << Di.cols() << std::endl;
    std::cout << "Size of sigmaz: " << sigmaz.rows() << " x " << sigmaz.cols() << std::endl;
    

saveit(tauz);

 }




