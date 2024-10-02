#include "iitpave2.h"

void init_LET(Pave& iitpave) {

    int N, n ;
    double H, Laml;
    VectorXd E, nu, zi, a(iitpave.Xl.rows()); a = 150.8 * VectorXd::Ones(iitpave.Xl.rows());   iitpave.a = a;
    MatrixXd Xp, Xl;
   
 
   // INPUT PARAMETERS
    N = iitpave.N; //Number of Bessel roots in integration
    n = iitpave.n; //Number of integration points in between each Bessel root
    Xp = iitpave.Xp; //Output coordinates
    Xl = iitpave.Xl; //Load posistion coordinates
   // a = iitpave.a; //Load radii
    E = iitpave.E; //Layer Young's moduli
    nu = iitpave.nu; //Layer Poisson's ratios
    zi = iitpave.zi; zi.conservativeResize(iitpave.zi.size()-1); // Layer interface depths (i.e. depth of layer n-1 out
                     // of n layers (layer n depth goes to infinity)
 
// -------------------------------------------------------------------------
// References
// -------------------------------------------------------------------------
// [1] Burmister, D. M. (1945). The general theory of stresses and 
//     displacements in layered systems. I. Journal of applied physics, 16(2),
//     89–94.
// [2] Huang, Y.H., 2003. Pavement Analysis Design, 2nd Edition, 
//     Prentice-Hall, New Jersey.
// [3] Ioannides, A. M., & Khazanovich, L. (1998). General formulation for 
//     multilayered pavement systems. Journal of transportation engineering, 
//     124(1), 82–90.
// [4] Andersen, S., Levenberg, E., & Andersen, M. B (2020). Efficient 
//     reevaluation of surface displacements in a layered elastic half-space. 
//     The International Journal of Pavement Engineering 21(4), 1-8. 
//     https://doi.org/10.1080/10298436.2018.1483502
// -------------------------------------------------------------------------

// General parameters
H    = zi(zi.size()-1);       iitpave.H    = H;     // Bottom depth of last finite layer
Laml = zi(0)/H; iitpave.Laml = Laml;  // Relative height of top layer

// // Integration points and weights
iitpave.XipWip  = numint_coeff(N, n, Xp, Xl, H , a) ;  //*filename XipWip sIntegration points and weights
     
// Evaluate coefficients of integration 
iitpave.ABCD = arb_func(E.size(),zi,E,nu,iitpave) ;

// Evaluate response 
LET_response(iitpave);

//testing
//cout << iitpave.sigz << endl;


int size = iitpave.sigx.size();

MatrixXd Response_Sheet(iitpave.sigx.size(), 15);

Response_Sheet << iitpave.ux, iitpave.uy, iitpave.uz, iitpave.sigx, iitpave.sigy, iitpave.sigz,
                  iitpave.sigxy, iitpave.sigyz, iitpave.sigxz, iitpave.epsx, iitpave.epsy, 
                  iitpave.epsz, iitpave.epsxy, iitpave.epsyz, iitpave.epsxz;

cout << Response_Sheet << endl;

// Create and open a CSV file
std::ofstream file("response_table.csv");

// Write the header, including the new terms
file << "Index, ux, uy, uz, sigx, sigy, sigz, sigxy, sigyz, sigxz, epsx, epsy, epsz, epsxy, epsyz, epsxz" << std::endl;

// Write the data to the CSV file
for (int i = 0; i < iitpave.sigx.size(); ++i) {
    file << i + 1 << ", "
         << iitpave.ux(i) << ", "
         << iitpave.uy(i) << ", "
         << iitpave.uz(i) << ", "
         << iitpave.sigx(i) << ", "
         << iitpave.sigy(i) << ", "
         << iitpave.sigz(i) << ", "
         << iitpave.sigxy(i) << ", "
         << iitpave.sigyz(i) << ", "
         << iitpave.sigxz(i) << ", "
         << iitpave.epsx(i) << ", "
         << iitpave.epsy(i) << ", "
         << iitpave.epsz(i) << ", "
         << iitpave.epsxy(i) << ", "
         << iitpave.epsyz(i) << ", "
         << iitpave.epsxz(i) << std::endl;
}

// Close the file
file.close();


    std::cout << "Response table saved to 'response_table.csv'." << std::endl;

    
}