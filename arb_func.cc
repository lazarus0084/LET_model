#include "iitpave2.h"

void arb_func(int n, const VectorXd& zi, const VectorXd& E, const VectorXd& nu, const Pave& iitpave) {
    
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
//for i+1 & j+1 coefficients

RowVectorXi ip = i.array() + 1 ;

RowVectorXi jp = j.array() + 1 ;
//for i-1 & j-1 coefficients

RowVectorXi im = i.array() - 1 ;

RowVectorXi jm = j.array() - 1 ;


// Define the equations
MatrixXd Ai = VectorXd::Ones(n-1)*ONES;

MatrixXd Bi_temp = (lam(i.array()) - lam(j.array())); Bi_temp = -Bi_temp*m;
MatrixXd Bi = Bi_temp.array().exp();

MatrixXd Ci_a =  2*nu(j.array())*ONES; Ci_a = 1 - Ci_a.array() ;Ci_a = Ci_a - lam(i.array())*m;
MatrixXd Ci = - Ci_a; 
            //-(1-2*nu(j)*ONES-lam(i)*m)  
MatrixXd Di_a =  2*nu(j.array())*ONES; Di_a = 1 - Di_a.array() ;Di_a = Di_a + lam(i.array())*m; //(1-2*nu(j)*ONES+lam(i)*m).*exp(-(lam(i)-lam(i-1))*m)         ...  % Di
MatrixXd Di = Di_a.array() * Bi.array() ;

MatrixXd Ai2 = -(lam(ip.array()) - lam(i.array())) * m ; Ai2 = -Ai2.array().exp();//.........Ai+1
MatrixXd Bi2 = - Ai; //......................................................................Bi+1

MatrixXd Ci2a = -2*nu(jp.array())*ONES - lam(i.array()) * m ; Ci2a = 1 + Ci2a.array() ;
MatrixXd Ci2b = -(lam(ip.array()) - lam(i.array())) * m ; Ci2b = Ci2b.array().exp();//.........
MatrixXd Ci2 = Ci2a.array() * Ci2b.array() ;

MatrixXd Di2a =  2*nu(jp.array())*ONES ; MatrixXd Di2b = lam(i)*m ;
MatrixXd Di2 = -1 + Di2a.array()  - Di2b.array() ;//.............................................. % Di+1
                              
//Equations based on the vertical stresses
// Concatenate matrices horizontally
MatrixXd sigmaz(Ai.rows(), Ai.cols() + Bi.cols() + Ci.cols() + Di.cols() + Ai2.cols() + Bi2.cols() + Ci2.cols() + Di2.cols());
sigmaz << Ai, Bi, Ci, Di, Ai2, Bi2, Ci2, Di2;

// Define the equations based on shear stresses
MatrixXd Ai_tau = VectorXd::Ones(n-1)*ONES;

MatrixXd Bi_tau_a = -(lam(i.array()) - lam(im)) *m ; 
MatrixXd Bi_tau = -Bi_tau_a.array().exp() ;//-exp(-(lam(i)-lam(i-1))*m)  

MatrixXd Ci_tau = 2 * nu(j.array()) * ONES + lam(i.array()) * m ;  //(2*nu(j)*ONES+lam(i)*m)     

MatrixXd Di_tau_a = 2*nu(j.array())*ONES - lam(i.array()) * m ;
MatrixXd Di_tau_b = - (lam(i.array()) - lam(im.array())) * m ;Di_tau_b = Di_tau_b.array().exp() ;
MatrixXd Di_tau = Di_tau_a.array() * Di_tau_b.array() ;

MatrixXd Ai_tau2 = -(lam(ip.array()) - lam(i.array())) * m;
Ai_tau2 = -Ai_tau2.array().exp() ;  // -exp(-(lam(i+1)-lam(i))*m)                                     ... % A{i+1}

MatrixXd Bi_tau2 = Ai_tau ;


MatrixXd Ci_tau2_a = -(2*nu(jp.array())*ONES + lam(i.array()) * m ) ;
MatrixXd Ci_tau2_b = - (lam(ip.array()) - lam(i.array())) * m ;Ci_tau2_b = Ci_tau2_b.array().exp() ;
MatrixXd Ci_tau2 = Ci_tau2_a.array() * Ci_tau2_b.array() ;

MatrixXd Di_tau2 = - (2 * nu(jp.array()) * ONES - lam(i.array()) * m) ; 


// Combine the matrices
MatrixXd taurz(Ai_tau.rows(), Ai_tau.cols() + Bi_tau.cols() + Ci_tau.cols() + Di_tau.cols() + Ai_tau2.cols() + Bi_tau2.cols() + Ci_tau2.cols() + Di_tau2.cols());

taurz << Ai_tau, Bi_tau, Ci_tau, Di_tau, Ai_tau2, Bi_tau2, Ci_tau2, Di_tau2;   

MatrixXd taurz1(taurz.rows(),taurz.cols()) ;
MatrixXd taurz2(taurz.rows(),taurz.cols()) ;
//...................................................Frictionless...............................................................................
if (iitpave.bond == "Frictionless") {

Ai_tau = VectorXd::Ones(n-1)*ONES;

Bi_tau_a = -(lam(i.array()) - lam(im)) *m ; 
Bi_tau = -Bi_tau_a.array().exp() ;//-exp(-(lam(i)-lam(i-1))*m)  

Ci_tau = 2 * nu(j.array()) * ONES + lam(i.array()) * m ;  //(2*nu(j)*ONES+lam(i)*m)     

Di_tau_a = 2*nu(j.array())*ONES - lam(i.array()) * m ;
Di_tau_b = - (lam(i.array()) - lam(im.array())) * m ;Di_tau_b = Di_tau_b.array().exp() ;
Di_tau = Di_tau_a.array() * Di_tau_b.array() ;

Ai_tau2 = Ai_tau2.setZero() ;

Bi_tau2 = Bi_tau2.setZero() ;

Ci_tau2 = Ci_tau2.setZero() ;

Di_tau2 = Di_tau2.setZero() ;

// Combine the matrices
//MatrixXd taurz1(Ai_tau.rows(), Ai_tau.cols() + Bi_tau.cols() + Ci_tau.cols() + Di_tau.cols() + Ai_tau2.cols() + Bi_tau2.cols() + Ci_tau2.cols() + Di_tau2.cols());

taurz1 << Ai_tau, Bi_tau, Ci_tau, Di_tau, Ai_tau2, Bi_tau2, Ci_tau2, Di_tau2;  

Ai_tau = Ai_tau.setZero() ;
Bi_tau = Bi_tau.setZero() ;
Ci_tau = Ci_tau.setZero() ;
Di_tau = Di_tau.setZero() ;

Ai_tau2 = -(lam(ip.array()) - lam(i.array())) * m;
Ai_tau2 = -Ai_tau2.array().exp() ;            
Bi_tau2 = VectorXd::Ones(n-1)*ONES;
Ci_tau2_a = -(2*nu(jp.array())*ONES + lam(i.array()) * m ) ;
Ci_tau2_b = - (lam(ip.array()) - lam(i.array())) * m ;Ci_tau2_b = Ci_tau2_b.array().exp() ;
Ci_tau2 = Ci_tau2_a.array() * Ci_tau2_b.array() ;
Di_tau2 = - (2 * nu(jp.array()) * ONES - lam(i.array()) * m) ; 

//MatrixXd taurz2(Ai_tau.rows(), Ai_tau.cols() + Bi_tau.cols() + Ci_tau.cols() + Di_tau.cols() + Ai_tau2.cols() + Bi_tau2.cols() + Ci_tau2.cols() + Di_tau2.cols());

taurz2 << Ai_tau, Bi_tau, Ci_tau, Di_tau, Ai_tau2, Bi_tau2, Ci_tau2, Di_tau2; 
}
 //........................................................................................................................................................................

//***Equations based on the vertical displacements****
MatrixXd Ai_uz = 1 + nu(j.array()).array() ;
Ai_uz = Ai_uz.array() / E(j.array()).array() ; 
Ai_uz = Ai_uz * ONES ;

MatrixXd Bi_uz_a =  -Ai_uz ;
MatrixXd Bi_uz= -(lam(i.array()) - lam(j.array())) * m ; Bi_uz = Bi_uz.array().exp() ;
Bi_uz = Bi_uz_a.array() * Bi_uz.array() ;  

MatrixXd Ci_uz = -4*nu(j.array())*ONES ; Ci_uz = 2 + Ci_uz.array(); Ci_uz = Ci_uz - lam(i.array()) *m ;
Ci_uz = - Ai_uz.array() * Ci_uz.array() ;

MatrixXd Di_uz_a = -4*nu(j.array())*ONES ; Di_uz_a = 2 + Di_uz_a.array(); Di_uz_a = Di_uz_a + lam(i.array()) *m ;
Di_uz_a = - Ai_uz.array() * Di_uz_a.array() ; 
MatrixXd Di_uz = -(lam(i.array()) - lam(j.array())) * m ; Di_uz = Di_uz.array().exp() ;
Di_uz = Di_uz_a.array() * Di_uz.array() ;

MatrixXd Ai2_uz_a = 1 + nu(jp.array()).array() ;
Ai2_uz_a = - Ai2_uz_a.array() / E(jp.array()).array() ; Ai2_uz_a = Ai2_uz_a * ONES ;
MatrixXd Ai2_uz = -(lam(ip.array()) - lam(i.array())) * m ; Ai2_uz = Ai2_uz.array().exp() ;
Ai2_uz = Ai2_uz.array() * Ai2_uz_a.array() ;

MatrixXd Bi2_uz = 1 + nu(jp.array()).array() ;
Bi2_uz = Bi2_uz.array() / E(jp.array()).array() ; 
Bi2_uz = Bi2_uz * ONES ;                        //(1+nu(j+1))./E(j+1)*ONES........... % B{i+1}


MatrixXd Ci2_uz_a = -4 * nu(jp.array()) * ONES;
Ci2_uz_a = 2 + Ci2_uz_a.array();
Ci2_uz_a = Ci2_uz_a - lam(i.array()) * m;
Ci2_uz_a = Bi2_uz.array() * Ci2_uz_a.array();

MatrixXd Ci2_uz = -(lam(ip.array()) - lam(i.array())) * m;
Ci2_uz = Ci2_uz.array().exp();
Ci2_uz = Ci2_uz_a.array() * Ci2_uz.array();

MatrixXd Di2_uz_a = 2 - (4 * nu(jp.array()) * ONES).array() ;
Di2_uz_a = Di2_uz_a + lam(i.array()) * m ;
MatrixXd Di2_uz = Bi2_uz.array() * Di2_uz_a.array() ;

// Combine the matrices
MatrixXd uz(Ai_uz.rows(), Ai_uz.cols() + Bi_uz.cols() + Ci_uz.cols() + Di_uz.cols() + Ai2_uz.cols() + Bi2_uz.cols() + Ci2_uz.cols() + Di2_uz.cols());

uz << Ai_uz, Bi_uz, Ci_uz, Di_uz, Ai2_uz, Bi2_uz, Ci2_uz, Di2_uz; 

//Equations based on the radial displacements

MatrixXd Ai_ur = Ai_uz ;

MatrixXd Bi_ur = - Bi_uz ;

MatrixXd Ci_ur = lam(i.array()) * m; Ci_ur = 1 + Ci_ur.array() ;
Ci_ur = Ai_ur.array() * Ci_ur.array() ;


MatrixXd Di_ur_a = 1 + nu(j.array()).array() ;
Di_ur_a = Di_ur_a.array() / E(j.array()).array() ; 
Di_ur_a = - Di_ur_a * ONES ;    
MatrixXd Di_ur_b = lam(i.array()) * m; Di_ur_b = 1 - Di_ur_b.array() ; Di_ur_b = Di_ur_a.array() * Di_ur_b.array() ;
MatrixXd Di_ur = -(lam(i.array()) - lam(j.array()))*m ; Di_ur = Di_ur.array().exp() ;
Di_ur = Di_ur_b.array() * Di_ur.array() ;
MatrixXd Ai2_ur_a = 1 + nu(jp.array()).array() ;
Ai2_ur_a = - Ai2_ur_a.array() / E(jp.array()).array() ; 
Ai2_ur_a = Ai2_ur_a * ONES ;
MatrixXd Ai2_ur_b =  -(lam(ip.array()) - lam(i.array())) * m ; Ai2_ur_b = Ai2_ur_b.array().exp() ;
MatrixXd Ai2_ur =  Ai2_ur_a.array() * Ai2_ur_b.array() ;

MatrixXd Bi2_ur = 1 + nu(jp.array()).array();
Bi2_ur = Bi2_ur.array() / E(jp.array()).array();
Bi2_ur = - Bi2_ur * ONES ;

MatrixXd Ci2_ur_a = 1 + nu(jp.array()).array();
Ci2_ur_a = -Ci2_ur_a.array() / E(jp.array()).array();
Ci2_ur_a = Ci2_ur_a * ONES; Ci2_ur_a = Ci2_ur_a.array() * (1 + (lam(i.array()) * m).array()).array();
MatrixXd Ci2_ur_b =   -(lam(ip.array()) - lam(i.array())) * m ; Ci2_ur_b = Ci2_ur_b.array().exp() ;        //   lam(i.array()) * m;
MatrixXd Ci2_ur = Ci2_ur_a.array() * Ci2_ur_b.array();

MatrixXd Di2_ur = 1 + nu(jp.array()).array();
Di2_ur = Di2_ur.array() / E(jp.array()).array();
Di2_ur = Di2_ur * ONES; 
Di2_ur = Di2_ur.array() * (1 - (lam(i.array()) * m).array()).array();
MatrixXd ur(Ai_ur.rows(), Ai_ur.cols() + Bi_ur.cols() + Ci_ur.cols() + Di_ur.cols() + Ai2_ur.cols() + Bi2_ur.cols() + Ci2_ur.cols() + Di2_ur.cols());

ur << Ai_ur, Bi_ur, Ci_ur, Di_ur, Ai2_ur, Bi2_ur, Ci2_ur, Di2_ur;

//Arrange the set of linear matrix system to solve for the parameters Ai, 
// Bi, Ci and Di
//Surface boundary conditions


//MatrixXd term1 = -m * (lam(1) - lam(0)) ; MatrixXd term2 = 1*0NES ;
RowVectorXd term1 = (-m * (lam(1) - lam(0))).array().exp() ;
RowVectorXd term2 = 1* ONES;
RowVectorXd term3 = -(1 - 2 * nu(0)) * (-m *(lam(1) - lam(0))).array().exp() ;
RowVectorXd term4 = 1 - (2 * nu(0) * ONES).array() ;
SparseMatrix<double> term5_temp(1, 4*Lm); term5_temp.reserve(0) ;
MatrixXd term5_dense = MatrixXd(term5_temp) ;
RowVectorXd term6 = (2 * nu(0)) * (-m *(lam(1) - lam(0))).array().exp() ;
RowVectorXd term7 = 2 * nu(0) * ONES ;


 // Create temporary dense matrix for BC0
    int rows = 2 * term1.rows();
    int columns = term1.cols() + term2.cols() + term3.cols() + term4.cols();
    MatrixXd BC0_temp(2 * term1.rows(), term1.cols() + term2.cols() + term3.cols() + term4.cols());

    BC0_temp << term1, term2, term3, term4, term1, -term2, term6, term7;

    // Convert dense matrix BC0_temp to a sparse matrix
    SparseMatrix<double> BC0_temp_sparse = BC0_temp.sparseView();
    SparseMatrix<double> SparseAllocation = spalloc(2, 4 * Lm, 0);

    // Create a new sparse matrix BC0 for the concatenation
    SparseMatrix<double> BC0(BC0_temp_sparse.rows(), BC0_temp_sparse.cols() + SparseAllocation.cols());

    // Copy BC0_temp_sparse into the new sparse matrix BC0
    for (int k = 0; k < BC0_temp_sparse.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(BC0_temp_sparse, k); it; ++it) {
            BC0.insert(it.row(), it.col()) = it.value();
        }
    }

    // Copy SparseAllocation into the new sparse matrix BC0 (starting from the last column of BC0_temp_sparse)
    for (int k = 0; k < SparseAllocation.outerSize(); ++k) {
        for (SparseMatrix<double>::InnerIterator it(SparseAllocation, k); it; ++it) {
            BC0.insert(it.row(), it.col() + BC0_temp_sparse.cols()) = it.value();
        }
    }

// All conditions in the intermediate layers
MatrixXd BCs(sigmaz.rows() + taurz.rows() + uz.rows() + ur.rows(), taurz.cols()) ;
BCs <<sigmaz,
      taurz,
      uz,
      ur;  
 if (iitpave.bond == "Frictionless"){
    BCs <<sigmaz,
          uz,
          taurz1,
          taurz2; 

 }

// When inserting m in order to evaluate BC0 and BCs (=BC), the first 
// Lm = length(m) columns of BC0 will be exp(-m.*...) with different m 
// values. The following Lm columns will contain 1 and -1. We want to
// reorganize such that the first columns contains the equations
// in BC0 with m(1), the following columns correspond to the equations 
// in BC0 with the m(2) value etc., i.e., we want to go from
// BC0 = [exp(-m(1).*...)  exp(-m(2).*...)  exp(-m(3).*...) .. ]
// to
// BC0 = [exp(-m(1).*...)  1  -(1-2*nu(1))*exp(-m(1)...) .. ]

// The same reorganization we want to do with the BCs matrix
// Vectors used to change the the column sequence are denoted indx_c

// indx(:)' = [0 0 0 .. 1 1 1 .. 2 2 2 .. Lm-1 Lm-1 Lm-1]: 
// 8 = number of equations (=columns) in BCs for one m-value
    // Create an Eigen VectorXd with linearly spaced values from start to end

int start; int end; int step; int totalElements;

// RowVectorXi vec = RowVectorXi::LinSpaced(Lm , 0, Lm - 1 ) ;
// MatrixXi indx_c = LinSeq(0, Lm,(8-1)*Lm +1) ;

start  = 1 ; step = Lm ; end  = (8 - 1) * Lm + 1 ;

totalElements =  (end - start)/ step ; totalElements = totalElements + 1 ;

RowVectorXi vec1 = RowVectorXi::LinSpaced(totalElements, start, end ) ;

RowVectorXi vec2 = RowVectorXi::LinSpaced(Lm , 0, Lm - 1 ) ;
 
MatrixXi indx_c =  repmat(vec1, 1, Lm) + repmat(vec2, 1, 8) ; 


// Also the rows of BCs should be reorganized in order for the lambda, nu 
// and E-values to be in the correct sequence. The index vector for this is
// denoted indx_r and is defined as

    // Define the matrix
    Eigen::MatrixXi indx_r(4, n-1);

    indx_r  <<   Eigen::RowVectorXi::LinSpaced(n-1, 1, n-1),
                 Eigen::RowVectorXi::LinSpaced(n-1, n, 2*n-2),
                 Eigen::RowVectorXi::LinSpaced(n-1, 2*n-1, 3*n-3),
                 Eigen::RowVectorXi::LinSpaced(n-1, 3*n-2, 4*n-4);

//Adjusting indices to 0
    indx_c.array() -= 1;
    indx_r.array() -= 1;

// Use indx_c and indx_r to reorganice BCs and BC0 and order them in a
// united matrix BC. Only the rows of BCs are to be reorganized. BC0 is fine.
// Select columns using indices in one statement
// Eigen::MatrixXi selected_columns = BC0(Eigen::all, indx_c);


//Eigen::MatrixXd selected = BC0(Eigen::all, indx_c); //BC0(:,indx_c)


           // Example row indices
    Eigen::VectorXi row_indices(3);
    row_indices << 0, 2, 3;

    // Example column indices
    Eigen::VectorXi column_indices(3);
    column_indices << 1, 3, 4;

    // Select rows and columns in one line using permutation matrices
    //Eigen::MatrixXd selected = BC0(row_indices, Eigen::all)(Eigen::all, column_indices);
           

//  Now we have a matrix BC that contains all the coefficients for each
//  integration point. BC has 4*n-2 rows (corresponding to the number of
//  unknowns) and 4*n (the number of unknowns plus 2 - later we reduce this 
//  to 4*n-2) columns per integration point, i.e. 4*n*Lm columns in
//  total. Each (4*n-2) x 4*n - referred to as a submatrix, that only 
//  represents a single integration point. We want to organize all submatrices 
//  in a diagonal matrix of dimension Lm*(4*n-2) x 4*n*Lm, with each submatrix 
//  decoupled from the others. 

//  Below and indx vector that inserts the first submatrix in BC into the
//  right positions in BCg is organized. The way to use indx is in the 
//  following way: BCg(indx) = BC(:)

//  indx for a single submatrix - step 1 of 5

vec1 = Eigen::RowVectorXi::LinSpaced(4*n-2, 1, 4*n-2) ;

MatrixXi indx = repmat( vec1 , 8, 1).transpose() ; // Indices on all cells in a single matrix. 
                                                                              //We have (4*n-2)*8 coefficients in a single 
                                                                              //sub matrix to fill in into the global matrix

// Additional content to indx is given (step 2 of 5). Two vectors a and c are 
// defined
vec1 = RowVectorXi::LinSpaced(8, 0, 7) ;
MatrixXi a = repmat(vec1,4*n-2,1);

start  = 0 ; step = 4 ; end  = 4*(n-2);
totalElements =  (end - start)/ step ; totalElements = totalElements + 1 ;
vec1 = RowVectorXi::LinSpaced(totalElements, start, end ) ;

VectorXi temp_var =  repmat(vec1,4,1).reshaped(4*vec1.size(),1);
VectorXi zerosVector = VectorXi::Zero(2); // Vector of zeros

    // Combine these vectors into a single vector

VectorXi arg(2 + 4*vec1.size()) ; // Size of the resulting vector
arg << zerosVector, temp_var;
MatrixXi c = repmat(arg, 8, 1);

//  Add vectors a and c to indx (step 3 of 5) - the content of indx will then
//  represent the indeces for the first submatrix in the BCs matrix for a 
//  single submatrix (the first one only)

//indx = indx(:) + (a(:)+c)*(4*n-2)*Lm;

indx = indx.reshaped(indx.size(),1) + (a.reshaped(a.size(),1) + c) * (4*n-2) *Lm ;

// Now indx is expanded to consider all Lm submatrices. A vector b is
// organized (step 4 of 5)
vec1 = RowVectorXi::LinSpaced(Lm, 0, Lm-1) ;
MatrixXi b = repmat(vec1, indx.size(), 1) ;

indx = repmat(indx, Lm, 1) + b.reshaped(b.size(),1)*(4*n*Lm+1)*(4*n-2) ;
// Insert BC content into BCg in the right positions
// tic
// BCg(indx) = BC; % corresponds to BCg(indx) = BC(:)
// time_BCg=toc
// Indices to subindeces

// [I,J] = ind2sub([Lm*(4*n-2) 4*n*Lm],indx);
// BCg   = sparse(I,J,BC,Lm*(4*n-2),4*n*Lm);

}



