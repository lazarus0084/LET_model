#include "iitpave2.h"

Pave LET_response(Pave& iitpave) {
    int N, n ;
    double H, Laml;
    VectorXd E, nu, zi, q, a(iitpave.Xl.rows()); a = 150.8 * VectorXd::Ones(iitpave.Xl.rows());   iitpave.a = a;
    MatrixXd Xp_temp,Xp,Xp_1, Xl,Xl1, XipWip,ABCD;
   
 
   // INPUT PARAMETERS
    N = iitpave.N; //Number of Bessel roots in integration
    n = iitpave.n; //Number of integration points in between each Bessel root
    Xp_temp = iitpave.Xp; //Output coordinates
    Xl1 = iitpave.Xl; //Load posistion coordinates
   // a = iitpave.a; //Load radii
    E = iitpave.E; //Layer Young's moduli
    nu = iitpave.nu; //Layer Poisson's ratios
    zi = iitpave.zi; zi.conservativeResize(iitpave.zi.size()-1); // Layer interface depths (i.e. depth of layer n-1 out
    q = iitpave.q;                 // of n layers (layer n depth goes to infinity)
    XipWip = iitpave.XipWip; // Integration points and weights
    ABCD = iitpave.ABCD;
    H = iitpave.H;
// -------------------------------------------------------------------------
// References
// -------------------------------------------------------------------------
// [1] Burmister, D. M. (1945). The general theory of stresses and 
//     displacements in layered systems. i. Journal of applied physics, 16(2),
//     89–94.
// [2] Huang, Y.H., 2003. Pavement Analysis Design, 2nd Edition, 
//     Prentice-Hall, New Jersey.
// [3] Ioannides, A. M., & Khazanovich, L. (1998). General formulation for 
//     multilayered pavement systems. Journal of transportation engineering, 
//     124(1), 82–90.
// -------------------------------------------------------------------------
// If we only have one layer, we generate two layers

// If zi has more entries than length(E)-1, these are removed. In principle
// the last entry in zi is infinite. However, the infinite value is not
// included in this code, as we do not operate with infinite number. We
// handle this in an optional way.

// Check and adjust the size of zi
    if (zi.size() <= E.size()) {
        // Remove elements from index E.size()-1 to end
        if (zi.size() > 0) {
            zi.conservativeResize(E.size()-1);
        }
    }
// Check if zi is empty
    if (zi.size() == 0) {
        cout << "The code is arranged such that we need minimum two layers!!" << std::endl;
    }
// Check if the lengths of E and nu are equal
    if (E.size() != nu.size()) {
        cout << "length(E) != length(nu) !!!" << std::endl;
    }

// Number of loads (xl) and deformation points (xd)
    int xl = Xl1.rows(); // Number of load points
    int xp = Xp_temp.rows(); // Number of deformation points
// Check if the length of q and a are correct
     if (q.size() < xl || a.size() < xl) {
      cout << "Number of q and/or a values is lower than the number of load coordinates (= number of loads), so all loads are given the same q- and a-values" << endl;
        q = VectorXd::Constant(xl, q(0));
        a = VectorXd::Constant(xl, a(0));
    }
// Increase Xl, q, alpha, Xd and z in order to estimate the radius from each
// of the loads to all of the points. The number of radii are: xl*xd
    Xl = repmat(Xl1, xp, 1);   // Extends with the number of evaluation
                             // points
    q  = repmat(q, xp, 1);    //-||-
    MatrixXd alpha = repmat(a, xp, 1).array()/H;  // -||-
    Xp_1   = repmat(Xp_temp, xl, 1);   // -||- <-- is reorganized below to 
                             // correspond to the right Xl content/sequence
// Reorganize Xd so the first xl rows correspond to the same deformation
// point, and following xl rows correspond to the next deformation point
// etc... see description of the vector r nine lines below to understand the 
// organization

// Built dx vector used to reorganize/sort Xd in the right order
    RowVectorXi Sequence  = RowVectorXi::LinSpaced(xp, 1, xp);
    MatrixXi dx1 = repmat(Sequence, xl, 1); // Assuming repmat is your own function

   // Now clear and reuse arg1 with new size
    Sequence.resize(xl);  // Resize to hold 'xl' elements

    Sequence = RowVectorXi::LinSpaced(xl, 0, xp * (xl - 1));

    RowVectorXi dx = dx1.reshaped(1,dx1.size()) + repmat(Sequence,1,xp); dx = dx.array() -1;

     Xp = Xp_1(dx , all);
     
     VectorXd z = Xp(Eigen::all, 2);
    

// // Built dx vector used to reorganize/sort Xd in the right order
// VectorXi arg1 = Eigen::RowVectorXi::LinSpaced(xp, 1, xp);
// VectorXi dx = repmat(arg1,xl,1);
// arg1 = Eigen::VectorXi::LinSpaced(xl, 0, xp * (xl - 1));
// dx = dx + repmat(arg1,xp,1); 
// dx = dx.array() -1; // adjusting to 0 based indexing
// Xp = Xp(dx,Eigen::all);  //Eigen::MatrixXd selected_rows = Xd(row_indices, Eigen::all);
// 
// VectorXd z  = Xp(Eigen::all,2);
// cout << dx << endl;

// Evaluate the radii between evaluation point and load
// r = [  evaluation point 1 and load 1
//        evaluation point 1 and load 2
//        evaluation point 1 and load 3
//                        .
//        evaluation point 1 and load N
//        evaluation point 2 and load 1
//        evaluation point 2 and load 2
//                        .
//                        .
//      last evaluation point and load N];


   VectorXd r = ((Xp.col(0) - Xl.col(0)).array().square() + (Xp.col(1) - Xl.col(1)).array().square()).sqrt();


// Introduce the parameter rho (notice: since r is a vector, rho is too)
   VectorXd rho = r.array()/H;
   rho = (rho.array() == 0).select(1e-20, rho); // fu a small value to avoid singularity in 

                   // calculation of stresses
// Integration points and weights (takes out the length(rho) rows in XipWip
// related to one load for one position)
// int lrho = rho.size();                % Number of load-to-displacement-points
// MatrixXd xip_z = XipWip(0*lrho+1:1*lrho,Eigen::all); % takes out the first length(rho) rows
// wip_z = XipWip(1*lrho+1:2*lrho,:); % takes out next
// MatrixXd xip_r = XipWip(2*lrho+1:3*lrho,:); % etc...
// wip_r = XipWip(3*lrho+1:4*lrho,:); % etc...



    int lrho = rho.size() ;  // Number of load-to-displacement points

// Take out the first lrho rows
//arg1 = VectorXi::LinSpaced(lrho, 0*lrho + 1, 1*lrho).array() -1;
MatrixXd xip_z = XipWip.block(0 * lrho, 0, lrho, XipWip.cols()); //xip_z = XipWip(0*lrho+1:1*lrho,:);

// Take out the next lrho rows
MatrixXd wip_z = XipWip.block(1 * lrho, 0, lrho, XipWip.cols());

// Take out the next lrho rows
MatrixXd xip_r = XipWip.block(2 * lrho, 0, lrho, XipWip.cols());

// Take out the next lrho rows
MatrixXd wip_r = XipWip.block(3 * lrho, 0, lrho, XipWip.cols()); //xip_r = XipWip(2*lrho+1:3*lrho,:); % etc...
// Number of columns in reorganized integration points and weights (= number
// of integration points per integral)
int nz = xip_z.cols(); 
// The depth z defines which E, nu, Ai, Bi, Ci and Di values should be used
// in the different required evaluations of the response. A
// horizontal vector called layer_no gives the information of which
// values should be used. The length(layer_no) corresponds to the number of
// radii, r, we have (i.e. the length of rho = r/H).

// Layer number where the point we use is given*
//VectorXd layer_no = sum(repmat(z,1,length(zi)) > repmat(zi,length(rho),1),2)'+1; 

    // Perform element-wise comparison
//MatrixXd comparison_result = (repmat(z, 1, z.size()).array() > repmat(zi, rho.size(), 1).array()).cast<double>();

RowVectorXd zi_reshaped = zi.transpose();

MatrixXd matrix1 = repmat(z, 1, zi.size()); MatrixXd matrix2 =  repmat(zi_reshaped, rho.size(), 1);

// Perform element-wise comparison
MatrixXd comparison_result = (matrix1.array() > matrix2.array()).cast<double>(); // Cast boolean result to double

RowVectorXd layer_no = comparison_result.rowwise().sum(); layer_no = layer_no.array() + 1; 
// *Note: if (Xd,z) is larger than (zi,1) the function gives 1+1
// = layer no. 2, if (Xd,z) is larger than (zi,2) the function gives 1+1+1 
// = layer no. 3 etc.

VectorXd Lami(rho.size()); Lami.setZero();
VectorXd Lami1 = Lami;


   for (int  i = 0; i < layer_no.size(); ++i) {
        // Adjust for 0-based indexing in C++ (layer_no[i] is assumed to be 1-based like in MATLAB)
        int layer_idx = layer_no[i] - 1;

        // First condition
        if (layer_idx >= zi.size()) {
            Lami[i] = zi(zi.size() - 1) / H * 1e20; // Multiply by a large number (like 1e20)
    
        } else {
            Lami[i] = zi(layer_idx) / H;

        }

        // Second condition for Lami1
        if (layer_idx == 0) { // (= lam_{i-1} = lam{0})
            Lami1[i] = 0;
  
        } else {
            Lami1[i] = zi(layer_idx - 1) / H;
      
        }
    }

  // Evaluating response

// Summation form
// length(rho) = number of point to load radii, nz = number of integration 
// points in every point to load integral
  
  MatrixXd Aa(rho.size(),nz); Aa.setZero();

  MatrixXd Bb, Cc, Dd;

  Bb = Aa; Cc = Aa; Dd = Aa;

  MatrixXd ABCDz(E.size(), xip_z.cols());

   for (int  j = 0; j < rho.size(); ++j) {
   
    ABCDz = arb_func_interp(E.size(), xip_z.row(j), ABCD);

    Aa.row(j) = ABCDz.row(4*(layer_no(j)-1) + 0);  Bb.row(j) = ABCDz.row(4*(layer_no(j)-1) + 1);

    Cc.row(j) = ABCDz.row(4*(layer_no(j)-1) + 2);  Dd.row(j) = ABCDz.row(4*(layer_no(j)-1) + 3);
            
}

VectorXd Nu, Ee ;

Nu = nu(layer_no.array()-1);
Ee = E(layer_no.array()-1);

// Define exponential functions for one-step Richardson extrapolation
double x1 = pow(2, -20);
MatrixXd Iz1 = (-x1 * xip_z.array().square()).exp();
MatrixXd Iz2 = (-(x1 / 2.0) * xip_z.array().square()).exp();
MatrixXd Ir1 = (-x1 * xip_r.array().square()).exp();
MatrixXd Ir2 = (-(x1 / 2.0) * xip_r.array().square()).exp();

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//                                                                                                uz
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



VectorXd zbyH =  z/H ;
VectorXd LamizbyH =  Lami - z/H ;
VectorXd Lami1zbyH = z/H - Lami1;

MatrixXd uzs = besselj(0, (xip_z.array() * repmat(rho, 1, nz).array())).array() * 

(((Aa.array() - (Cc.array() * ((2 - 4*repmat(Nu,1,nz).array()) 
- (xip_z.array() * repmat(zbyH,1,nz).array())).array()).array()).array()
* ((-xip_z.array()*repmat(LamizbyH,1,nz).array()).array()).array().exp() ) -
 
 ((Bb.array() + (Dd.array() * ((2 - 4*repmat(Nu,1,nz).array()) 
+ (xip_z.array() * repmat(zbyH,1,nz).array())).array()).array()).array()
* ((-xip_z.array()*repmat(Lami1zbyH,1,nz).array()).array()).array().exp() )).array();

// Improve the convergence for points residing close to the surface

MatrixXd Iuz1 =  (-H * q.array() * alpha.array() * (1 + Nu.array()).array()).array()/ (Ee.array()).array()
                * ( (uzs.array() * ((((besselj(1,xip_z.array()*repmat(alpha,1,nz).array())).array()
                / (xip_z.array())).array()).array() *  wip_z.array() 
                * Iz1.array()).array()).rowwise().sum() ).array();

MatrixXd Iuz2 =  (-H * q.array() * alpha.array() * (1 + Nu.array()).array()).array()/ (Ee.array()).array()
                * ( (uzs.array() * ((((besselj(1,xip_z.array()*repmat(alpha,1,nz).array())).array()
                / (xip_z.array())).array()).array() *  wip_z.array() 
                * Iz2.array()).array()).rowwise().sum() ).array();

MatrixXd uz = (4*Iuz2.array() - Iuz1.array()).array()/ 3 ;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//                                                                                                ur
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



MatrixXd urs = besselj(1, xip_r.array()*repmat(rho,1,nz).array()).array()            *

     ((Aa.array() + (Cc.array() * (1 + (xip_r.array() * repmat(zbyH, 1, nz).array()).array()).array()).array()).array() * (((-xip_r.array()*repmat(LamizbyH,1,nz).array()).array()).array().exp() ).array()  
                                                                                         + 
     (Bb.array() - (Dd.array() * (1 - (xip_r.array() * repmat(zbyH, 1, nz).array()).array()).array()).array()).array() *  ((-xip_r.array() * repmat(Lami1zbyH, 1, nz).array()).array().exp() ).array() ).array();


MatrixXd Iur1 = ((H* q.array() * alpha.array() * (1 + Nu.array())).array()/Ee.array()).array() * ( (urs.array() * ((((besselj(1,xip_r.array()*repmat(alpha,1,nz).array())).array()
                / (xip_r.array())).array()).array() *  wip_r.array() 
                * Ir1.array()).array()).rowwise().sum() ).array(); 



MatrixXd Iur2 = ((H* q.array() * alpha.array() * (1 + Nu.array())).array()/Ee.array()).array() * ( (urs.array() * ((((besselj(1,xip_r.array()*repmat(alpha,1,nz).array())).array()
                / (xip_r.array())).array()).array() *  wip_r.array() 
                * Ir2.array()).array()).rowwise().sum() ).array(); 

MatrixXd ur = (4*Iur2 - Iur1).array()/3 ;



// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//                                                                                     STRESSES 
// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

                                                                                      //sigma_z


MatrixXd sigma_zs =   (( -xip_z.array() * besselj(0, xip_z.array()*repmat(rho,1,nz).array()).array()).array()   *

((Aa.array() -  (  Cc.array() * ( (1 - 2*repmat(Nu,1, nz).array()) - (xip_z.array() * repmat(zbyH, 1, nz).array()) ).array()).array()).array() *  (((-xip_z.array()*repmat(LamizbyH,1,nz).array()).array()).array().exp() ).array()

                                                                                                                  +
(Bb.array() +  (  Dd.array() * ( (1 - 2*repmat(Nu,1, nz).array()) + (xip_z.array() * repmat(zbyH, 1, nz).array()) ).array()).array()).array() * ((-xip_z.array() * repmat(Lami1zbyH, 1, nz).array()).array().exp()).array()).array()).array();


// Improve the convergence for points residing close to the surface

MatrixXd Isz1 = q.array() * alpha.array() * (((sigma_zs.array() * besselj(1, xip_z.array() * repmat(alpha, 1, nz).array()).array()).array()/xip_z.array()).array()*wip_z.array()*Iz1.array()).rowwise().sum() ;

MatrixXd Isz2 = q.array() * alpha.array() * (((sigma_zs.array() * besselj(1, xip_z.array() * repmat(alpha, 1, nz).array()).array()).array()/xip_z.array()).array()*wip_z.array()*Iz2.array()).rowwise().sum() ;

MatrixXd sigma_z = (4*Isz2 - Isz1).array()/3 ;

//                                                                                         sigma_r


MatrixXd sigma_rs = ( ((xip_r.array() * besselj(0, xip_r.array() * repmat(rho,1, nz).array()).array()) 

    - besselj(1,xip_r.array() * repmat(rho, 1, nz).array()).array()/repmat(rho, 1, nz).array()).array() *

     ((Aa.array() + (  Cc.array() * ( 1 + (xip_r.array() * repmat(zbyH, 1, nz).array()).array() ).array()).array()).array() *  (((-xip_r.array()*repmat(LamizbyH,1,nz).array()).array()).array().exp() ).array()
   
                                                                                   +    
    (Bb.array() - (  Dd.array() * ( 1 - (xip_r.array() * repmat(zbyH, 1, nz).array()).array() ).array()).array()).array() *  (((-xip_r.array()*repmat(Lami1zbyH,1,nz).array()).array()).array().exp() ).array()).array()).array()
                                                                                   + 
    ( (2*repmat(Nu, 1, nz).array() * xip_r.array() * besselj(0, xip_r.array() * repmat(rho, 1, nz).array()).array()).array()  *

(Cc.array() * ((-xip_r.array()*repmat(LamizbyH, 1, nz).array()).array().exp()).array() - Dd.array() * ((-xip_r.array()*repmat(Lami1zbyH, 1, nz).array()).array().exp()).array()).array() ).array();



//      Improve the convergence for points residing close to the surface

MatrixXd Isr1 = q.array() * alpha.array() * (((sigma_rs.array() * besselj(1, xip_r.array() * repmat(alpha, 1, nz).array()).array()).array()/xip_r.array()).array()*wip_r.array()*Iz1.array()).rowwise().sum() ;

MatrixXd Isr2 = q.array() * alpha.array() * (((sigma_rs.array() * besselj(1, xip_r.array() * repmat(alpha, 1, nz).array()).array()).array()/xip_r.array()).array()*wip_r.array()*Iz2.array()).rowwise().sum() ;

MatrixXd sigma_r = (4*Isr2 - Isr1).array()/3 ;






MatrixXd sigma_thetas = ( (besselj(1, xip_r.array() * repmat(rho, 1, nz).array()).array()/repmat(rho, 1, nz).array()).array() *

     ((Aa.array() + (  Cc.array() * ( 1 + (xip_r.array() * repmat(zbyH, 1, nz).array()).array() ).array()).array()).array() *  (((-xip_r.array()*repmat(LamizbyH,1,nz).array()).array()).array().exp() ).array()
   
                                                                                   +    
    (Bb.array() - (  Dd.array() * ( 1 - (xip_r.array() * repmat(zbyH, 1, nz).array()).array() ).array()).array()).array() *  (((-xip_r.array()*repmat(Lami1zbyH,1,nz).array()).array()).array().exp() ).array()).array()).array()
                                                                                   + 
    ( (2*repmat(Nu, 1, nz).array() * xip_r.array() * besselj(0, xip_r.array() * repmat(rho, 1, nz).array()).array()).array()  *

(Cc.array() * ((-xip_r.array()*repmat(LamizbyH, 1, nz).array()).array().exp()).array() - Dd.array() * ((-xip_r.array()*repmat(Lami1zbyH, 1, nz).array()).array().exp()).array()).array() ).array();


//      Improve the convergence for points residing close to the surface

 Isr1 = q.array() * alpha.array() * (((sigma_thetas.array() * besselj(1, xip_r.array() * repmat(alpha, 1, nz).array()).array()).array()/xip_r.array()).array()*wip_r.array()*Iz1.array()).rowwise().sum() ;

 Isr2 = q.array() * alpha.array() * (((sigma_thetas.array() * besselj(1, xip_r.array() * repmat(alpha, 1, nz).array()).array()).array()/xip_r.array()).array()*wip_r.array()*Iz2.array()).rowwise().sum() ;

MatrixXd sigma_theta = (4*Isr2 - Isr1).array()/3 ;



//%%%%% tau_rz %%%%%



MatrixXd tau_rzs =  (xip_z.array() * besselj(1, xip_z.array()*repmat(rho,1,nz).array()).array()).array() *

(((Aa.array() + (Cc.array() * ((2*repmat(Nu,1,nz)).array() + (xip_z.array() * repmat(zbyH,1 ,nz).array()).array()).array()).array()).array()*
((-xip_z.array() * repmat(LamizbyH, 1, nz).array() ).array().exp()).array()).array()
                                                                       -  
((Bb.array() - (Dd.array() * ((2*repmat(Nu,1,nz)).array() - (xip_z.array() * repmat(zbyH,1 ,nz).array()).array()).array()).array()).array()*
((-xip_z.array() * repmat(Lami1zbyH, 1, nz).array() ).array().exp()).array()).array() ).array() ;

//Improve the convergence for points residing close to the surface

MatrixXd Itr1 = q.array() * alpha.array() * (((tau_rzs.array() * besselj(1, xip_z.array() * repmat(alpha, 1, nz).array()).array()).array()/xip_z.array()).array()*wip_z.array()*Iz1.array()).rowwise().sum() ;

MatrixXd Itr2 = q.array() * alpha.array() * (((tau_rzs.array() * besselj(1, xip_z.array() * repmat(alpha, 1, nz).array()).array()).array()/xip_z.array()).array()*wip_z.array()*Iz2.array()).rowwise().sum() ;

MatrixXd tau_rz = (4*Itr2 - Itr1).array()/3 ;


// -------------------------------------------------------------------------
// Transformation between (r, theta) and (x, y) coordinates
// -------------------------------------------------------------------------
// Below is a transformation of each response (evaluated above from
// the (r, theta) coordinates to the (x, y) coordinates). After the
// transformation, the responses are added together.
//
//                         p  (load point)
//   (z is downwards)
//    z -> x              ^ (deformation vector ur; pointing towards the load)
//    |                 /
//  y v                /
//                    o  (deformation point)
//
// Note: the r-axis goes in the direction


 //SparseMatrix<double> T = spalloc(3*r.size(), 3*r.size(), (3*3 -4)*r.size());

// MatrixXd T(3 * r.size(), 3 * r.size() ); T.setZero(); //T = sparse(3*length(r),3*length(r)


//Organize content of the T matrix

 VectorXd COS = (Xp.col(0) - Xl.col(0)).array()/r.array() ;

 VectorXd SIN = (Xp.col(1) - Xl.col(1)).array()/r.array() ;

// Assuming r, COS, and SIN are Eigen matrices or vectors
// Define a small tolerance for floating-point comparison
 double epsilon = 1e-12; // Adjust the tolerance as needed

// Use the tolerance to check if the values in r are "close" to 0
   COS = (r.array().abs() < epsilon).select(1, COS);

   SIN = (r.array().abs() < epsilon).select(0, SIN);

//  % Address the content into the matrix

    MatrixXd T(3 * r.size(), 3 * r.size()); T.setZero(); T(3 * r.size() - 1, 3 * r.size() - 1) = (3*3-4) * r.size();
    VectorXd T_temp(3*r.size() * 3*r.size()); T_temp.setZero();

    int start = 1;  // Start at 1
    int step = 3*3*r.size() + 3;  // Step size is 147
    int end = T_temp.size() ;  // End size (2304)

    // Calculate how many steps are needed
    int numElements = (T_temp.size() - 3 * 3 * r.size()) / (3*3*r.size()) + 1;  // Number of elements in the sequence

    // Generate the sequence using LinSpaced
    VectorXi indices = VectorXi::LinSpaced(numElements, start, start + (numElements - 1) * step).array() -1;
    
   // indices = indices.array() -1;
   T_temp( VectorXi::LinSpaced(numElements, start, start + (numElements - 1) * step).array() -1) =  COS ;
   start = 2;
   T_temp( VectorXi::LinSpaced(numElements, start, start + (numElements - 1) * step).array() -1) = -SIN ;
   
   start = 3 * r.size() + 1;
   T_temp( VectorXi::LinSpaced(numElements, start, start + (numElements - 1) * step).array() -1) =  SIN ;
   
   start = 3 * r.size() + 2;
   T_temp( VectorXi::LinSpaced(numElements, start, start + (numElements - 1) * step).array() -1) =  COS ;
   
   start = 2 * 3*  r.size() + 3;
   T_temp( VectorXi::LinSpaced(numElements, start, start + (numElements - 1) * step).array() -1) = VectorXd::Ones(COS.size()); 

   T = T_temp.reshaped(3 * r.size(), 3 * r.size());

// Organize displacement vector
    VectorXd u1(3 * r.size(), 1); u1.setZero();
    start = 1;  // Start at 1
    step = 3;  // Step size is 147
    end = 3 * r.size();  // End size (2304)

        // Generate the sequence using LinSpaced
   // VectorXi indices = VectorXi::LinSpaced( (end - start) / step + 1, start, start + ((end - start) / step) * step);

    u1(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1) = ur;

    start = 3;
    //u1(2:3:end) = 0; % Displacements transverse to the radius (u_theta = 0)
    u1(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1) = uz;
    
    // Organize stress matrix

    MatrixXd sig1(3 * r.size(), 3 * r.size()); sig1.setZero(); sig1(3 * r.size() - 1, 3 * r.size() - 1) = (3*3-4) * r.size();
    VectorXd sig1_temp(3*r.size() * 3*r.size()); sig1_temp.setZero();
    

    // Flattening r and theta responses
    VectorXd Flat_sigma_r = sigma_r.reshaped(sigma_r.size(), 1);
    VectorXd Flat_sigma_z = sigma_z.reshaped(sigma_z.size(), 1);
    VectorXd Flat_sigma_theta = sigma_theta.reshaped(sigma_theta.size(), 1);
    VectorXd Flat_tau_rz = tau_rz.reshaped(tau_rz.size(), 1);
                
    start = 1;  // Start at 1
    step = 3 * 3* r.size() + 3;  // Step size is 147
    end =  sig1_temp.size() ;

    // Generate the sequence using LinSpaced
    sig1_temp(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1) 
    = Flat_sigma_r.head((VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1).size());

    start = 3*r.size() + 2 ;  // Start at 1
    
    sig1_temp(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1) 
    = Flat_sigma_theta.head((VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1).size());

    start =2*3*r.size() + 3 ;  // Start at 1

    sig1_temp(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1) 
    = Flat_sigma_z.head((VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1).size());

    start = 3 ;  /// Start at 1

    sig1_temp(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1) 
    = Flat_tau_rz.head((VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1).size());

    
    start =2*3*r.size() + 1 ;  // Start at 1

    sig1_temp(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1) 
    = Flat_tau_rz.head((VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1).size());

    sig1 = sig1_temp.reshaped(3 * r.size(), 3 * r.size());

    // Transform deformations in to (x,y,z)-coordinates
    u1 = T.transpose() * u1;
    MatrixXd u  = u1;
    
    // Transform stresses into (x,y,z)-coordinates
    sig1 = T.transpose() * sig1 * T;
    MatrixXd sigm  = sig1;

    // Transform 3x3 matrices into 6x1 vectors
    VectorXi dos(6*xp*xl); dos.setZero();

    start = 1;  
    step = 9*xp*xl + 3;
    end =  9*xp*xl*r.size();

    VectorXi dos_pop1 = VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step) ;

    
    start = 1;  
    step = 6;
    end =  dos.size() - 5;

    dos(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1) 
    = dos_pop1.head((end - start) / step + 1);




    start = 2;  
    step = 6;
    end =  dos.size() - 4;
    
    
    dos(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1) 
    = dos(VectorXi::LinSpaced((dos.size() - 1) / 6 + 1, 1, 1 + ((dos.size() - 1) / 6) * 6).array() -1).array() + 3*xp*xl+1;

    start = 3;  
    step = 6;
    end =  dos.size() - 3;
    
    dos(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1) 
    = dos(VectorXi::LinSpaced((dos.size() - 6) / 6 + 1, 1, 1 + ((dos.size() - 6) / 6) * 6).array()-1).array() + 2*3*xp*xl+2;
    
    start = 4;  
    step = 6;
    end =  dos.size() - 2;
    
    dos(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1) 
    =  dos(VectorXi::LinSpaced((dos.size() - 6) / 6 + 1, 1, 1 + ((dos.size() - 6) / 6) * 6).array()-1).array() + 1;
      
    start = 5;  
    step = 6;
    end =  dos.size() - 1;
    
    dos(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1) 
    =  dos(VectorXi::LinSpaced((dos.size() - 6 ) / 6 + 1, 2, 2 + ((dos.size() - 6) / 6) * 6).array()-1).array() + 1;

    start = 6;  
    step = 6;
    end =  dos.size();
    
    dos(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() - 1) 
    =  dos(VectorXi::LinSpaced((dos.size() - 6) / 6 + 1, 1, 1 + ((dos.size() - 6) / 6) * 6).array()-1).array() + 2;

    dos = dos.array() -1;
    
    VectorXd sigm_flat = sigm.reshaped() ;
    VectorXd sigv = sigm_flat(dos);

  // Organize matrix Txyz for adding displacements together
    
    VectorXi unos(3*xp*xl); unos.setZero();

    start = 1;  step = 9 * xp; end = 9*xl*xp ;
    
    unos(VectorXi::LinSpaced(xl,1,xl).array() -1 ) =  VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step);

    unos(VectorXi::LinSpaced(xl, xl + 1, 2 * xl).array() -1 ) =   unos(VectorXi::LinSpaced(xl,1,xl).array() -1 ).array() + 3*xp + 1;  

    unos(VectorXi::LinSpaced(xl, 2 * xl + 1, 3 * xl).array() -1 ) =   unos(VectorXi::LinSpaced(xl,1,xl).array() -1 ).array() + 2 *3 * xp + 2;  


   for (int i = 1; i < xp ; ++i) {

   VectorXi unos_rhs = unos(VectorXi::LinSpaced(3 * xl, 1, 3 * xl).array() -1).array() + i *(3 * xp * 3 * xl + 3);

    unos(VectorXi::LinSpaced(3 * xl, 3 * xl + 1 + 3 * xl * (i - 1), 3 * xl + 3 * xl * i).array() -1)
    = unos_rhs ;
  }


    // Get row and column indices from linear indices
    MatrixXi sub_indices = ind2sub(3*xp,3*r.size(), unos);
    sub_indices = sub_indices.array() -1;

     // Define the sparse matrix (5x5 example)
    SparseMatrix<double> Txyz = spalloc(3 * xp, 3 * r.size(), unos.size()) ;

      // List of triplets to hold the non-zero elements
    vector<Triplet<double>> tripletList;

    // Populate the triplet list with values
    for (int i = 0; i < sub_indices.rows() ; ++i) {
        tripletList.emplace_back(sub_indices(i,0), sub_indices(i,1), 1.0);  // Insert element 1 at (row_ind, col_ind)
    }


    // Populate the sparse matrix with the triplet list
    Txyz.setFromTriplets(tripletList.begin(), tripletList.end());

 VectorXi tres(6*xp*xl); tres.setZero();
 start = 1;  step = 9 * 4 * xp; end = 9 * 4 * xl * xp ;
    
 tres(VectorXi::LinSpaced(xl,1,xl).array() -1) =  VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step);   

    
 tres(VectorXi::LinSpaced(xl,xl+1,2*xl).array() -1) =  tres(VectorXi::LinSpaced(xl,1,xl).array() -1).array() + 6 * xp + 1;
 
    
 tres(VectorXi::LinSpaced(xl,2*xl+1,3*xl).array() -1) =  tres(VectorXi::LinSpaced(xl,1,xl).array() -1).array() + 2 * 6 * xp + 2;

 tres(VectorXi::LinSpaced(xl,3*xl+1,4*xl).array() -1) =  tres(VectorXi::LinSpaced(xl,1,xl).array() -1).array() + 3 * 6 * xp + 3;

 tres(VectorXi::LinSpaced(xl,4*xl+1,5*xl).array() -1) =  tres(VectorXi::LinSpaced(xl,1,xl).array() -1).array() + 4 * 6 * xp + 4;

 tres(VectorXi::LinSpaced(xl,5*xl+1,6*xl).array() -1) =  tres(VectorXi::LinSpaced(xl,1,xl).array() -1).array() + 5 * 6 * xp + 5;

  VectorXi tres_rhs(6 * xl);

   for (int i = 1; i < xp ; ++i) {

   tres_rhs = tres(VectorXi::LinSpaced(6 * xl, 1, 6 * xl).array() -1).array() + i *(6 * xp * 6 * xl + 6);

    tres(VectorXi::LinSpaced(6 * xl, 6*xl+1+6*xl*(i-1), 6*xl+6*xl*i).array() -1) = tres_rhs ;

  }

    // Get row and column indices from linear indices
    MatrixXi sub_indices2 = ind2sub(6*xp, 6*r.size(), tres);
    sub_indices2 = sub_indices2.array() - 1;

     // Define the sparse matrix (5x5 example)
    SparseMatrix<double> Mxyz = spalloc(6 * xp, 6 * r.size(), tres.size()) ;

      // List of triplets to hold the non-zero elements
    vector<Triplet<double>> tripletList2;

    // Populate the triplet list with values
    for (int i = 0; i < sub_indices2.rows() ; ++i) {
        tripletList2.emplace_back(sub_indices2(i,0), sub_indices2(i,1), 1.0);  // Insert element 1 at (row_ind, col_ind)
    }


    // Populate the sparse matrix with the triplet list
    Mxyz.setFromTriplets(tripletList2.begin(), tripletList2.end());

// Final response vector

 // Displacements
u          = Txyz*u;

start = 1;  step = 3; end = u.size() -2 ;
    
iitpave.ux  =  u(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() -1, Eigen::all);   

start = 2;  step = 3; end = u.size() -1 ;
    
iitpave.uy  =  u(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() -1, Eigen::all);   

start = 3;  step = 3; end = u.size() ;
    
iitpave.uz =  u(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() -1, Eigen::all);   

// Stresses

MatrixXd sig = -Mxyz*sigv ;     // minus inserted to fullfil inwards positive


start = 1;  step = 6; end = sig.size() -5 ;
    
iitpave.sigx  =  sig(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() -1, Eigen::all);   

start = 2;  step = 6; end = sig.size() -4 ;
    
iitpave.sigy  =  sig(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() -1, Eigen::all);   

start = 3;  step = 6; end = sig.size() -3 ;
    
iitpave.sigz =  sig(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() -1, Eigen::all);   

start = 4;  step = 6; end = sig.size() -2 ;
    
iitpave.sigxy =  sig(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() -1, Eigen::all);

start = 5;  step = 6; end = sig.size() -1 ;
    
iitpave.sigyz =  sig(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() -1, Eigen::all); 

start = 6;  step = 6; end = sig.size() ;
    
iitpave.sigxz =  sig(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() -1, Eigen::all);  

// Strains
start = 1;  step = xl; end = Ee.size() ;

VectorXd Eel = Ee(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() -1); // Reduce vector to evaluation points only

start = 1;  step = xl; end = Ee.size() ;

VectorXd Nul = Nu(VectorXi::LinSpaced((end - start) / step + 1, start, start + ((end - start) / step) * step).array() -1); // Reduce vector to evaluation points only

iitpave.epsx =  (Eel.array().inverse()).array() * (iitpave.sigx.array() - (Nul.array()*(iitpave.sigy + iitpave.sigz).array()).array()).array() ;

iitpave.epsy =  (Eel.array().inverse()).array() * (iitpave.sigz.array() - (Nul.array()*(iitpave.sigx + iitpave.sigy).array()).array()).array() ;

iitpave.epsxy  =(1+Nul.array()).array() * (Eel.array().inverse()).array() * iitpave.sigxy.array();

iitpave.epsyz  =(1+Nul.array()).array() * (Eel.array().inverse()).array() * iitpave.sigyz.array();

iitpave.epsxz  =(1+Nul.array()).array() * (Eel.array().inverse()).array() * iitpave.sigxz.array();

return iitpave;


}


