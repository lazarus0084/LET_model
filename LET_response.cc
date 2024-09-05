#include "iitpave2.h"

void LET_response(Pave& iitpave) {
    cout << "control at LET_response" << endl;
    int N, n ;
    double H, Laml;
    VectorXd E, nu, zi, q, a(iitpave.Xl.rows()); a = 150.8 * VectorXd::Ones(iitpave.Xl.rows());   iitpave.a = a;
    MatrixXd Xp, Xl, XipWip,ABCD;
   
 
   // INPUT PARAMETERS
    N = iitpave.N; //Number of Bessel roots in integration
    n = iitpave.n; //Number of integration points in between each Bessel root
    Xp = iitpave.Xp; //Output coordinates
    Xl = iitpave.Xl; //Load posistion coordinates
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
    int xl = Xl.rows(); // Number of load points
    int xp = Xp.rows(); // Number of deformation points
// Check if the length of q and a are correct
     if (q.size() < xl || a.size() < xl) {
      cout << "Number of q and/or a values is lower than the number of load coordinates (= number of loads), so all loads are given the same q- and a-values" << endl;
        q = VectorXd::Constant(xl, q(0));
        a = VectorXd::Constant(xl, a(0));
    }
// Increase Xl, q, alpha, Xd and z in order to estimate the radius from each
// of the loads to all of the points. The number of radii are: xl*xd
    Xl = repmat(Xl, xp, 1);   // Extends with the number of evaluation
                             // points
    q  = repmat(q, xp, 1);    //-||-
    MatrixXd alpha = repmat(a, xp, 1).array()/H;  // -||-
    Xp    = repmat(Xp, xl, 1);   // -||- <-- is reorganized below to 
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

    Sequence = Eigen::RowVectorXi::LinSpaced(xl, 0, xp * (xl - 1));

    RowVectorXi dx = dx1.reshaped(1,dx1.size()) + repmat(Sequence,1,xp); dx = dx.array() -1;

     MatrixXd Xp1 = Xp(dx , all);
     
     VectorXd z = Xp1(Eigen::all, 2);

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


   for (size_t i = 0; i < layer_no.size(); ++i) {
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

    cout << Lami1 << endl;
    MatrixXd dense = Lami1;
    saveit(dense);

}


