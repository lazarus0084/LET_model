#include "iitpave2.h"


MatrixXd numint_coeff(int N, int n, MatrixXd Xp, MatrixXd Xl, double H,VectorXd a){


    VectorXi row_indices, column_indices ;
   
 // Number of load points (xl) and deformation points (xd)
    int xl = Xl.rows(); // Number of load points (rows in Xl)
    int xp = Xp.rows(); // Number of deformation points (rows in Xd)
 
    MatrixXd alpha = repmat(a, xp, 1).array()/H;
    Xl    = repmat(Xl, xp, 1);
    Xp    = repmat(Xp, xl, 1);  //Xd is reorganized below to correspond to the

   MatrixXi LinearSeqence = RowVectorXi::LinSpaced(xp,1,xp);
 // RowVectorXd temp_arg_rowvec = LinSeq(1,xp).transpose();
    MatrixXi dx1 = repmat(LinearSeqence, xl, 1);
    dx1 = dx1.reshaped();
    dx1 = dx1.transpose();
  
   
   LinearSeqence = RowVectorXi::LinSpaced(xl, 0, xp * (xl - 1));

   MatrixXi dx = dx1 + repmat(LinearSeqence,1,xp); dx = dx.array() -1;

   row_indices = dx.transpose();

  // Selecting specific rows and all columns using Eigen::all
   MatrixXd Xp_sub = Xp(row_indices, Eigen::all);
 
   MatrixXd r  = ((Xp_sub.col(0).array() - Xl.col(0).array()).square() +
              (Xp_sub.col(1).array() - Xl.col(1).array()).square()).sqrt();
          
   
  //Introduce the parameter rho (notice: since r is a vector, rho is also  a vector)
    VectorXd rho = r/H; 
    
    MatrixXd B0r = besselroots(0, N, 1);
    MatrixXd B1r = besselroots(1, N, 1);

    int n3 = rho.size();       //Number of radii
    int n4 = alpha.size();     //Number of loads: NB: n3 = n4, due to repmat(..)
                                 //commands above

 // Find non-zero radii/rho
    vector<int> rho_Non0_temp;
    for (int i = 0; i < rho.size(); ++i) {
        if (rho(i) != 0) {
            rho_Non0_temp.push_back(i);
        }
    }
    //Convert std::vector<int> to Eigen::VectorXi for easy Matrix manipulations
    VectorXi rho_Non0 = Eigen::Map<Eigen::VectorXi>(rho_Non0_temp.data(), rho_Non0_temp.size());

// Organize matrix with integration points:
// We want to know the m0 and m1-values for the roots, i.e., where
// J0(m0*r) = 0, J1(m1*r) = 0, and J1(m1*a) = 0. These are
// evaluated by taking the above found roots B0r and B1r as:
// m0 = B0r/r, m1 = B1r/r, m1 = J1/a. The reason why we need the specific m0
// and m1 values is that m-values are included in the integration as
// separate values. In case r = 0, we have a special case where J0(m0*r) =
// constant and J1(m1*r) = constant. So we do not need parameters for this
// in the integration. We only need integration points for the J1(m1*a).
// So we won't need the m values for the J0(m*r) and J1(m*r) functions, only
// for J1(m1*a), since it will never be 0. But below, we need the matrix
// size to be consistent. So we put the values to be zero. This does not
// complicate the computations, since (as shown below) we only choose the
// first N nonzero roots in the integration. And we will always have N
// nonzero roots since the radius, a, will always be nonzero. But we still
// need to organize zeros in the roots below before we can pick N nonzero
// values.
//
// The first xl rows represent the integration points to evaluate the
// deformation of one point due to the xl loads. If a radius is zero,
// division by zero exists.
//
// Define roots matrices with zeros at B0r spaces and B1r roots at the
// remaining spaces.
 
//Concatenate the zerosMatrix and divisionResult horizontally

    MatrixXd roots_z(n3, 2*N);  
    roots_z << MatrixXd::Zero(n3, N), repmat(B1r,n4,1).array() / repmat(alpha,1,N).array(); 

// Insert B0r roots in places where rho > 0 (if rho = 0, B0r/rho = NaN).
   VectorXd rho_sub = rho(rho_Non0);

    if (rho_Non0.size() != 0) {
    //roots_z(rho_Non0, VectorXi::LinSpaced(N, 0, N-1)) = 
 
    roots_z.block(rho_Non0(0), 0, rho_Non0.size(), N) = (repmat(B0r,rho_Non0.size(),1)).array()/(repmat(rho_sub,1,N)).array();
    
    }

// roots_z = sort(roots_z')' Sort each row of roots_z in ascending order

//MatrixXd sorted_roots_z = roots_z.transpose(); // Transpose
for (int i = 0; i < roots_z.rows(); ++i) {
    VectorXd row = roots_z.row(i);
    std::sort(row.data(), row.data() + row.size()); // Sort row
    roots_z.row(i) = row;
}

//roots_z = sorted_roots_z.transpose(); // Transpose back
// Initialize a vector with values from 1 to n3*N
    VectorXd values = VectorXd::LinSpaced(n3 * N, 1,  n3 * N);
    Map<Eigen::MatrixXd> indx(values.data(), n3, N);

    // Find non-zero radii/rho
    vector<int> rowIndices;
    for (int i = 0; i < rho.size(); ++i) {
        if (rho(i) == 0) {
            rowIndices.push_back(i);
        }
    }
    
 VectorXi columnIndices = VectorXi::LinSpaced(indx.cols(), 0, indx.cols()-1);
 VectorXd B(columnIndices.size());B.setOnes(); B = B*n3*N;

 for (int i = 0; i < rowIndices.size(); ++i) {
        indx.row(rowIndices[i]) += B.transpose();
    }
rowIndices.clear();assert(rho.size() == alpha.size());
double epsilon = 1e-8;
for (int i = 0; i < rho.size(); ++i) {
       if (fabs(rho(i) - alpha(i)) < epsilon) {
    rowIndices.push_back(i);
}
}

B = VectorXd::LinSpaced(N, 0, N-1);; B = B*n3;

 for (int i = 0; i < rowIndices.size(); ++i) {
       indx.row(rowIndices[i]) += B.transpose(); 

    }

VectorXd arg1 = VectorXd::Zero(n3);
MatrixXd seq_z = JoinAllRows(arg1,linearIndexing(roots_z,indx));

// Distribution of n integration points and weights of the first interval 
auto result = lookup_gauss(n, seq_z.col(0),seq_z.col(1));
MatrixXd xip_z1 = result.first;  
MatrixXd wip_z1 = result.second;


//Distribution of 20 integration points and weights of the second 

result = lookup_gauss(20,seq_z.col(1) ,seq_z.col(2));
MatrixXd xip_z2 = result.first;  
MatrixXd wip_z2 = result.second;

//Distribution of 10 integration points and weights of the third interval

result = lookup_gauss(10, seq_z.col(2), seq_z.col(3));
MatrixXd xip_z3 = result.first;  
MatrixXd wip_z3 = result.second;



// Extract a range of columns (e.g., columns 2 to 5, zero-based indexing)
    int startColumn = 3; // reduce 1 from the val in  MATLAB
    int numColumns = seq_z.cols()- startColumn- 1;



    // Perform the block operation
    MatrixXd seq_z1 = seq_z.block(0, startColumn, seq_z.rows(), numColumns).transpose();
   
    startColumn = 4; // reduce 1 from the val in  MATLAB
    numColumns = seq_z.cols()- startColumn;
    MatrixXd seq_z2 = seq_z.block(0, startColumn, seq_z.rows(), numColumns).transpose();
    

     //Distribution of 5 integration points and weights of the remaining intervals
    VectorXd Arg1 = Map<VectorXd>(seq_z1.data(), seq_z1.size());
    VectorXd Arg2 = Map<VectorXd>(seq_z2.data(), seq_z2.size());
    result = lookup_gauss(5, Arg1,Arg2);
    MatrixXd xip_z = result.first; 
    MatrixXd wip_z = result.second; 

 
    xip_z.transposeInPlace();
    xip_z = Eigen::Map<Eigen::MatrixXd>(xip_z.data(), xip_z.size() / n3, n3);xip_z.transposeInPlace();
    
 
    wip_z.transposeInPlace();
    wip_z = Eigen::Map<Eigen::MatrixXd>(wip_z.data(), wip_z.size() / n3, n3);wip_z.transposeInPlace();

    // Collect all integration points Concatenate matrices horizontally
    MatrixXd xip_r(xip_z1.rows(),xip_z1.cols() + xip_z2.cols() + xip_z3.cols() + xip_z.cols());
    xip_r << xip_z1, xip_z2, xip_z3, xip_z;
     // Collect all integration points Concatenate matrices horizontally
    MatrixXd wip_r(wip_z1.rows(),wip_z1.cols() + wip_z2.cols() + wip_z3.cols() + wip_z.cols());
    wip_r << wip_z1, wip_z2, wip_z3, wip_z;

    MatrixXd XipWip(2*xip_r.rows()+2*wip_r.rows(), xip_r.cols());

    XipWip << xip_r,wip_r,xip_r,wip_r;
    
    return XipWip;
  
}
