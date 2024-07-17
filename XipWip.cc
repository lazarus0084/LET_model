#include "iitpave2.h"


MatrixXd numint_coeff(int N, int n, MatrixXd Xp, MatrixXd Xl, double H,float a){
   
 // Number of load points (xl) and deformation points (xd)
    int xl = Xl.rows(); // Number of load points (rows in Xl)
    int xp = Xp.rows(); // Number of deformation points (rows in Xd)
 
    MatrixXd alpha = repmat((a/H), xp, 1);
    Xl    = repmat(Xl, xp, 1);
    Xp    = repmat(Xp, xl, 1);  //Xd is reorganized below to correspond to the
    
   RowVectorXd temp_arg_rowvec = LinSeq(1,xp).transpose();
   MatrixXd dx = repmat(temp_arg_rowvec,xl,1);
 
   MatrixXd dx1 = dx + repmat(LinSeq(0,xp,xp*(xl-1)),xp,1).transpose();
   dx = dx1;
   MatrixXd Xpp = getSubMatrix(Xp, dx.transpose()); // using transpose() for row selection takes a vector input (column)
   Xp = Xpp;
   MatrixXd Arr = GetColumn(Xp, ":", 1) - GetColumn(Xl, ":", 1);
   VectorXd r = ((Xp.col(0) - Xl.col(0)).array().square() + (Xp.col(1) - Xl.col(1)).array().square()).sqrt();
   //Introduce the parameter rho (notice: since r is a vector, rho is also  a vector)
   VectorXd rho = r/H; 
   MatrixXd B0r = besselroots(0, N, 1);
   MatrixXd B1r = besselroots(1, N, 1);
   int n3 = rho.size();       //Number of radii
   int n4 = alpha.size();     //Number of loads: NB: n3 = n4, due to repmat(..)
                                 //commands above

    // Find non-zero radii/rho
    vector<int> rho_Non0;
    for (int i = 0; i < rho.size(); ++i) {
        if (rho(i) != 0) {
            rho_Non0.push_back(i);
        }
    }
    MatrixXd mat1;
    mat1.resize(n3, N); mat1.setZero();
    MatrixXd mat2 = repmat(B1r,n4,1).array()/repmat(alpha,1,N).array();
    MatrixXd roots_z = JoinAllRows(mat1,mat2);
   
    
    VectorXd arg1 = rho(rho_Non0);
    VectorXi columnIndices = VectorXi::LinSpaced(N, 0, N-1);
    //rowIndices = rho_Non0 ;
    if (rho_Non0.size() != 0) {
      roots_z(rho_Non0, columnIndices) = (repmat(B0r,rho_Non0.size(),1)).array()/(repmat(arg1,1,N)).array();
    }
     // Sort each row of roots_z in ascending order
    for (int i = 0; i < roots_z.rows(); ++i) {
        VectorXd row = roots_z.row(i);
        sort(row.data(), row.data() + row.size());
        roots_z.row(i) = row;
    }
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
    
 columnIndices = VectorXi::LinSpaced(indx.cols(), 0, indx.cols()-1);
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

arg1 = VectorXd::Zero(n3);
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

    //Reorganize integration points and weights - Every row
    //(after reorganization) now refers to a single radius/point-to-load case
    xip_z.transposeInPlace();     //First we transform the matrices
    
         
    Arg1 = Map<VectorXd>(xip_z.data(), xip_z.size());
    xip_z = Map<MatrixXd>(Arg1.data(),n3 ,Arg1.size() / n3 );
    
     // Concatenate matrices horizontally
    MatrixXd xip_zConctd(xip_z1.rows(), xip_z1.cols() + xip_z2.cols() + xip_z3.cols() + xip_z.cols());
    xip_zConctd << xip_z1, xip_z2, xip_z3, xip_z;
    xip_z = xip_zConctd;
    MatrixXd xip_r = xip_zConctd;
   
     
    //same of wip_z
    wip_z.transposeInPlace();     //First we transform the matrices
         
    Arg1 = Map<VectorXd>(wip_z.data(), wip_z.size());

            // Perform the reshaping in a single line
    wip_z = Map<MatrixXd>(Arg1.data(),n3 , Arg1.size() / n3 );

    // Concatenate wip_z1, wip_z2, wip_z3, and wip_z horizontally
    MatrixXd wip_zConctd(wip_z1.rows(), wip_z1.cols() + wip_z2.cols() + wip_z3.cols() + wip_z.cols());
    wip_zConctd << wip_z1, wip_z2, wip_z3, wip_z;
    wip_z = wip_zConctd;
    MatrixXd wip_r = wip_zConctd;


   //Concatenate matrices vertically
    MatrixXd XipWip_returning(xip_z.rows() + wip_z.rows() + xip_r.rows() + wip_r.rows(), xip_z.cols());
    XipWip_returning << xip_z,
              wip_z,
              xip_r,
              wip_r;
    // std::cout << "Size of xip_z: " << xip_z.rows() << "x" << xip_z.cols() << std::endl;
    // std::cout << "Size of wip_z: " << wip_z.rows() << "x" << wip_z.cols() << std::endl;
    // std::cout << "Size of xip_r: " << xip_r.rows() << "x" << xip_r.cols() << std::endl;
    // std::cout << "Size of wip_r: " << wip_r.rows() << "x" << wip_r.cols() << std::endl;
    // Display size of XipWip
    //cout << "Size of XipWip: " << XipWip.rows() << " rows x " << XipWip.cols() << " columns." << endl;

  
return XipWip_returning;
  
}
