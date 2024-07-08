#include "iitpave2.h"

MatrixXd numint_coeff(int N, int n, MatrixXd Xp, MatrixXd Xl, double H,float a){
   myfunctions my;
 // Number of load points (xl) and deformation points (xd)
    int xl = Xl.rows(); // Number of load points (rows in Xl)
    int xp = Xp.rows(); // Number of deformation points (rows in Xd)
 
    MatrixXd alpha = my.repmat((a/H), xp, 1);
    Xl    = my.repmat(Xl, xp, 1);
    Xp    = my.repmat(Xp, xl, 1);  //Xd is reorganized below to correspond to the
    
   RowVectorXd temp_arg_rowvec = my.LinSeq(1,xp).transpose();
   MatrixXd dx = my.repmat(temp_arg_rowvec,xl,1);
 
   MatrixXd dx1 = dx + my.repmat(my.LinSeq(0,xp,xp*(xl-1)),xp,1).transpose();
   dx = dx1;
   MatrixXd Xpp = my.getSubMatrix(Xp, dx.transpose()); // using transpose() for row selection takes a vector input (column)
   Xp = Xpp;
   MatrixXd Arr = my.GetColumn(Xp, ":", 1) - my.GetColumn(Xl, ":", 1);
   VectorXd r = ((Xp.col(0) - Xl.col(0)).array().square() + (Xp.col(1) - Xl.col(1)).array().square()).sqrt();
   //Introduce the parameter rho (notice: since r is a vector, rho is also  a vector)
   VectorXd rho = r/H; 
   MatrixXd B0r = besselroots(0, N, 1);
   MatrixXd B1r = besselroots(1, N, 1);
   double n3 = rho.size();        //Number of radii
   double n4 = alpha.size();     //Number of loads: NB: n3 = n4, due to repmat(..)
                                //commands above

    // Find non-zero radii/rho
    vector<double> rho_Non0_temp;
    for (int i = 0; i < rho.size(); ++i) {
        if (rho(i) != 0) {
            rho_Non0_temp.push_back(i+1);
        }
    }
    VectorXd rho_Non0(rho_Non0_temp.size());
    for (size_t i = 0; i < rho_Non0_temp.size(); ++i) {
        rho_Non0(i) = rho_Non0_temp[i];
    }
    MatrixXd mat1;mat1.resize(n3, N); mat1.setZero(); //(rho(rho_Non0) selects indices of 
    MatrixXd mat2 = my.repmat(B0r,n4,1);
    MatrixXd mat3 = my.repmat(my.ExtractVecElements(rho, rho_Non0),1,N);
    //MatrixXd roots_z = mat2.array()*mat3.array();
  
       
    


   
    



    //MatrixXd roots_z = my.JoinAllRows(zerosMatrix , (my.repmat(B0r,rho_Non0.size(),1).array()* my.repmat(Extracted,1,N).array()));
       #if 0
      MatrixXd mat1;mat1.resize(n3, N); mat1.setZero();
      MatrixXd mat2 = my.repmat(B1r,n4,1);
  


    MatrixXd mat1;mat1.resize(n3, N); mat1.setZero();
    MatrixXd mat2 = my.repmat(B1r,n4,1).array()/my.repmat(alpha,1,N).array();
    MatrixXd roots_z = my.JoinAllRows(mat1,mat2);
    VectorXd Extracted = my.ExtractVecElements(rho, rho_Non0);
    MatrixXd m11 = my.repmat(B1r,n4,1);
    MatrixXd m2 = my.repmat(Extracted,1,N);
    MatrixXd m3 = m1.array()/m2.array();
    
    // Create a matrix of zeros with size n3 x N
    MatrixXd zerosMatrix = MatrixXd::Zero(n3, N);


    // Check if rho_Non0 is not empty
    if (!rho_Non0.empty()) {
          my.repmat(B0r,rho_Non0.size(),1)./repmat()
        // Add your code here if rho_Non0 is not empty
    }

     roots_z(rho_Non0,1:N) = repmat(B0r,[length(rho_Non0),1])...
    ./repmat(rho(rho_Non0),[1,N]);
    //roots_z = [zeros(n3,N) repmat(B1r,[n4,1])./repmat(alpha,[1,N])];

 
   
// Organize matrix with integration points:
// We want to know the m0 and m1-values for the roots, i.e. where
// J0(m0*r) = 0, J1(m1*r) = 0 and J1(m1*a) = 0.  These are
// evaluated by taking the above found roots B0r and B1r and as:
// m0 = B0r/r, m1 = B1r/r, m1 = J1/a. The reason why we need the specific m0
// and m1 values is that m-values are included in the integration as
// separate values. In case r = 0, we have a special case where J0(m0*r) =
// constant and J1(m1*r) = constant. So we do not need parameters for this
// in the integration. We only need integration points for the J1(m1*a).
// So we won't need the m values for the J0(m*r) and J1(m*r) functions. Only
// for J1(m1*a), since it will never be 0. But, below we need the matrix
// size to be consistent. So we put the values to be zero. This does not
// complicate the computations, since (as shown below) we only choose the
// first N nonzero roots in the integration. And we will always have N
// nonzero roots since the radius, a, will always be nonzero. But we still
// need to organize zeros in the roots below, before we can pick N nonzero
// values.
//
// The first xl rows represents the integration points to evaluate the
// deformation of one point due to the xl loads. If a radius is zero,
// division by zero exist

// Define roots matrices with zeros at B0r spaces and B1r roots at the
// remaining spaces
roots_z = [zeros(n3, N) repmat(B1r, [n4, 1])./repmat(alpha, [1, N])];

// Insert B0r roots in places where rho > 0 (if rho = 0, B0r/rho = NaN).
if isempty(rho_Non0) != 1
    roots_z(rho_Non0, 1:N) = repmat(B0r, [length(rho_Non0), 1])...
    ./repmat(rho(rho_Non0), [1, N]);
end

// Organize roots in each row in ascending order
roots_z = sort(roots_z')';

// Select the points including up to the first N nonzero zero-value points
// (if we include more than N nonzero zero-value points, we might miss some
// zero points in between <-- because we divide roots by alpha and rho, the
// m-value ranges are different.)

// Define indx vector to select the nonzero roots
indx = zeros(n3, N); // Organize zero matrix
indx(:) = 1:n3*N;    // Fill up matrix with indices that assume that no
                     // radius is zero (meaning that roots are present)

// Replace the indices where the radius is zero with new indices where
// integration points are nonzero
indx(find(rho==0), :) = indx(find(rho==0), :) + n3*N;

// A special case can appear where alpha = rho, i.e., where the deformation
// point is at the edge of the loading. In this case the
// roots_r = [Br/rho Br/alpha], where Br/rho = Br/alpha. So we have the same
// integration points appearing twice. These are arranged side by side after
// we have used the 'sort'-function. So in this case we have to replace the
// index with another index to ensure we do not take out the same values
// twice. This is done below where we add a vector to each of these rows
// given as [0 1 2 3 .... N-2 N-1]*n3, i.e., the columns are shifted (we
// pick out every second column instead of the first N columns).

   
cout << "Matrix non 0 rho:\n" << rho_Non0 << endl;
 = [zeros(n3,N) repmat(B1r,[n4,1])./repmat(alpha,[1,N])];
 

  

  


 
  cout << "Value rho:" << endl << rho <<  endl;
  
    // Print rho_Non0 vector
    cout << "rho_Non0 vector:" << endl;
    for (int i = 0; i < rho_Non0_temp.size(); ++i) {
        cout << rho_Non0_temp(i) << endl;
    }
  //rho(rho_Non0)
  vector<double> ExtractedElements;

    // Adjust the indices and push them to the ExtractedElements vector
    for (int i = 0; i < rho_Non0_temp.size(); ++i) {
        ExtractedElements(i) = rho(rho_Non0_temp(i - 1));
    }

    // Print the ExtractedElements vector
    cout << "Extracted elements with adjusted indices: ";
    for (int i = 0; i < ExtractedElements.size(); ++i) {
        cout << ExtractedElements[i] << " ";
    }
     cout << "Size of roots_z : " << roots_z.rows() << "x" << roots_z.cols() << endl;
    //cout << "roots_z : " << roots_z << endl;
    cout << "Element at index (5, 6): " << roots_z(5, 6) << endl;
    cout << "Size of roots_z : " << roots_z.rows() << "x" << roots_z.cols() << endl;
   
   #endif

  cout << "Size of mat2: " << mat2.rows() << " x " << mat2.cols() << endl;
  cout << "Size of mat3: " << mat3.rows() << " x " << mat3.cols() << endl;
   
   
   return MatrixXd::Identity(1, 1);

}

