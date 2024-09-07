#include "iitpave2.h"

MatrixXd arb_func_interp(int n, const RowVectorXd& xip, const MatrixXd& ABCD) {

//-------------------------------------------------------------------------
// DESCRIPTION:
// This function utilizes a cubic spline interpolation scheme to return 
// interpolated values for sought integration points (xip) based on the 
// calculated solutions for the coefficients of integration (arbitrary 
// functions).

// INPUT PARAMETERS:
// n:    Number of layers (including half space layer), i.e. n = length(E)
// xip:  Sought integration points
// ABCD: Arbitrary function values (soultion points) to be interpolated.
//-------------------------------------------------------------------------
    // Initialize the output matrix with the same number of rows as ABCD and columns as xip
    MatrixXd ABCDi(4*n,xip.cols()); //(ABCD.rows(), xip.size());
    
    RowVectorXd m = ABCD.bottomRows(1);

    // Extract the relevant rows from ABCD
    MatrixXd sub_ABCD = ABCD.topRows(ABCD.rows()-1);

    // Interpolate each row
    for (int i = 0; i < 4*n; ++i) { 
        // Linear interpolation of the i-th row of selectedABCD
        ABCDi.row(i) = linearInterp(m, sub_ABCD.row(i), xip);
    }
    
    return ABCDi;
}

