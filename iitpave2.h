#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <stdexcept>
#include <vector>
#include <limits>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <complex>
#include <algorithm>
#include <cassert>
#include <numeric> 

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using namespace std;
using namespace Eigen;

/*Class declarations*/
class Layer
{
    public:
        double ElasticModulus, PoissonRatio, Thickness;

        Layer() //Constructor
        {
            ElasticModulus = 0.0;
            PoissonRatio = 0.0;
            Thickness = NAN;
        }
};

class Load
{
    public:
        double load, pressure, x0, y0;

        Load()
        {
            load = 0.0;
            pressure = 0.0;
            x0 = 0.0;
            y0 = 0.0;
        }

        double radius()
        {
            try
            {
                return sqrt(load / (M_PI * pressure));
            }
            catch(exception& e)
            {
                cout << "Error in calculating radius! Load = " << load << ", pressure = " << pressure << ". Program exiting.\n";
                cout << "Error message: " << e.what() << endl;
            }
            
        }
};

class Point
{
    public:
        double x0, y0, z0;

    Point()
    {
        x0 = 0.0;
        y0 = 0.0;
        z0 = 0.1;
    }
};

class Pave{
    public:
        string analysis, bond; 
        
        //Numerical parameters
        int N;  //Number of Bessel zero points in numerical integration
        int n; //Number of Gauss points points between zero points.
        double Laml;
        double H;//Bottom depth of last bottom layer
        Eigen::VectorXd kh;  //= [1e6 1e6];  Partial Bonded  

       
        //Pavement material properties (minimum two layers required)
        Eigen::VectorXd zi;        //% Depth of first n-1 layers from the 
                                  //surface [mm]: last z = inf, and should not 
                                 // be added NB: zi(i) > zi(i-1) > z(i-2)...
                                // Layer Young's moduli [MPa]
        Eigen::VectorXd nu;    //Layer Poisson's ratio [-] 
        Eigen::VectorXd E;      //Layer Young's moduli [MPa]
        
        //Load configuration - Single case
        Eigen::VectorXd q;     //Load pressure [MPa] (uniform vertical pressure)
        VectorXd a;    //Load radii [mm] (circular load)
        Eigen::MatrixXd Xl;  //Load positions [mm]: [x1 y1; x2 y2;..xi yi]
        
        // Location of evaluation points: [x1 y1 z1; x2 y2 z2;..] - Single case
       Eigen::MatrixXd Xp; 

       // integration points and weights
      MatrixXd XipWip;
      MatrixXd ABCD; //Layer coffecient matrix A,B,C & D

       //Output Stress 
       Eigen::VectorXd sigz; 
       
       //Output  Strains
       Eigen::VectorXd epsx; 
       Eigen::VectorXd epsy; 
       Eigen::VectorXd epsz; 

       void PopulateLayers(const std::vector<Layer>& layers) {
        // Resize nu to match the size of layers
         nu.resize(layers.size());
         E.resize(layers.size());
         zi.resize(layers.size());

        
        for (int i = 0; i < layers.size(); ++i) {
            nu[i] = layers[i].PoissonRatio;
            E[i] = layers[i].ElasticModulus;
            zi[i] = layers[i].Thickness;
        }
        return;
        }
      void PopulateLoads(const std::vector<Load>& loads) {
        // Resize nu to match the size of layers
       
        
        Xl.resize(loads.size(), 2);
        q.resize(loads.size());

       for (int i = 0; i < loads.size(); ++i) {
       Xl(i, 0) = loads[i].x0; // Corrected to use parentheses and 0-based indexing
       Xl(i, 1) = loads[i].y0; // Corrected to use parentheses and 0-based indexing
       q(i)     = loads[i].pressure;

}
    return;
    }
    
       void PopulatePoints(const std::vector<Point>& points) {
    Xp.resize(points.size(), 3); // Assuming Xd is an Eigen::MatrixXd

    for (int i = 0; i < points.size(); ++i) {
        Xp(i, 0) = points[i].x0;
        Xp(i, 1) = points[i].y0;
        Xp(i, 2) = points[i].z0;
    }


        return;
        }
};


/*Function prototypes*/
void readInputFile(int &numLayers, int &numLoads, int &numPoints, vector<Layer> &layers, vector<Load> &loads, vector<Point> &points,Pave &Gen, string inputfilename);
void init_LET(Pave& iitpave);
MatrixXd numint_coeff(int N, int n, MatrixXd Xp, MatrixXd Xl, double H ,VectorXd a);
MatrixXd besselroots(int o, int N, int k);
MatrixXd arb_func(int n, const VectorXd& zi, const VectorXd& E, const VectorXd& nu, const Pave& iitpave);
void LET_response(Pave& iitpave);
//*MATLAB Function prototypes*/
MatrixXd repmat(int mat, int rows, int cols);
MatrixXd repmat(double mat, int rows, int cols);
MatrixXd repmat(const VectorXd& mat, int rows, int cols);
MatrixXd repmat(const RowVectorXd& mat, int rows, int cols);
MatrixXd repmat(const MatrixXd& mat, int rows, int cols);
MatrixXi repmat(const RowVectorXi& mat, int rows, int cols);
MatrixXi repmat(const VectorXi& mat, int rows, int cols) ;
MatrixXi repmat(const MatrixXi& mat, int rows, int cols) ;
Eigen::SparseMatrix<double> repmat(const Eigen::SparseVector<double>& mat, int rows, int cols);

MatrixXd LinSeq(int j, int k);
MatrixXd LinSeq(int j, int i, int k);

MatrixXd getSubMatrix(const MatrixXd& Mat, const VectorXd& indices);
MatrixXd GetColumn(const MatrixXd& Xp, const std::string& row_index, int column_index);
MatrixXd JoinAllRows(const MatrixXd& mat1, const MatrixXd& mat2);
MatrixXd JoinAllRows(const VectorXd& vec1, const MatrixXd& mat2);

void sortColumns(MatrixXd& matrix);

VectorXd ExtractVecElements(const VectorXd& rho, const VectorXd& rho_non0);
void saveit(const Eigen::MatrixXd& mat);
void saveit(const Eigen::MatrixXi& mat) ;
MatrixXd linearIndexing(const MatrixXd &A, const MatrixXd &B) ;
pair<MatrixXd, MatrixXd> lookup_gauss(int n, const VectorXd& a, const VectorXd& b);

Eigen::SparseMatrix<double> spalloc(int m, int n, int nzmax);
SparseMatrix<double> H_SparseConcatenate(const SparseMatrix<double>& A, const SparseMatrix<double>& B) ;
SparseMatrix<double> V_SparseConcatenate(const SparseMatrix<double>& A, const SparseMatrix<double>& B)  ;
SparseMatrix<double> extractSubmatrix(const Eigen::SparseMatrix<double>& matrix,    const Eigen::VectorXi& rowIndices,   const Eigen::VectorXi& colIndices);
 
MatrixXi ind2sub(int rows, int columns, const VectorXi& linear_indices) ;

SparseMatrix<double> extractSparseSubMat(const SparseMatrix<double>& BC0, 
                                         const VectorXi& row_indices, 
                                         const VectorXi& column_indices) ;
MatrixXd arb_func_interp(int n, const RowVectorXd& xip, const MatrixXd& ABCD);
RowVectorXd linearInterp(const RowVectorXd& x, const RowVectorXd& y, const RowVectorXd& xi);
 