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
        Eigen::VectorXd a;    //Load radii [mm] (circular load)
        Eigen::MatrixXd Xl;  //Load positions [mm]: [x1 y1; x2 y2;..xi yi]
        
        // Location of evaluation points: [x1 y1 z1; x2 y2 z2;..] - Single case
       Eigen::MatrixXd Xp; 


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
class myfunctions{ //MATLAB functions
    public :
// Function overload for int
MatrixXd repmat(int mat, int rows, int cols) {
    return MatrixXd::Constant(rows, cols, static_cast<double>(mat));
}

// Function overload for double
MatrixXd repmat(double mat, int rows, int cols) {
    return MatrixXd::Constant(rows, cols, mat);
}

// Function overload for VectorXd
MatrixXd repmat(const VectorXd& mat, int rows, int cols) {
    int size = mat.size();
    MatrixXd result(size * rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.block(i * size, j, size, 1) = mat;
        }
    }

    return result;
}


// Function overload for RowVectorXd
MatrixXd repmat(const RowVectorXd& mat, int rows, int cols) {
    int size = mat.size();
    MatrixXd result(rows, cols * size);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.block(i, j * size, 1, size) = mat;  // Use block operation to assign the row vector
        }
    }

    return result;
}

// Function overload for MatrixXd
MatrixXd repmat(const MatrixXd& mat, int rows, int cols) {
    int numRows = mat.rows();
    int numCols = mat.cols();
    MatrixXd result(numRows * rows, numCols * cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.block(i * numRows, j * numCols, numRows, numCols) = mat;
        }
    }

    return result;
}

// Function to generate a sequence j:k
MatrixXd LinSeq(int j, int k) {
    int m = k - j + 1;
    MatrixXd result(m, 1);
    for (int idx = 0; idx < m; ++idx) {
        result(idx, 0) = j + idx;
    }
    return result;
}

// Function to generate a sequence j:i:k
MatrixXd LinSeq(int j, int i, int k) {
    int m = static_cast<int>((k - j) / i) + 1;
    MatrixXd result(m, 1);
    for (int idx = 0; idx < m; ++idx) {
        result(idx, 0) = j + idx * i;
    }
    return result;
}
MatrixXd getSubMatrix(const MatrixXd& Mat, const VectorXd& indices) {
    int numRows = indices.size();
    int numCols = Mat.cols();
    MatrixXd submatrix(numRows, numCols);
    submatrix.setZero();

    VectorXd adjustedIndices = indices.array() - 1; // Adjust indices by subtracting 1

    // Extract rows based on adjustedIndices
    for (int i = 0; i < numRows; ++i) {
        int rowIdx = static_cast<int>(adjustedIndices(i));
        if (rowIdx >= 0 && rowIdx < Mat.rows()) {
            submatrix.row(i) = Mat.row(rowIdx);
        } else {
            cerr << "Row index " << rowIdx << " out of bounds." << endl;
        }
    }

    return submatrix;
}
MatrixXd GetColumn(const MatrixXd& Xp, const std::string& row_index, int column_index) {
    int rows = Xp.rows();
    int cols = Xp.cols();
    MatrixXd result;

    if (row_index == ":") {
        // Get all rows for a specific column
        result.resize(rows, 1);
        result = Xp.col(column_index - 1);  // Adjust column index by -1
    } else {
        // Get a specific row and column
        int row = std::stoi(row_index) - 1; // Adjust row index by -1
        result.resize(1, 1);
        result(0, 0) = Xp(row, column_index - 1);  // Adjust column index by -1
    }

    return result;
}




};



/*Function prototypes*/
void readInputFile(int &numLayers, int &numLoads, int &numPoints, vector<Layer> &layers, vector<Load> &loads, vector<Point> &points,Pave &Gen, string inputfilename);

void init_LET(Pave& iitpave);
MatrixXd numint_coeff(int N, int n, MatrixXd Xp, MatrixXd Xl, double H ,float a);
MatrixXd besselroots(int o, int N, int k);
