#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <stdexcept>
#include <vector>
#include <limits>
#include <math.h>
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using namespace std;

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

/*Function prototypes*/
void readInputFile(int &numLayers, int &numLoads, int &numPoints, vector<Layer> &layers, vector<Load> &loads, vector<Point> &points, string inputfilename);