#include "iitpave2.h"

int main(int argc,char* argv[])
{
    int numLayers, numLoads, numPoints;
    vector<Layer> layers;
    vector<Load> loads;
    vector<Point> points;
    
    //Read inputs
    readInputFile(numLayers, numLoads, numPoints, layers, loads, points, argv[1]);

    //Tests
    cout << layers[layers.size()-1].PoissonRatio << endl;
    cout << points[0].x0 << endl;
    cout << loads[0].radius() << endl;

    return 0;

}