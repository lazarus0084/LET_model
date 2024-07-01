#include "iitpave2.h"

int main(int argc,char* argv[])
{
    int numLayers, numLoads, numPoints;
    vector<Layer> layers;
    vector<Load> loads;
    vector<Point> points;
    Pave iitpave;
    
   
    
    //Read inputs
    readInputFile(numLayers, numLoads, numPoints, layers, loads, points, iitpave, argv[1]);
    //iitpave object populating functons
    iitpave.PopulateLayers(layers);
    iitpave.PopulateLoads(loads);
    iitpave.PopulatePoints(points);

    
    
    init_LET(iitpave);//filename LET
  
    
 

    
    return 0;
    
 }
    
   
    

    

  

    

    
 
    
    

