#include "iitpave2.h"

void readInputFile(int &numLayers, int &numLoads, int &numPoints, vector<Layer> &layers, vector<Load> &loads, vector<Point> &points, string inputfilename)
{
    ifstream ifile;
    string buf, heading;
    Layer tempLayer;
    Load tempLoad;
    Point tempPoint;

    try
    {
        ifile.open(inputfilename, ios::in);

        if(!ifile.is_open())
        {
            throw runtime_error("Failed to open " + inputfilename + ". Is the file and path correct?\n");
        }

        /*Read file*/
        //Read *LAYERS section
        getline(ifile, heading);
        if(heading.compare("*LAYERS") != 0)
        {
            throw runtime_error("*LAYERS section not found input file! It must be the first section of the file.\n");
        }
        getline(ifile, buf);
        sscanf(buf.c_str(), "%d", &numLayers);

        for(int i = 0; i < numLayers; i++)
        {
            if(i != numLayers - 1) //Except for last layer
            {
                getline(ifile, buf);
                sscanf(buf.c_str(), "%lf,%lf,%lf", &tempLayer.ElasticModulus, &tempLayer.PoissonRatio, &tempLayer.Thickness);
                layers.push_back(tempLayer);
            }
            else //Last layer - no thickness
            {
                getline(ifile, buf);
                sscanf(buf.c_str(), "%lf,%lf", &tempLayer.ElasticModulus, &tempLayer.PoissonRatio);
                layers.push_back(tempLayer);
            }
        }

        //Read *LOADS section
        getline(ifile, heading);
        if(heading.compare("*LOADS") != 0)
        {
            throw runtime_error("*LOADS section not found input file! It must be the second section of the file.\n");
        }
        getline(ifile, buf);
        sscanf(buf.c_str(), "%d", &numLoads);

        for(int i = 0; i < numLoads; i++)
        {
            getline(ifile, buf);
            sscanf(buf.c_str(), "%lf,%lf,%lf,%lf", &tempLoad.load, &tempLoad.pressure, &tempLoad.x0, &tempLoad.y0);
            loads.push_back(tempLoad);
        }

        //Read *POINTS section
        getline(ifile, heading);
        if(heading.compare("*POINTS") != 0)
        {
            throw runtime_error("*POINTS section not found input file! It must be the third section of the file.\n");
        }
        getline(ifile, buf);
        sscanf(buf.c_str(), "%d", &numPoints);

        for(int i = 0; i < numPoints; i++)
        {
            getline(ifile, buf);
            sscanf(buf.c_str(), "%lf,%lf,%lf", &tempPoint.x0, &tempPoint.y0, &tempPoint.z0);
            points.push_back(tempPoint);
        }

        cout << "Finished reading input file.\n";
    }
    catch(exception& e)
    {
        cout << "Program encountered an error reading the input file. The error message is: " << e.what() << "\nProgram exiting.\n";
        exit(1);
    }
    
}