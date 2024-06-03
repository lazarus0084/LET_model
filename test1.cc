#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<math.h>
#include<limits>

using namespace std;

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

int main()
{
    ifstream ifile;
    int imax;
    string buf, heading;
    vector<Layer> layers;
    Layer tempLayer;

    ifile.open("input1.txt", ios::in);

    if(!ifile.is_open())
    {
        cout << "File not found!\n";
        return 0;
    }

    //Read *LAYERS section
    getline(ifile, heading);
    if(!heading.compare("*LAYERS"))
    {
        cout << "LAYERS section not found!\n";
        return 0;
    }
    getline(ifile, buf);
    sscanf(buf.c_str(), "%d", &imax);

    for(int i = 0; i < imax; i++)
    {
        if(i != imax - 1)
        {
            getline(ifile, buf);
            sscanf(buf.c_str(), "%lf,%lf,%lf", &tempLayer.ElasticModulus, &tempLayer.PoissonRatio, &tempLayer.Thickness);
            layers.push_back(tempLayer);
        }
        else
        {
            getline(ifile, buf);
            sscanf(buf.c_str(), "%lf,%lf", &tempLayer.ElasticModulus, &tempLayer.PoissonRatio);
            layers.push_back(tempLayer);

        }
    }

    cout << layers[imax-1].PoissonRatio << endl;

    return 0;

}