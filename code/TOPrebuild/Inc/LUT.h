#ifndef LUT_h
#define LUT_h

#include <string>
#include <vector>
#include <map>
#include <fstream>

using namespace std;
class LUT {
  public:
    LUT(string,int nbins=128,int dimensions=1);
    ~LUT();
    double GetValue(int i){ return table.at(i); }
    double GetY(double x);
    int GetNbins(){ return nbins; }
    string GetName(){ return name; }

    int nbins;
    int dimensions;
    string name;
    std::vector<double> table;
    std::map<double,double> table2D;
};

#endif
