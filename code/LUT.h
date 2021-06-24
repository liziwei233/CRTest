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
LUT::LUT(string a,int n,int d):nbins(n),dimensions(d),name(a)
{
  table.clear();
  table2D.clear();
  ifstream Input(name.c_str());
  if(dimensions==1)
  {
    double tmp;
    for(int i=0;i<nbins;i++)
    {
      Input>>tmp;
      table.push_back(tmp);
    }
    Input.close();
  }
  else if(dimensions==2)
  {
    double tmpx,tmpy;
    for(int i=0;i<nbins;i++)
    {
      Input>>tmpx>>tmpy;
      table2D[tmpx]=tmpy;
    }
    Input.close();
  }
}

LUT::~LUT(){}

double LUT::GetY(double x)
{
  std::map<double,double>::iterator itlow=table2D.lower_bound(x);
  if(itlow == table2D.begin()) return itlow->second;
  if(itlow == table2D.end()) return (--itlow)->second;
  else {
    double x1=itlow->first, y1=itlow->second;
    itlow--;
    double x2=itlow->first, y2=itlow->second;
    return y1+(y2-y1)*(x-x1)/(x2-x1);
  }
}
#endif
