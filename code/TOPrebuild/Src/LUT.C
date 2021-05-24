#include "LUT.h"

using namespace std;
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
