#include <string>
#include <vector>
#include <map>
#include <fstream>

#include "DircData.h"
#include "DircRecSvc.h"
#include "DircRecAlg.h"

using namespace TMath;
using namespace std;

int main(int argc,char** argv)
{
  string input  = "data/Pi.root";
  string output = "result/DircRec.root";
  if(argc!=1) input  = argv[1];
  if(argc>=3) output = argv[2];
  cout<<"Input  file: "<<input <<endl;
  cout<<"Output file: "<<output<<endl;
  DircRecAlg* Dirc = new DircRecAlg(input,output);
  Dirc -> Loop();
  delete Dirc;
  return 0;
}
