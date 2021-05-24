using namespace TMath;

double A1=0.6961663, B1=0.0684043;
double A2=0.4079426, B2=0.1162414;
double A3=0.8974794, B3=9.8961610;

double GetNp(double WaveLength) //unit of WaveLength: nm
{
  WaveLength = WaveLength/1000.; //nm to um
  return Sqrt(1.+A1*WaveLength*WaveLength/(WaveLength*WaveLength-B1*B1)
      +A2*WaveLength*WaveLength/(WaveLength*WaveLength-B2*B2)
      +A3*WaveLength*WaveLength/(WaveLength*WaveLength-B3*B3));
}

double GetDnp(double WaveLength) //unit of WaveLength: nm
{
  return (GetNp(WaveLength+0.001)-GetNp(WaveLength-0.001))/0.002;
}

double GetNg(double WaveLength) //unit of WaveLength: nm
{
  double np=GetNp(WaveLength);
  return np/(1+WaveLength/np*GetDnp(WaveLength));
}

double GetEnergy(double WaveLength) //unit of WaveLength: nm
{
  return 1240.7/WaveLength;
}

double GetWaveLength(double Energy) //unit of Energy: eV
{
  return 1240.7/Energy;
}

void refractiveIndex()
{
  const int num = 400;
  double WaveLength[num],Energy[num],Np[num],Ng[num];

  double Ebegin = GetEnergy(800.), Eend = GetEnergy(150.);
  double dE = (Eend-Ebegin)/num;

  ofstream out("n.txt");
  out.setf(ios::fixed);

  for(int i=0;i<num;i++)
  {
    Energy[i]     = Ebegin + dE*i;
    WaveLength[i] = GetWaveLength(Energy[i]);
    Np[i]         = GetNp(WaveLength[i]);
    Ng[i]         = GetNg(WaveLength[i]);

    out<<setprecision(8)<<Np[i]<<" "<<Ng[i]<<std::endl;
  }
  out.close();
  cout<<GetNp(300)<<" "<<GetNp(403)<<" "<<GetNp(650)<<endl;
  return;
}
