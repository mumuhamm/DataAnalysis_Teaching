#include <iostream>
#include <cmath>
#include "TH1D.h"
#include "TRandom.h"
#include "TCanvas.h"
using namespace std;



class box3D{
private:
  double len;
  double brdth;
  double hei;
  double px;
  double py;
public:

  void initData1(double l, double b, double h)
  { 
    len=l;
    brdth=b;
    hei=h;
  }

  void initData2(double momx, double momy)
  {
    px=momx;
    py=momy;
  }

  double area(){return len*brdth;}
  double volume(){return len*brdth*hei;}
  double transmom(){return sqrt(px*px+py*py);}
};






int main()
{
  TH1D * myhist = new TH1D("myhist","p_{T};p_{T} (GeV); Events", 100, 0, 50);
  TRandom * r = new TRandom();
  double px, py, pT;
  box3D mybox;
  for (int i =0;i<1000;i++)
    { px = r->Uniform(0., 20.);
      py = r->Uniform(0., 20.);
      mybox.initData2(px, py);
      pT = mybox.transmom();
      myhist->Fill(pT);
    }
  TCanvas *c = new TCanvas();
  myhist->DrawNormalized();
  c->SaveAs("plot.png");
}
