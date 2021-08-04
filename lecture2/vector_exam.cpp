#include <iostream>
#include <vector>
#include <list>
#include <cstdlib>
#include <math.h>
#include "TH1D.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TMath.h"
#include "perceptron.h"
#include <stdio.h>
#include <fstream>
using namespace std;


class box2d
{
private:
  double a;
  int b;
public:
  box2d(double k, int m);
  void baal(double g, int h);
};

  box2d::box2d(double k, int m)
  {
    a=k; 
    b=m;
  }
 // void box2d::baal(double g, int h){return g*h;}

int main()
{
  box2d gh(2.3, 5);
   perceptron rk(0.3,  34);
  vector<vector<double>> mymatrix
    {
      {89.3, 56., 34},
	{45.3, 2.5, 3.2},
	  {3.4, 3.6, 6.}
    };
  for (int i =0 ; i<3; ++i){
  std::cout<<"size of the 2D vector :"<< mymatrix[i].size()<<"\n";
  }}
