#include <iostream>
#include <cmath>
#include "TString.h"//includes Root way of making strings
#include "TTree.h"
#include <TTree.h>
#include <TChain.h>

#include "TMath.h"
#include "TLorentzVector.h"
using namespace std;

float fact(int n)
{ int prod = 1;
  for (int j =1 ; j<=n; j++)
    {   prod =prod*j;
        
    }
 return prod ;
}
float RandomFloat(float a, float b) {
   float random = ((float) rand()) / (float) RAND_MAX;
   float diff = b - a;
   float r = random * diff;
   return a + r;
}

void sine()

{ 
       


  TProfile * sin_plot = new TProfile("sin_plot", "sin(x); x  (rad) ; sin(x) (a.u.)", 100, 0, M_PI, -1, +1);
       float xmin = 0; 
       float xmax = M_PI;
       float sinx = 0;
   for (int j = 0 ; j<=1000; j++)
      {
      
       float xvar = RandomFloat( xmin, xmax);
       
      for( int i = 0 ; i < 5 ; i ++ )
	 {
             
             
             sinx =  (pow(-1, i)*pow(xvar, (2*i + 1)))/fact((2*i + 1));
             std::cout<<sinx<<"\t"<<i<<"\n";
             sin_plot->Fill(xvar, sinx, 1);
	 }
      }
       TCanvas * bplot = new TCanvas("bplot", "plot of sin x ", 600, 600);
       sin_plot->Draw("colz");

}
