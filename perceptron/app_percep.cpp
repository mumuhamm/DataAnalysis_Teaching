#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <math.h>
#include "TH1D.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TMath.h"

#include "perceptron.h"


using namespace std;
//using namespace perceptron;

vector< vector<double> > getIrisX();
vector<double> getIrisy();


vector<double> getIrisy()
{
   vector<double> y;
   
   ifstream inFile;
   inFile.open("data/y.data");
   string sampleClass;
   for (int i = 0; i < 100; i++)
      {
      inFile >> sampleClass;
      if (sampleClass == "Iris-setosa")
         {
         y.push_back(-1);
         }
      else
         {
         y.push_back(1);
         }
      }
   
   return y;
}

vector< vector<double> > getIrisX()
{
   ifstream af;
   ifstream bf;
   ifstream cf;
   ifstream df;
   af.open("data/a.data");
   bf.open("data/b.data");
   cf.open("data/c.data");
   df.open("data/d.data");
   
   vector< vector<double> > X;
   
   for (int i = 0; i < 100; i++)
      {
      char scrap;
      int scrapN;
      af >> scrapN;
      bf >> scrapN;
      cf >> scrapN;
      df >> scrapN;
      
      af >> scrap;
      bf >> scrap;
      cf >> scrap;
      df >> scrap;
      double a, b, c, d;
      af >> a;
      bf >> b;
      cf >> c;
      df >> d;
      X.push_back(vector < double > {a, b, c, d});
      }

   af.close();
   bf.close();
   cf.close();
   df.close();
   
   return X;
}


int main()
{
   std::cout<<"My name is ali"<<"\n";
   vector< vector<double> > X = getIrisX();
    vector<double> y = getIrisy();
    vector<double> test1;
    test1.push_back(5.0);
    test1.push_back(3.3);
    test1.push_back(1.4);
    test1.push_back(0.2);
    
    vector<double> test2;
    test2.push_back(6.0);
    test2.push_back(2.2);
    test2.push_back(5.0);
    test2.push_back(1.5);
    for (int i = 0; i< X.size(); i++){for (int j=0; j<X.size(); j++){ std::cout<<" X matrix val : "<< X[i][j]<<" "<<"\n";}}
    for (int i = 0; i < y.size(); i++){ cout << y[i] << " "; }cout << endl;
    
    perceptron calculate(0.1, 14);
    calculate.perceptron::fit(X, y);
    calculate.printErrors();
    cout << "Now Predicting: 5.0,3.3,1.4,0.2(CorrectClass=-1,Iris-setosa) -> " << calculate.predict(test1) << endl;
    cout << "Now Predicting: 6.0,2.2,5.0,1.5(CorrectClass=1,Iris-virginica) -> " << calculate.predict(test2) << endl;
   
   /* system("PAUSE");
    return 0;*/
}

