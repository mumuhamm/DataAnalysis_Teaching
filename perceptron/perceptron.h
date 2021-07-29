#ifndef _PERCEPTRON_H
#define _PERCEPTRON_H
#pragma once

#include <iostream>
#include <vector>
#include <list>
#include <cstdlib>
#include <cmath>
#include "TH1D.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TMath.h"
#include <cstdio>
using namespace std;

class perceptron
{
private:
   double m_eta;
   int m_epochs;
   vector <double> m_w;
   vector <double> m_errors;
   
public:
   
   
   perceptron(double eta, int epochs);
   double netInput(vector<double> X);
   int predict(vector<double> X);
   void fit( vector<vector<double>>  X, vector<double> y);
   void printErrors();
   void exportWeights(string filename);
   void importWeights(string filename);
   void printWeights();

};
#include "perceptron.cpp"
#endif
