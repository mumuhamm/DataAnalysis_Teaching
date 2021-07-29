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

perceptron::perceptron(double eta, int epochs)
{
   m_epochs = epochs;
   m_eta = eta;
}


void perceptron::fit(vector<vector<double>> X , vector<double> y)
{
   for (int i = 0; i < X[0].size() + 1; i++) // X[0].size() + 1 -> I am using +1 to add the bias term
      {
      m_w.push_back(0);
      }
   for (int i = 0; i < m_epochs; i++)
      {
      double errors = 0;
      for (int j = 0; j < X.size(); j++)
         {
         double update = m_eta * (y[j] - predict(X[j]));
         for (int w = 1; w < m_w.size(); w++)
            {
               m_w[w] += update * X[j][w - 1];
            }
         m_w[0] = update;
         errors += update != 0 ? 1 : 0;
         }
      m_errors.push_back(errors);
      }
}




double perceptron::netInput(vector<double> X)
{
      // Sum(Vector of weights * Input vector) + bias
   double probabilities = m_w[0];
   for (int i = 0; i < X.size(); i++)
      {
      probabilities += X[i] * m_w[i + 1];
      }
   return probabilities;
}


int perceptron::predict(vector<double> X)
{
   return netInput(X) > 0 ? 1 : -1; //Step Function
}

void perceptron::printErrors()
{
   cout << "error: ";
   for (int i = 0; i < m_errors.size(); i++){
   std::cout<<m_errors[i]<<"\n";
   }
   
}

 
void perceptron::exportWeights(string filename)
{
   ofstream outFile;
   outFile.open(filename);
   
   for (int i = 0; i < m_w.size(); i++)
      {
      outFile << m_w[i] << endl;
      }
   
   outFile.close();
}

void perceptron::importWeights(string filename)
{
   ifstream inFile;
   inFile.open(filename);
   
   for (int i = 0; i < m_w.size(); i++)
      {
      inFile >> m_w[i];
      }
}

void perceptron::printWeights()
{
   cout << "weights: ";
   for (int i = 0; i < m_w.size(); i++)
      {
      cout << m_w[i] << " ";
      }
   cout << endl;
}



