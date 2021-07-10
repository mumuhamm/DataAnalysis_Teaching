#include <iostream>
#include <cmath>
#include <cstdlib>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Cholesky>

using namespace std;
//using namespace RooFit;
using namespace Eigen;

int main()
{

  typedef Eigen::Matrix<float, 3, 3> MyMatrix33f;
  typedef Eigen::Matrix<float, 3, 1> MyVector3f;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MyMatrix;
  MyMatrix33f a;
  MyVector3f v;
  MyMatrix m(10,15);
  a = MyMatrix33f::Zero(); // fill matrix elements with zeros
  a = MyMatrix33f::Identity(); // fill matrix as Identity matrix
  v = MyVector3f::Random(); // fill matrix elements with random values
  std::cout<<"Get the output of the identity matrix \n"
	   <<a<<"\n";
  std::cout<<"Get the output of the random vector \n"
	   <<"["<<v<<"]\n";



  int data_x[] = {1,2,3,4};
  Eigen::Map<Eigen::RowVectorXi> vx(data_x,4);

  std::cout<<"Get the output of the row vector with four inputs\n"
           <<"["<<vx<<"]\n";
  
  std::vector<float> data = {1,2,3,4,5,6,7,8,9};
  Eigen::Map<MyMatrix33f> ax(data.data());

  
  std::cout<<"Get the output of the filled data in to the 33f case \n"
           <<"["<<ax<<"]\n";
  
}
