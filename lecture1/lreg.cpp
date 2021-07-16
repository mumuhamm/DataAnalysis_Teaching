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

  typedef float DType;
  using Matrix = Eigen::Matrix<DType, Eigen::Dynamic, Eigen::Dynamic> ;
  typedef Eigen::Matrix<float, 10, 10> MyMatrix10f;
   MyMatrix10f x;
   MyMatrix10f y;
   x = MyMatrix10f::Random();
   y = MyMatrix10f::Random();
   
  int n = 100000;
  //Matrix x(n,1);
  //Matrix y(n,1);
  Eigen::LeastSquaresConjugateGradient<Matrix> gd;
  gd.setMaxIterations(10000);
  gd.setTolerance(0.001) ;
  gd.compute(x);
  auto b = gd.solve(y);
  auto fu = b.transpose().array();
   std::cout<<"coefficient b as it get solved "<< b <<"\n";
   std::cout<<"coefficient fu as it get solved "<< fu <<"\n";


  /* Eigen::MatrixXf new_x(5, 2);
   new_x << 1, 1, 1, 2, 1, 3, 1, 4, 1, 5;
   auto new_y = new_x.array().rowwise() * b.transpose().array();
*/

}
