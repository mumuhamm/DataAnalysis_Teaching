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
  //MyMatrix m(10,15);
  a = MyMatrix33f::Zero(); // fill matrix elements with zeros
  a = MyMatrix33f::Identity(); // fill matrix as Identity matrix
  v = MyVector3f::Random(); // fill matrix elements with random values
  std::cout<<"Get the output of the identity matrix \n"
	   <<a<<"\n";
  std::cout<<"Get the output of the random vector \n"
	   <<"["<<v<<"]\n";



  int data_x[] = {1,2,3,4};
  Eigen::Map<Eigen::RowVectorXi> vx(data_x,4);
  std::vector<float> data = {1,2,3,4,5,6,7,8,9};
  Eigen::Map<MyMatrix33f> ax(data.data());

   std::cout<<"Get the output of the row vector with four inputs\n"
   <<"["<<vx<<"]\n";
   std::cout<<"Get the output of the filled data in to the 33f case \n"
           <<"["<<ax<<"]\n";

  auto amat = Eigen::Matrix2d::Random();
  auto bmat = Eigen::Matrix2d::Random();
  auto result_add = amat + bmat;
   
 
  
  auto result_mul = amat.array() * bmat.array(); // element wise multiplication
  auto result_div = amat.array() / bmat.array();
  auto result_mat_mul = amat * bmat; // matrix multiplication
  auto bmat_change = bmat.array() * 4;
   
   
   std::cout<<result_add<<"\n";
   std::cout<<result_mul<<"\n";
   std::cout<<result_div<<"\n";
   std::cout<<result_mat_mul<<"\n";
   std::cout<< bmat_change<<"\n";
   
   
   Eigen::MatrixXf m(4,4);
   m <<  1, 2, 3, 4,
   5, 6, 7, 8,
   9,10,11,12,
   13,14,15,16;
   cout << m<<"\n";
   cout << "m.leftCols(2) =" << endl << m.leftCols(2) << endl << endl;
   cout << "m.bottomRows<2>() =" << endl << m.bottomRows<2>() << endl << endl;
   m.topLeftCorner(1,3) = m.bottomRightCorner(3,1).transpose();
   cout << "After assignment, m = " << endl << m << endl;
   cout << "=========" <<"\n";
   cout << "Block in the middle" <<"\n";
   cout << m.block<2,2>(1,1) << endl << endl;
   for (int i = 1; i <= 3; ++i)
      {
      cout << "Block of size " << i << "x" << i << endl;
      cout << m.block(0,0,i,i) << endl << endl;
      }
   m.block(1,1,2,2) *= 4; // change values in original matrix
   std::cout<<m<<"\n";
   
   
   
   Eigen::Array22f k;
   k << 1,2,
   3,4;
   Array44f l = Array44f::Constant(0.6);
   cout << "Here is the array a:" << endl << l << endl << endl;
   l.block<2,2>(1,1) = k;
   cout << "Here is now a with m copied into its central 2x2 block:" << endl << l << endl << endl;
   l.block(0,0,2,3) = l.block(2,1,2,3);
   cout << "Here is now a with bottom-right 2x3 block copied into top-left 2x3 block:" << endl << l << endl << endl;
   
   
   Eigen::ArrayXf g(6);
   g << 1, 2, 3, 4, 5, 6;
   cout << "g.head(3) =" << endl << g.head(3) << endl << endl;
   cout << "g.tail<3>() = " << endl << g.tail<3>() << endl << endl;
   g.segment(1,4) *= 2;
   cout << "after 'g.segment(1,4) *= 2', g =" << endl << g << endl;
   
  
}
