/*The xtensor library is a C++ library for numerical analysis with multidimensional array expressions. Containers of xtensor are inspired by NumPy, the Python array programming library. ML algorithms are mainly described using Python and NumPy, so this library can make it easier to move them to C++. The following container classes implement multidimensional arrays in the xtensor library.*/
//Muhammad Alibordi

#include <iostream>
#include "/usr/local/opt/xtensor/include/xtensor/xarray.hpp"
#include "/usr/local/opt/xtensor/include/xtensor/xio.hpp"
#include "/usr/local/opt/xtensor/include/xtensor/xview.hpp"
#include "/usr/local/opt/xtensor/include/xtensor/xstrided_view.hpp"
#include "/usr/local/opt/xtensor/include/xtensor/xfixed.hpp"
#include "/usr/local/opt/xtensor/include/xtensor/xadapt.hpp"
#include "/usr/local/opt/xtensor/include/xtensor/xrandom.hpp"
using namespace std;

int main(int argc, char* argv[])
{
  xt::xarray<double> arr1
  {{1.0, 2.0, 3.0},
      {2.0, 5.0, 7.0},
	{2.0, 5.0, 7.0}};
  xt::xarray<double> arr2
  {5.0, 6.0, 7.0};
  xt::xarray<double> res = xt::view(arr1, 1) + arr2;
  std::cout << res << std::endl;
   
   std::vector<size_t> shape = { 3, 2, 4 };
   xt::xarray<double, xt::layout_type::row_major> a(shape);
   std::cout << a << std::endl;
   
   std::array<size_t, 3> shapex = { 3, 2, 4 };
   xt::xtensor<double, 3> ax(shapex);
   std::cout << ax << std::endl;
   
   xt::xtensor_fixed<double, xt::xshape<3, 2, 4>> afix;
   std::cout << afix << std::endl;
   
   auto e = xt::ones<double>({2, 3});
   std::cout << e << std::endl;
   
   // Evaluated versions--ones
   using fixed_tensor = xt::xtensor_fixed<double, xt::xshape<2, 3>>;
   xt::xarray<double>     O0 = xt::ones<double>({2, 3});
   xt::xtensor<double, 2> O1 = xt::ones<double>({2, 3});
   fixed_tensor           O2 = xt::ones<double>({2, 3});
   std::cout << O1<<O2<<O0 << std::endl;
   
   // Evaluated versions -- zeros
   //using fixed_tensor = xt::xtensor_fixed<double, xt::xshape<2, 3>>;
   xt::xarray<double>     Z0 = xt::zeros<double>({2, 3});
   xt::xtensor<double, 2> Z1 = xt::zeros<double>({2, 3});
   fixed_tensor           Z2 = xt::zeros<double>({2, 3});
   std::cout <<"----------------------------"<< std::endl;
   std::cout <<Z0<<Z1<<Z2<< std::endl;
   
   // Evaluated versions -- Empty
   xt::xarray<double>::shape_type sh0 = {2, 3};
   auto E0 = xt::empty<double>(sh0); // E0 is xt::xarray<double>
   xt::xtensor<double, 2>::shape_type sh1 = {2, 3};
   auto E1 = xt::empty<double>(sh1); // E1 is xt::xtensor<double, 2>
   xt::xshape<2, 3> sh2;
   auto E2 = xt::empty<double>(sh2); // E2 is xt::xtensor_fixed<double, xt::xshape<2, 3>>
   std::cout <<"----------------------------"<< std::endl;
   std::cout <<E0<<E1<<E2<< std::endl;
   
   
   std::vector<float> data{1,2,3,4};
   std::vector<size_t> shapek{2,2};
   auto data_x = xt::adapt(data, shapek);
   std::cout <<"----------------------------"<< std::endl;
   std::cout <<data_x<< std::endl;
   
   auto aarr = xt::random::rand<double>({2,2});
   auto barr = xt::random::rand<double>({2,2});
   auto carr = aarr + barr;
   std::cout <<"----------------------------"<< std::endl;
   std::cout <<carr<< std::endl;
  return 0;
}
