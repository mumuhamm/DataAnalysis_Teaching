#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>

using namespace std;

class xval{
  float y ;
  string name ;

private: 
  
public:
   xval(){y = 3.4; name = "look";}
  void sety(float& a){y = a;}

};
  int main(){
    xval x;
    float k = 4;
    x.sety(k);
    std::cout<< k<<"\n";
  }
