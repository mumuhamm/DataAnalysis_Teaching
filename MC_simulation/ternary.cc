#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>


using namespace std;

int main(){
   
   srand((unsigned) time(0));
   int a, b ;
   for (int index = 0; index < 10; index++) {
      std::cout<<" print the generated random values of a and b"<<"\n";
      a = (rand() % 80) + 1;
      b = (rand() % 80) + 1;
      cout <<"a  :\t" << a <<" and   b  :\t" <<b<<"\n";
      for( unsigned k = (rand() % 5) + 1 ; k<=(a < b ? a : b) ; ++k){
          std::cout<<" print the value of the iterator "<< k <<"\n";
          }
      
   }
  //return(a < b ? a : b);

}
