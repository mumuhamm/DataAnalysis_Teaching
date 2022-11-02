#include<iostream>
#include<cmath>
using namespace std;
float f(float(x))
{
  return (sin(x)-cos(x));
}
float g(float(x))
{
  return (cos(x)+sin(x));
}

int main()
{
  float a,b,d,i,n;
  cout<<" Given Trancedental equation is sin(x)-cos(x) "<<endl;
  cout<<" Enter the initial guess of the root "<<endl;
  cin>>a;
  cout<<" Enter the number of Iterations "<<endl;
  cin>>i;
  if(f(a) != 0.0)
    {
      for(n=1;n<=i;n++)
	{
	  b=a-(f(a)/g(a));
	  if(f(b)==0)
	    {
	      cout<<" i = "<<n<<" Root of the given equation is "<<b<<endl;
	    }
	  else
	    {
	      d=fabs((b-a)/b);
	      cout<<" i = "<<n<<" a = "<<a<<" b = "<<b<<" |E| = "<<d*100<<" % "<<endl;
	      a=b;
	    }
	}
      cout<<" Root of the given equation is "<<b<<endl;
    }
  else
    {
      cout<<" Root of the given equation is "<<a<<endl;
    }
  return 0;
}
