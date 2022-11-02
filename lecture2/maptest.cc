#include <iostream>
#include <iterator>
#include <string>
#include <map>
  
using namespace std;
  
int main()
{
  
  // empty map container
  std::map<int, std::string> gquiz1;
   string str1;
   string str2;
   string str3;
   string str4;
   string str5;
   string str6;
   string str7;
  
  // insert elements in random order
  gquiz1.insert(pair<int, std::string>(1, "str1"));
  gquiz1.insert(pair<int, std::string>(2, "str2"));
  gquiz1.insert(pair<int, std::string>(3, "str3"));
  gquiz1.insert(pair<int, std::string>(4, "str4"));
  gquiz1.insert(pair<int, std::string>(5, "str5"));
  gquiz1.insert(pair<int, std::string>(6, "str6"));
  gquiz1.insert(pair<int, std::string>(7, "str7"));
  
  // printing map gquiz1
  map<int, std::string>::iterator itr;
  cout << "\nThe map gquiz1 is : \n";
  cout << "\tKEY\tELEMENT\n";
  for (itr = gquiz1.begin(); itr != gquiz1.end(); ++itr) {
    cout << '\t' << itr->first
	 << '\t' << itr->second.c_str() << '\n';
  }
  cout << endl;
  
  // assigning the elements from gquiz1 to gquiz2
  map<int, std::string> gquiz2(gquiz1.begin(), gquiz1.end());
  
  // print all elements of the map gquiz2
  cout << "\nThe map gquiz2 after"
       << " assign from gquiz1 is : \n";
  cout << "\tKEY\tELEMENT\n";
  for (itr = gquiz2.begin(); itr != gquiz2.end(); ++itr) {
    cout << '\t' << itr->first
	 << '\t' << itr->second << '\n';
  }
  cout << endl;
  
  // remove all elements up to
  // element with key=3 in gquiz2
  cout << "\ngquiz2 after removal of"
    " elements less than key=3 : \n";
  cout << "\tKEY\tELEMENT\n";
  gquiz2.erase(gquiz2.begin(), gquiz2.find(3));
  for (itr = gquiz2.begin(); itr != gquiz2.end(); ++itr) {
    cout << '\t' << itr->first
	 << '\t' << itr->second << '\n';
  }
  
  // remove all elements with key = 4
  int num;
  num = gquiz2.erase(4);
  cout << "\ngquiz2.erase(4) : ";
  cout << num << " removed \n";
  cout << "\tKEY\tELEMENT\n";
  for (itr = gquiz2.begin(); itr != gquiz2.end(); ++itr) {
    cout << '\t' << itr->first
	 << '\t' << itr->second << '\n';
  }
  
  cout << endl;
  
  // lower bound and upper bound for map gquiz1 key = 5
  cout << "gquiz1.lower_bound(5) : "
       << "\tKEY = ";
  cout << gquiz1.lower_bound(5)->first << '\t';
  cout << "\tELEMENT = "
       << gquiz1.lower_bound(5)->second << endl;
  cout << "gquiz1.upper_bound(5) : "
       << "\tKEY = ";
  cout << gquiz1.upper_bound(5)->first << '\t';
  cout << "\tELEMENT = "
       << gquiz1.upper_bound(5)->second << endl;
  
  return 0;
}
