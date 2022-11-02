#include <iostream>
#include <unordered_map>
#include <algorithm>
 
class myclass
{
  template<typename K, typename V>
  void operator()(const std::pair<K, V> &p) {
    std::cout << "{" << p.first << ": " << p.second << "}\n";
  }
} ob;
 
template<typename K, typename V>
void print(const std::pair<K, V> &p) {
  std::cout << "{" << p.first << ": " << p.second << "}\n";
}
 
template<typename K, typename V>
void print_map(std::unordered_map<K, V> const &m)
{
  // specify a lambda expression
  std::for_each(m.begin(),
                m.end(),
                [](const std::pair<int, char> &p) {
		  std::cout << "{" << p.first << ": " << p.second << "}\n";
                });
 
  // or pass an object of a class overloading the ()operator
  // std::for_each(m.begin(), m.end(), ob);
 
  // or specify a function
  // std::for_each(m.begin(), m.end(), print<K, V>);
}
 
int main()
{
  std::unordered_map<int, char> m = {
    {1, 'A'},
    {2, 'B'},
    {3, 'C'}
  };
 
  print_map(m);
 
  return 0;
}
