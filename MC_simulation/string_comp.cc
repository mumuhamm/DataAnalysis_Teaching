#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <dirent.h>


using namespace std;
void comp(const char* lhs, const char* rhs)
{
   int rc = strcmp(lhs, rhs);
   const char *rel = rc < 0 ? "precedes" : rc > 0 ? "follows" : "equals";
   printf("[%s] %s [%s]\n", lhs, rel, rhs);
}

int main()
{

   DIR *dir;
   struct dirent *diread;
   vector<char*> files;
   if ((dir = opendir("/")) != nullptr) {
     while ((diread = readdir(dir)) != nullptr) {
        unsigned len = strlen(diread->d_name);
         std::cout<<len<<"\n";
         files.push_back(diread->d_name);
      }
      closedir (dir);
   }
   for (auto bara : files) cout << bara << " |==| ";
   cout << endl;
   
  /*  else {
      perror ("opendir");
      return EXIT_FAILURE;
   }
   
   

   
   const char* libus = "liliput";
   comp(libus, "muskan");
   comp(libus, "amraj");
   comp(libus, "lili");
   char key[]= "apple";
   char buffer[80];
   do{
      printf("Guess my favourite fruite?");
      fflush(stdout);
      scanf("%79s",buffer);
   }
   while (strcmp(key,buffer) !=0);
   puts("Correct Answer!");
   return 0;*/
}
