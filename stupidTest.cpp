#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "stupidTest.h"
#include <chrono>
#include <ctime>
#include <ratio>
using namespace std::chrono;
int main(int argc, char const *argv[]) {
  high_resolution_clock::time_point t1 = high_resolution_clock::now();
  printf("%s\n", "Hello");
  for(int i = 0; i < z; i++)
  {
    for(int j = 0; j < e; j++)
    {
      for(int l = 0; l < q; l++)
      {
        for(int n = 0; n < 4; n++)
        {
          double* pointer = test.access(i,j,l,n);
          *pointer = 4.0;
        }
      }
    }
  }
  high_resolution_clock::time_point t2 = high_resolution_clock::now();
  duration<double> span1 = duration_cast<duration<double>>(t2-t1);
  int test2[z][e][q];
  for(int i = 0; i < z; i++)
  {
    for(int j = 0; j < e; j++)
    {
      for(int l = 0; l < q; l++)
      {
        test2[i][j][l] = 4;
      }
    }
  }
  high_resolution_clock::time_point t3 = high_resolution_clock::now();
  cout<< test2[3][2][6] << endl;
  high_resolution_clock::time_point t4 = high_resolution_clock::now();
  duration<double> span2 = duration_cast<duration<double>>(t4-t3);
  cout << "Vector took: " << span1.count() << " seconds." <<endl;
  cout << "Array took: " << span2.count() << " seconds." <<endl;
  return 0;
}
