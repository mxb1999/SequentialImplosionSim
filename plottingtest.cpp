#include <iostream>
#include <fstream>
#include "./python3.6/matplotlib-cpp/matplotlibcpp.h"
namespace plt = matplotlibcpp;
int main(int argc, char const *argv[]) {
  std::vector<double> test1 = {1,2,3,4,5,6,7,8,9,10};
  std::vector<double> test2 = {1,22,33,44,55,66,77,88,99,100};
  plt::figure();
  plt::plot(test1,test2);
  plt::show();
  return 0;
}
//g++ -Wall -Werror -std=c++17 plottingtest.cpp -I/usr/include/python3.6 -I/usr/include/python3.6/numpy -lpython3.6m
// g++ -Wall -Werror -std=c++17 plottingtest.cpp -I/usr/include/python3.6 -I/usr/include/python3.6/numpy
//g++ -Wall -Werror -std=c++17 plottingtest.cpp -I/usr/include/python2.7
//g++ -Wall -Werror -std=c++17 plottingtest.cpp -I/usr/include/python3.6 -L//usr/lib/python3.6/config-3.6m-x86_64-linux-gnu/libpython3.6.so -lpython3.6
//g++ -Wall -Werror -std=c++17 plottingtest.cpp -I/usr/include/python3.6 -lpython3.7
