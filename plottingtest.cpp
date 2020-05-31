#include <iostream>
#include <fstream>
#include "./python3.6/matplotlib-cpp/matplotlibcpp.h"
namespace plt = matplotlibcpp;
int main(int argc, char const *argv[]) {
  std::vector<double> test = {1,2,3,4,5,6,7,8,9,10};
  plt::plot(test);
  plt::show();
  return 0;
}
