#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
using namespace std;
class interpFunc
{
private:
  double* x;
  double* y;
  int n;
public:
  int binSearch(double tgt, double* in, int index, int size);
  double fitLine(double tgt, int index);
  double getInterp(double tgt);
  double* xData(){return x;};
  double* yData(){return y;};
  int getSize(){return n;};
  interpFunc(double* inX, double* inY, int size)
  {
    x = inX;
    y = inY;
    n = size;
  };
};
