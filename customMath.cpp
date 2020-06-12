#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include "customMath.h"
#include <limits>

int binSearch(double tgt, double* in, int index, int size)
{
  double left = in[index];
  double mid = in[index + size/2];


  /*
  cout << endl;
  cout << "Size:" << size << endl;
  cout << "Left: " << left << endl;
  cout << "Target:" << tgt << endl;
  cout << "Middle: " << mid << endl;
  */
  if(size <= 2)
  {
    return index;
  }

  if(tgt > left && tgt < mid)
  {
      return binSearch(tgt, in, index,size/2);
  }
  if(tgt > left && tgt > mid)
  {
      return binSearch(tgt, in, index+size/2 -1,size/2);
  }

  return -1;
};
double interp(double* xArr, double* yArr, double target, int xSize)
{

  if(target <= xArr[0])
  {
      return yArr[0];
  }
  if(target >= xArr[xSize-1])
  {

      return yArr[xSize-1];
  }
  int index = binSearch(target, xArr, 0, xSize);
  //out << "index: " << index << endl;
  double m = (yArr[index+1]-yArr[index])/(xArr[index+1]-xArr[index]);
  double b = yArr[index]-m*xArr[index];
  //cout << "xArr[index]: " << xArr[index] << " || yArr[index]: "  << yArr[index] << endl;
  //cout << "xArr[index+1]: " << xArr[index+1] << " || yArr[index+1]: "  << yArr[index+1] << endl;
  return m*target+b;
};
bool areEqual(double a, double b, int precision)
{
  /*
  cout << scientific;
  cout << endl;
  cout << "fabs(a-b)"<< std::fabs(b-a)<<endl;
  cout << "eps*fabs(a+b)"<< 1e-10 * std::fabs(a+b)<<endl;
  */
  return (std::fabs(b-a) <= 1e-10 * std::fabs(a+b) * precision) || (std::fabs(a-b) <= std::numeric_limits<double>::min());
};
int compareDub(double a, double b)
{
  if(areEqual(a,b,1))
  {
    return 0;
  }
  if(a > b)
  {
    return 1;
  }
  return -1;
};


void span(double* target, double start, double stop, int num)
{
  float increment = (stop-start)/(num - 1);
  for(int i = 0; i < num; i++)
  {
    target[i] = start + (increment * i);
  }
}
