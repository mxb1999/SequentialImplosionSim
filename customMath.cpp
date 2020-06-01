#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include "customMath.h"

int binSearch(double tgt, double* in, int index, int size)
{
  double comp = in[index+size/2];
  double plus = in[index+size/2+1];
  double minus = in[index+size/2-1];
  if(tgt > comp)
  {
    if(tgt > plus)
    {
      return binSearch(tgt, in, index+size/2,size/2);
    }else
    {
      return index+size/2;
    }
  }
  if(tgt < comp)
  {
    if(tgt < minus)
    {
      return binSearch(tgt, in, index,size/2);
    }else
    {
      return index+size/2 - 1;
    }
  }
  return index+size/2;
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
  cout << "index: " << index << endl;
  double m = (yArr[index]-yArr[index+1])/(xArr[index]-xArr[index+1]);
  double b = yArr[index]-m*xArr[index];
  cout << "xArr[index]: " << xArr[index] << " || yArr[index]: "  << yArr[index] << endl;
  cout << "xArr[index+1]: " << xArr[index+1] << " || yArr[index+1]: "  << yArr[index+1] << endl;
  return m*target+b;
}
