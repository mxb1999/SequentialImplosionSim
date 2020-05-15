#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include "customMath.h"

int interpFunc::binSearch(double tgt, double* in, int index, int size)
{
  if(tgt > *(in+size-1) || tgt < *(in))
  {
    std::cout << "-1 Returned" << '\n';
    return -1;
  }
  double comp = *(in+index+size/2);
  double plus = *(in+index+size/2+1);
  double minus = *(in+index+size/2-1);
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
double interpFunc::fitLine(double tgt, int index)
{
  double m = (*(y+index+1)-*(y+index))/(*(x+index+1)-*(x+index));
  double b = *(y+index)-m*(*(x+index));
  return m*tgt + b;
};
double interpFunc::getInterp(double tgt)
{
  int index = binSearch(tgt, this->x, 0, this->n);
  if(index == -1)
  {
    std::cout << "Value not in range of function" << '\n';
    return -1.0
  }
  return fitLine(tgt, index);
};
