#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include "customMath.h"
#include <limits>
//implementations of custom math functions
//binary search algorithm, searches given array in log(n) time
int binSearch(double tgt, double* in, int index, int size)
{
  double left = in[index];
  double mid = in[index + size/2];
  if(size == 0)
  {
    return index;
  }

  if(tgt > left && tgt < mid)
  {
    return binSearch(tgt, in, index,size/2);
  }
  if(tgt > left && tgt > mid)
  {
    return binSearch(tgt, in, index+size/2,size - 1);
  }

  return -1;
};

//returns the interpolated value from the given sorted arrays
//Note for Prof. Sefkow: This is the array performing differently from Yorick
//It functions identically to other implementations of interp functions online
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
  double m = (yArr[index+1]-yArr[index])/(xArr[index+1]-xArr[index]);
  double b = yArr[index]-m*xArr[index];
  return m*target+b;
};
//used for interpolating all values of a given array, relies on previous interp function
//xsize is the length of "xArr" and "yArr", size is the length of the target array
double* interpArr(double* xArr, double* yArr, double* target, int xsize, int size)
{
  double* result = new double[size];
  for(int i = 0; i < size; i++)
  {
    result[i] = interp(xArr, yArr, target[i],xsize);

  }
  return result;
}
//fills an array of size num with evenly spaced double values between start and stop
void span(double* target, double start, double stop, int num)
{
  float increment = (stop-start)/(num - 1);
  for(int i = 0; i < num; i++)
  {
    target[i] = start + (increment * i);
  }
}
