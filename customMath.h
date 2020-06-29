#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
using namespace std;
//custom math functions used in simulation code
extern int binSearch(double tgt, double* in, int index, int size);
extern double interp(double* xArr, double* yArr, double target, int xSize);
extern double* interpArr(double* xArr, double* yArr, double* target, int xsize, int size);
extern bool areEqual(double a, double b);
extern int compareDub(double a, double b);
extern void span(double* target, double start, double stop, int num);
