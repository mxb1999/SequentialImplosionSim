#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <bits/stdc++.h>

#include "customMath.h"

using namespace std;

int main(int argc, char const *argv[]) {
  double testerX[10];
  double testerY[10];
  for(int i = 0; i < 10; i++)
  {
    testerX[i] = (double)rand()/RAND_MAX * 10;
    testerY[i] = (double)rand()/RAND_MAX * 10;
  }
  sort(testerX, testerX+10);
  sort(testerY,testerX+10);
  cout<<"X Points: ";
  for(int i = 0; i < 10; i++)
  {
    cout << testerX[i];
    cout<<", ";
  }
  cout<<endl;
  cout<<"Y Points: ";
  for(int i = 0; i < 10; i++)
  {
    cout << testerY[i];
    cout<<", ";
  }
  cout<<endl;
  interpFunc testFunc = interpFunc(testerX, testerY, 10);
  cout << testFunc.getInterp((testerX[0] + testerX[9])/2) << '\n';
  return 0;
}
