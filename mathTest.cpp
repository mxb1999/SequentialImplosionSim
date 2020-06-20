#include <iostream>
#include "customMath.h"
#include <cmath>

using namespace std;

int main(int argc, char const *argv[]) {
  double xarr[10] = {1,3,5,7,9,11,13,15,17,19};
//  double yarr[10] = {10,20,30,40,50,60,70,80,90,100};
  double tarr[10] = {2,4,6,8,10,12,14,16,18,20};
  for(int i = 0; i < 10; i++)
  {
    cout << binSearch(tarr[i], xarr, 0, 10) << endl;
  }
  return 0;
}
