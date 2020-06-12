#include "implSim.h"

using namespace std;

int main(int argc, char const *argv[]) {
  printf("%s\n", "Check 1");
  initialize();
  //printf("%s\n", "Check 2");
  launchRays();
//  printf("%s\n", "Check 3");
  cbet();
  
  updateH5();
  /*
  for(int i = 0; i < nx; i++)
  {
    for(int j = 0; j < nz;j++)
    {
      cout<< "edep[i][j][0]: " << edep[i][j][0]<<endl;
      cout<< "edep[i][j][1]: " << edep[i][j][1]<<endl;
    }
  }
  */
  return 0;
}
