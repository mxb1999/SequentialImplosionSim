#include "implSim.h"

using namespace std;
//main function called in program
int main(int argc, char const *argv[]) {
  /*
    Basic steps of the program:
    1. Initialize the Arrays
    2. Launch rays and track for overlapping beams
    3. Perform CBET calculations and update arrays
    4. Update/write HDF5 file with desired output arrays
    5. Plot arrays (performed by matplotting.py)
  */
  initialize();
  launchRays();
  cbet();
  updateH5();
  return 0;
}
