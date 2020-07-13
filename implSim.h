#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <cmath>
#include "H5Cpp.h"
#include "customMath.h"
#include "simConst.h"
#include "vector.h"
using namespace std;

#include "H5ArrayType.h"
#ifndef IMPLSIM_H_
#define IMPLSIM_H_
//note: I* and D* denote whether the vector is of type double or type int

    //Functions
  void initialize();
  void launchRays();
  void launch_ray_XZ(double x_init, double z_init, double kx_init, double kz_init, int raynum, double urayinit);
  void cbet();
  void updateH5();
  //Values needed throughout simulation

 extern double maxDev;
 //Launch Ray Values
 extern double cs;
 extern int ray1num;
 extern int beam;
 extern double maxDev;
 //Pointers for necessary arrays
 extern double* rayx;//nt
 extern double* rayz;//nt
 extern double* amp_norm;//nt
 extern double** intersections; //nx nz
 extern int** marked; //nbeams numstored nx nz
 extern double** dedendx; //nx nz
 extern double** dedendz; //nx nz
 extern double* x; //nx
 extern double* z; //nz
 extern double** eden; //nx nz

 extern double*** edep; //nx+2 nz+2 nbeams
 extern int*** present; //nx nz nbeams
 extern double** machnum; //nx nz
 extern int**** boxes; //nbeams nrays nx*3 2
 extern double*** W1_storage; //nx nz numstored
 extern double*** W2_storage; //nx nz numstored
 extern double** u_flow; //nx nz
 extern double*** dkx; //nbeams nrays 2
 extern double*** dkz; //nbeams nrays 2
 extern double*** dkmag; //nbeams nrays 2
 extern double** W1;//nx nz
 extern double** W2;//nx nz
 extern double** W1_init;//nx nz
 //Launch_Ray_XZ specific arrays (all have a length of nt)


 //CBET specific arrays
 extern double** W2_init;//nx nz
 extern double** W1_new;//nx nz
 extern double** W2_new;//nx n
 extern double*** i_b; //nbeams nx nz
 extern double*** i_b_new;//nbeams nx nz
 extern double** wpe; //nx nz
 extern int*** truemark;//nx nz nbeams nrays needs to store all of the ray intersections. Boxes was previously used in order to store the location, and marked a certain number of rays. To generalize the algorithm, another must be used.
 extern double*** crossesz; //nbeams nrays ncrossings
 extern double*** crossesx; //nbeams nrays ncrossings
 extern int*** ints; //nbeams nrays ncrossings
 //arrays used only for plotting
 extern double** i_bplot;//nx nz
 extern double** i_b_newplot;//nx nz
 extern double** edenplot; //the array is eden/ncrit,  nx nz
 extern double** edepplot; //nx nz
#endif
