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

#include "H5ArrayType.h"
#ifndef IMPLSIM_H_
#define IMPLSIM_H_
    //Functions
  void initialize();
  void launchRays();
  void launch_ray_XZ(double x_init, double z_init, double kx_init, double kz_init);
  void cbet();
  int markedAccess(int dim1, int dim2, int dim3, int dim4);
  void markedWrite(int dim1, int dim2, int dim3, int dim4, int val);
  void updateH5();
  //Values needed throughout simulation
  template<typename G>
  std::vector<std::vector<std::vector<G>>> vector3d(int a, int b, int c, G val = G{})
  {
    return std::vector<std::vector<std::vector<G>>>(a, std::vector<std::vector<G>>(b, std::vector<G>(c, val)));
  }
  extern int beam;
  extern int raynum;
  extern int thisx;
  extern int thisz;
  //Launch Ray Values
  extern int thisx_0;
  extern int thisx_00;
  extern int thisz_0;
  extern int thisz_00;
  extern int finalt;
  extern double ztarg;
  extern double slope;
  extern double xtarg;
  extern double k;
  extern double knorm;
  extern double cs;
  extern int ray1num;
  //Pointers for necessary arrays
  extern double* uray; //nt
  extern double* rayx;//nt
  extern double* rayz;//nt
  extern double* amp_norm;//nt
  extern double** intersections; //nx nz
  extern int** marked; //nx nz numstored nbeams
  extern double** dedendx; //nx nz
  extern double** dedendz; //nx nz
  extern double* x; //nx
  extern double* z; //nz
  extern double** eden; //nx nz

  extern double*** edep; //nx+2 nz+2 nbeams
  extern int*** present; //nx nz nbeams
  extern double** machnum; //nx nz
  extern int**** boxes; //nbeams nrays nx*3 2
  extern bool**** boxTrack;
  extern double*** W1_storage; //nx nz nrays
  extern double*** W2_storage; //nx nz nrays
  extern double** u_flow; //nx nz
  extern double*** dkx; //nbeams nrays 2
  extern double*** dkz; //nbeams nrays 2
  extern double*** dkmag; //nbeams nrays 2
  extern double** W1;//nx nz
  extern double** W2;//nx nz
  extern double** W1_init;//nx nz
  //Launch_Ray_XZ specific arrays (all have a length of nt)
  extern double* myx; //nt
  extern double* mytime;//nt
  extern double* myz;//nt
  extern double* mykx;//nt
  extern double* mykz;//nt
  extern double* myvx;//nt
  extern double* myvz;//nt
  extern double* amplitude_norm;//nt
  extern double* markingx;//nt
  extern double* markingz;//nt
  extern double* nuei;//nt
  //CBET specific arrays
  extern double** W2_init;//nx nz
  extern double** W1_new;//nx nz
  extern double** W2_new;//nx n
  extern double** i_b1;//nx nz
  extern double** i_b2;//nx nz
  extern double** i_b1_new;//nx nz
  extern double** i_b2_new;//nx nz
  extern double** wpe; //nx nz
  extern double*** crossesz; //nbeams nrays ncrossings
  extern double*** crossesx; //nbeams nrays ncrossings
  extern int*** ints; //nbeams nrays ncrossings
  //arrays used only for plotting
  extern double** i_bplot;//nx nz
  extern double** i_b_newplot;//nx nz
  extern double** edenplot; //the array is eden/ncrit,  nx nz
  extern double** edepplot; //nx nz
#endif
