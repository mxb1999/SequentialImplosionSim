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
void updateH5();
//Values needed throughout simulation
inline double injected = 0.0;
inline double intensity = 2.0e15;
inline double uray_mult = intensity*(courant_mult)*pow(double(rays_per_zone),-1.0);
inline H5::H5File* store;
inline double timeKeep;
inline ofstream file;
inline int beam;
inline int raynum;
inline int thisx;
inline int thisz;
//Launch Ray Values
inline int thisx_0;
inline int xcnt;  // DELETE
inline int zcnt; // DELETE
inline int thisx_00;
inline int thisz_0;
inline int thisz_00;
inline int finalt;
inline double ztarg;
inline double slope;
inline double xtarg;
inline double k;
inline double knorm;
inline double cs;
inline int ray1num;
//Pointers for necessary arrays
inline double* uray; //nt
inline double* rayx;//nt
inline double* rayz;//nt
inline double* amp_norm;//nt
inline double** intersections; //nx nz
inline int**** marked; //nx nz numstored nbeams
inline double** dedendx; //nx nz
inline double** dedendz; //nx nz
inline double* x; //nx
inline double* z; //nz
inline double** eden; //nx nz
inline bool**** markedTrack; //nx nz numstored nbeams

inline double*** edep; //nx+2 nz+2 nbeams
inline int*** present; //nx nz nbeams
inline double** machnum; //nx nz
inline int**** boxes; //nbeams nrays nx*3 2
inline bool**** boxTrack;
inline double*** W1_storage; //nx nz numstored
inline double*** W2_storage; //nx nz numstored
inline double** u_flow; //nx nz
inline double*** dkx; //nbeams nrays 2
inline double*** dkz; //nbeams nrays 2
inline double*** dkmag; //nbeams nrays 2
inline double** W1;//nx nz
inline double** W2;//nx nz
inline double** W1_init;//nx nz
//Launch_Ray_XZ specific arrays (all have a length of nt)
inline double* myx; //nt
inline double* mytime;//nt
inline double* myz;//nt
inline double* mykx;//nt
inline double* mykz;//nt
inline double* myvx;//nt
inline double* myvz;//nt
inline double* amplitude_norm;//nt
inline double* markingx;//nt
inline double* markingz;//nt
inline double* nuei;//nt
//CBET specific arrays
inline double** W2_init;//nx nz
inline double** W1_new;//nx nz
inline double** W2_new;//nx n
inline double** i_b1;//nx nz
inline double** i_b2;//nx nz
inline double** i_b1_new;//nx nz
inline double** i_b2_new;//nx nz
inline double** wpe; //nx nz
inline double*** crossesz; //nbeams nrays ncrossings
inline double*** crossesx; //nbeams nrays ncrossings
inline int*** ints; //nbeams nrays ncrossings
//arrays used only for plotting
inline double** i_bplot;//nx nz
inline double** i_b_newplot;//nx nz
inline double** edenplot; //the array is eden/ncrit,  nx nz
inline double** edepplot; //nx nz
#endif
