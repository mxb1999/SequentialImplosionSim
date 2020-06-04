#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include "simConst.h"


//Functions
void initialize();
void launchRays();
void launch_ray_XZ();
void cbet();

//Timing Arrays
double elapsed[3];
double elapsed0[3];
double elapsedtotal[3];
double cat01[3];
double cat02[3];
double cat03[3];
double cat04[3];
double cat05[3];
double cat06[3];
double cat07[3];
double cat08[3];
double cat09[3];
double cat10[3];
double cat11[3];
double cat12[3];
double cat13[3];
double cat14[3];
double cat15[3];
double cat16[3];
double cat17[3];
double cat18[3];
double cat19[3];
double cat20[3];
//Values needed throughout simulation
double injected = 0.0;
double intensity = 2.0e15;
double uray_mult = intensity*(courant_mult)*pow(double(rays_per_zone),-1.0);

double timeKeep;
int beam;
int raynum;
int ncrossings = nx*3;
int thisx;
int thisz;
//Launch Ray Values
int thisx_0;
int thisx_00;
int thisz_0;
int thisz_00;
int finalt;
double ztarg;
double slope;
double xtarg;
//Pointers for necessary arrays used across multiple functions
double* uray;
double* rayx;
double* rayz;
double* amp_norm;
double** intersections; //nx nz
int**** marked; //nx nz numstored nbeams
double** dedendx; //nx nz
double** dedendz; //nx nz
double** x; //nx nz
double** z; //nx nz
double** eden; //nx nz
double*** edep; //nx+2 nz+2 nbeams
int*** present; //nx nz nbeams
double** machnum; //nx nz
int**** boxes; //nbeams nrays nx*3 2
//___________________________________________
double*** W1_storage; //nx nz numstored
double*** W2_storage; //nx nz numstored
double** u_flow; //nx nz
double*** dkx; //nbeams nrays 2
double*** dkz; //nbeams nrays 2
double*** dkmag; //nbeams nrays 2
double** W1;//nx nz
double** W2;//nx nz
double** W1_init;//nx nz
double** W2_init;//nx nz
double** W1_new;//nx nz
double** W2_new;//nx nz
double** i_b1;//nx+2 nz+2
double** i_b2;//nx+2 nz+2
double** wpe; //nx nz
double*** crossesz; //nbeams nrays ncrossings
double*** crossesx; //nbeams nrays ncrossings
int*** ints; //nbeams nrays ncrossings
