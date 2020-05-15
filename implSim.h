#include <iostream>
#include <omp.h>
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <cmath>
#include "simConst.h"


//Functions
void initialize();
void launchRays();


//Timing Arrays
double elapsed[3] ;
double elapsed0[3] ;
double elapsedtotal[3] ;
double cat01[3] ;
double cat02[3] ;
double cat03[3] ;
double cat04[3] ;
double cat05[3] ;
double cat06[3] ;
double cat07[3] ;
double cat08[3] ;
double cat09[3] ;
double cat10[3] ;
double cat11[3] ;
double cat12[3] ;
double cat13[3] ;
double cat14[3] ;
double cat15[3] ;
double cat16[3] ;
double cat17[3] ;
double cat18[3] ;
double cat19[3];
double cat20[3];
//Values needed throughout simulation
double injected = 0.0;
double intensity = 2.0e15;
double uray_mult = intensity*(courant_mult)*pow(double(rays_per_zone),-1.0);

//Pointers for necessary arrays used across multiple functions
double** eden;
double** dedendz;
double** dedendx;
double** x;
double** z;
double** edep_x;
double** edep_z;
double** edep;
double** etemp;
double** wpe;
double** eden_norm;
double** machnum;
double* phase_x;
double* pow_x;
double* myx;
double* myz;
double* mykx;
double* mykz;
double* myvx;
double* myvz;
double* uray;
double* amplitude_norm;
double* nuei;
double* mytime;
double** finalts;
double*** mysaved_x;
double*** mysaved_z;
double**** marked;
double*** present;
int* markingx;
double* xbounds;
double* xbounds_double;
int* markingz;
double* zbounds;
double* zbounds_double;
double*** crossesx;
double*** crossesz;
double**** boxes;
double*** ints;
double* x0;
double* z0;
double* kx0;
double* kz0;
double** intersections;
double*** W1_storage;
double*** W2_storage;
