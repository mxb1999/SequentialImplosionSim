#include <vector>
#include <string>
#include <iostream>



const int nx=200; const float xmin = -5.0e-4; const float xmax=5.0e-4; const float dx = (xmax-xmin)/(nx-1);
const int nz=200; const float zmin = -5.0e-4; const float zmax=5.0e-4; const float dz = (zmax-zmin)/(nz-1);
const int nbeams = 2;
const int rays_per_zone = 5 ;
const double c =29979245800.0;
const double sigma = 1.7e-4;
const double e0 =8.85418782e-12;
const double me =9.10938356e-31;
const double pi =3.14159265359;
const double ec = 1.60217662e-19;
const double lambda = 1.053e-4/3.0;	// wavelength of light, in cm. This is frequency-tripled "3w" or "blue" (UV) light

const double freq = c/lambda;		// frequency of light, in Hz
const double omega = 2*pi*freq;	// frequency of light, in rad/s
const double ncrit = 1e-6*(pow(omega,2.0)*me*e0/pow(ec,2.0));
const double beam_max_z = 3.0e-4; const double beam_min_z = -3.0e-4;
const int nrays=int(rays_per_zone*(beam_max_z-beam_min_z)/dz)+0;
const float courant_mult = 0.2; // 0.37 // 0.25 // 0.36 // 0.22;
const double dt=courant_mult*min(dx,dz)/c;
const int nt=int(pow(courant_mult,-1.0)*max(nx,nz)*2.0);
const int numstored = int(5*rays_per_zone);
/*Defining Constants*/
