#include <iostream>
#include <omp.h>
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <cmath>
#include "customMath.h"
#include "implSim.h"
using namespace std;

void span(double* target, double start, double stop, int num)
{
  float increment = (stop-start)/num;
  for(int i = 0; i <= num; i++)
  {
    *(target + i) = start + (increment*i);
  }
}
void initialize()
{
  cout << "Setting initial conditions for ray tracker..." <<endl;
  cout << "nrays per beam is"<< nrays <<endl;

  eden = new double*[nx];
  dedendx = new double*[nx];
  dedendz = new double*[nx];
  x = new double*[nx];
  z = new double*[nx];
  edep_x = new double*[nx+2];
  edep_z = new double*[nx+2];
  etemp = new double*[nx];
  machnum = new double*[nx];
  wpe = new double*[nx];
  eden_norm = new double*[nx];

  for(int i = 0; i < nx;i++)
  {
    for(int j = 0; j < nz;i++)
    {

        *(eden[i]+j) = fmax(0.0,((0.3*ncrit-0.1*ncrit)/(xmax-xmin))*(*(x[i]+j)-xmin)+(0.1*ncrit));
        *(etemp[i]+j) = 100.0;
        *(machnum[i]+j) = fmax(0.0,(((-0.4)-(-2.4))/(xmax-xmin))*(*(x[i]+j)-xmin))+(-2.4);
        *(wpe[i]+j) = sqrt(*(eden[i]+j)*1e6*pow(ec,2.0)/me*e0);
        *(eden_norm[i]+j) = *(eden[i]+j);
    }
  }
  for (int i=0;i<=nx-1;++i)
  {
	   for (int j=0;j<=nz-1;++j)
     {
       *(dedendz[i]+j) = (*(eden[i]+j+1)-*(eden[i]+j)-*(eden[i]+j)/(*(z[i]+j+1))-*(z[i]+j));
       *(dedendx[i]+j) = (*(eden[i]+j+1)-*(eden[i]+j))/(*(z[i]+j+1)-*(z[i]+j));
	   }
  }
  double storagex[nx];
  span(storagex,zmin, zmax, nx);
  dedendx[nx-1]=dedendx[nx-2];
  for(int i = 0; i < nx; i++)
  {
    for(int j = 0; j < nz; j++)
    {
        *(x[j]+i) = storagex[j];
    }
    span(z[i],zmin, zmax, nx);
  }
  for(int i = 0; i < nx; i++)
  {
    *(dedendz[i]+nz-1)=*(dedendz[i]+nz-2);
  }
  phase_x = new double[nrays];
  span(phase_x, beam_min_z,beam_max_z,nrays);
  pow_x = new double[nrays];
  for(int i = 0; i < nrays; i++)
  {
    pow_x[i] = exp(-1*pow(pow((phase_x[i]/sigma),2.0),(4.0/2.0)));
  }
  myx= new double[nt]; myz= new double[nt];
  mykx=new double[nt]; mykz=new double[nt];
  myvx=new double[nt]; myvz=new double[nt];
  uray=new double[nt]; nuei=new double[nt];
  amplitude_norm=new double[nt];
  for(int i = 0; i < nt;i++)
  {
    uray[i] = 1.0;
    nuei[i] = 1.0;
    amplitude_norm[i] = 1.0;
  }
   mytime=new double[nt];
   span(mytime,dt,nt*dt,nt);
   double injected = 0.0;

  finalts = new double*[nrays];//nrays, nbeams
  mysaved_x = new double**[nt];//nt, nrays, nbeams
  mysaved_z = new double**[nt];//nt , nrays, nbeams
  marked=new double***[nx];//nx, nz, numstored, nbeams
  present=new double**[nx];//nx nz nbeams
  markingx = new int[nt];
  markingz = new int[nt];
  xbounds = new double[nt];
  zbounds = new double[nt];
  xbounds_double = new double[nt];
  zbounds_double = new double[nt];
  int ncrossings = nx*3;	// Maximum number of potential grid crossings by a ray
  crossesx = new double**[nbeams];//nbeams nrays ncrossings
  crossesz = new double**[nbeams];//nbeams nrays ncrossings
  boxes = new double***[nbeams];//nbeams nrays ncrossings 2
  ints = new double**[nbeams];//nbeams nrays ncrossings
  x0 =  new double[nrays];
  z0 = new double[nrays];
  kx0 = new double[nrays];
  kz0 = new double[nrays];

  span(z0, beam_min_z, beam_max_z, nrays);
  for(int i = 0; i < nrays;i++)
  {
    kx0[i] = 1.0;
    kz0[i] = -0.1;
    x0[i] = xmin-(dt/courant_mult*c*0.5);
    z0[i] = z0[i]+offset-(dz/2)-(dt/courant_mult*c*0.5);
  }
}


void launchRays()
{
  cout << "Tracking Rays" << endl;
  int if_ray_tracker_diagnostic = 1;
  interpFunc powVPhase = new interpFunc(pow_x, phase_x, nrays);
  int beam = 1;
  cout << "BEAMNUM is" << beam << endl;
  for(int i = 0; i < nrays;i++)
  {
    int raynum = i;
    uray[i] = uray_mult*powVPhase.getInterp(z0[n]);
    injected += uray[i];

  }
}
void launch_ray_XZ(double x_init, double z_init, double kx_init, double kz_init)
{

}
