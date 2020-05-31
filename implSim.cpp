#include <iostream>
#include <omp.h>
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>
#include <cmath>
#include "customMath.h"
#include "implSim.h"
#include "arrCheck.h"
using namespace std;
void span(double* target, double start, double stop, int num)
{
  float increment = (stop-start)/(num - 1);
  for(int i = 0; i < num; i++)
  {
    *(target + i) = start + (increment * i);
  }
}
void initialize()
{
  uray = new double[nt];
  intersections = new double*[nx]; //nx nz
  marked = new int***[nx]; //nx nz numstored nbeams
  dedendx = new double*[nx]; //nx nz
  dedendz = new double*[nx]; //nx nz
  x = new double*[nx]; //nx nz
  z = new double*[nx]; //nx nz
  eden = new double*[nx]; //nx nz
  edep = new double**[nx+2]; //nx+2 nz+2 nbeams
  present = new int**[nx]; //nx nz nbeams
  machnum = new double*[nx]; //nx nz
  boxes = new int***[nbeams]; //nbeams nrays nx*3 2
  W1_storage = new double**[nx]; //nx nz numstored
  W2_storage = new double**[nx]; //nx nz numstored
  u_flow = new double*[nx]; //nx nz
  dkx = new double**[nbeams]; //nbeams nrays 2
  dkz = new double**[nbeams]; //nbeams nrays 2
  dkmag = new double**[nbeams]; //nbeams nrays 2
  W1 = new double*[nx];//nx nz
  W2 = new double*[nx];//nx nz
  W1_init = new double*[nx];//nx nz
  W2_init = new double*[nx];//nx nz
  W1_new = new double*[nx];//nx nz
  W2_new = new double*[nx];//nx nz
  i_b1 = new double*[nx+2];//nx+2 nz+2
  i_b2 = new double*[nx+2];//nx+2 nz+2
  wpe = new double*[nx]; //nx nz
  crossesz = new double**[nbeams]; //nbeams nrays ncrossings
  crossesx = new double**[nbeams]; //nbeams nrays ncrossings
  ints = new int**[nbeams]; //nbeams nrays ncrossings
  for(int i = 0; i < nx+2; i++)
  {
    if(i < nx)
    {
      intersections[i] = new double[nz];
      x[i] = new double[nz];
      z[i] = new double[nz];
      eden[i] = new double[nz];
      machnum[i] = new double[nz];
      dedendx[i] = new double[nz];
      dedendz[i] = new double[nz];
      present[i] = new int*[nz];
      marked[i] = new int**[nz];
      W1_storage[i] = new double*[nz];
      W2_storage[i] = new double*[nz];
      u_flow[i] = new double[nz];
      W1[i] = new double[nz];
      W2[i] = new double[nz];
      W1_init[i] = new double[nz];
      W2_init[i] = new double[nz];
      W1_new[i] = new double[nz];
      W2_new[i] = new double[nz];
      wpe[i] = new double[nz];
      marked[i] = new int**[nz];
    }
    i_b1[i] = new double[nz+2];
    i_b2[i] = new double[nz+2];
    edep[i] = new double*[nz+2];
    for(int j = 0; j < nz+2; j++)
    {
      if(j < nz && i < nx)
      {
        marked[i][j] = new int*[numstored];
        present[i][j] = new int[nbeams];
        W1_storage[i][j] = new double[numstored];
        W2_storage[i][j] = new double[numstored];
        for(int m = 0; m < numstored;m++)
        {
          marked[i][j][m] = new int[nbeams];
        }
      }
      edep[i][j] = new double[nbeams];

    }
  }
  for(int i = 0; i < nbeams; i++)
  {
    boxes[i] = new int**[nrays];
    dkx[i] = new double*[nrays];
    dkz[i] = new double*[nrays];
    dkmag[i] = new double*[nrays];
    crossesz[i] = new double*[nrays];
    crossesx[i] = new double*[nrays];
    ints[i] = new int*[nrays];
    for(int j = 0; j < nrays; j++)
    {
      dkx[i][j] = new double[2];
      dkz[i][j] = new double[2];
      dkmag[i][j] = new double[2];
      crossesz[i][j] = new double[ncrossings];
      crossesx[i][j] = new double[ncrossings];
      boxes[i][j] = new int*[nx*3];
      ints[i][j] = new int[ncrossings];
      for(int m = 0; m < nx*3;m++)
      {
        boxes[i][j][m] = new int[2];
      }
    }
  }
  cout << "Setting initial conditions for ray tracker..." <<endl;
  cout << "nrays per beam is"<< nrays <<endl;
  double zstore[nx];
  span(zstore, zmin, zmax, nx);

  for(int i = 0; i < nx; i++)
  {
    for(int j = 0; j < nz; j++)
    {
      z[i][j] = zstore[i];
      //printf("%s%f%s%d%s%d\n", "XStore: ", x[i][j], "||i: ", i, "||j: ", j);
    }
    span(x[i],xmin, xmax,nz);
  }


  for(int i = 0; i < nx;i++)
  {
    for(int j = 0; j < nz;j++)
    {
      eden[i][j] = max(0.0,((0.3*ncrit-0.1*ncrit)/(xmax-xmin))*(x[i][j]-xmin)+(0.1*ncrit));
      machnum[i][j] = max(0.0,(((-0.4)-(-2.4))/(xmax-xmin))*(x[i][j]-xmin)+(-2.4));
    }
  }
  printf("%s\n", "Initialize Check 2");

  for(int i = 0; i < nx; i++)
  {
    for(int j = 0; j < nz; j++)
    {
      dedendz[i][j] = 0;
      dedendx[i][j] = 0;
    }
  }

}



void cbet()
{
  double constant1 = (pow(estat,2.0))/(4*(1.0e3*me)*c*omega*kb*Te*(1+3*Ti/(Z*Te)));
  cout << "Calculating CBET gains" << endl;
  double cs = 1e2*sqrt(ec*(Z*Te_eV+3.0*Ti_eV)/mi_kg);	// acoustic wave speed, approx. 4e7 cm/s in this example


  for(int i = 0; i < nx;i++)
  {
    for(int j = 0; j < nz;j++)
    {
      *(i_b1[i]+j) = *(*(edep[i]+j)+1);
      *(i_b2[i]+j) = *(*(edep[i]+j)+2);
      u_flow[i][j] = *(machnum[i]+j)*cs;
    }
  }
  	// acoustic wave speed, approx. 4e7 cm/s in this example
    int ray1num;

  for(int i = 0; i < nbeams-1;i++)
  {
    for(int j = 0; j < nrays; j++)
    {
      for(int m = 0; m < ncrossings;m++)
      {
        if(*(*(*(boxes[i]+ j)+ m)+ 1) == 0 || *(*(*(boxes[i]+ j)+ m)+ 2) == 0 )
        {
          break;
        }
        //marked array has dimensions nx nz numstored nbeams
        int ix = *(*(*(boxes[i]+ j)+ m));
        int iz = *(*(*(boxes[i]+ j)+ m)+ 1);
          if(*(intersections[ix]+iz) != 0)
          {
            vector<int> nonzeros1;
            vector<int> nonzeros2;
            int numrays1 = 0;
            int numrays2 = 0;
            for(int s = 0; s < numstored;s++)
            {
              if(marked[ix][iz][s][0] != 0)
              {
                nonzeros1.push_back(s);
                numrays1++;
              }
              if(marked[ix][iz][s][1] != 0)
              {
                nonzeros2.push_back(s);
                numrays2++;
              }
            }

            int marker1[numrays1];
            int marker2[numrays2];
            int mark2Copy[numrays2];
            int max = fmax(numrays1, numrays2);
            for(int l = 0; l < max; l++)
            {

              if(l < numrays1)
              {
                marker1[l] = marked[ix][iz][nonzeros1[l]][0];
              }

              if(l < numrays2)
              {
                marker2[l] = marked[ix][iz][nonzeros2[l]][1];
                mark2Copy[l] = marked[ix][iz][nonzeros2[l]][1];
              }
            }

          for(int r = 0; r < numrays1;r++)
          {
            if(marker1[r] == j)
            {
              ray1num = r;
              break;
            }
          }
          for(int n = 0; n < numrays2; n++)
          {
            for(int q = 0; q <ncrossings; q++)
            {
              int ix2 = *(*(*(boxes[i+1]+ mark2Copy[n])+ q) + 1);
              int iz2 = *(*(*(boxes[i+1]+ mark2Copy[n])+ q) + 2);
              if ( ix == ix2 && iz == iz2 )
              {
                mark2Copy[n] = q;
                break;
              }
            }
          }
        int n2limit = fmin(*(*(present[ix]+iz)+1),numrays2);

        for ( int n2 = 1; n2 < n2limit; n2++)
        {
          double ne = *(eden[ix]+iz);
          double epsilon = 1.0-ne/ncrit;
          double kmag = (omega/c)*sqrt(epsilon);         // magnitude of wavevector

          double kx1 = kmag*(*(*(dkx[i]+j)+m))/(*(*(dkmag[i]+j)+m)+1.0e-10);
          double kx2 = kmag*(*(*(dkx[i+1] + marker2[n2]) + mark2Copy[n2]))/(*(*(dkmag[i+1] + marker2[n2]) + mark2Copy[n2])+1.0e-10);
          double kz1 = kmag*(*(*(dkz[i] + j) + m))/(*(*(dkmag[i] + j) + m)+1.0e-10);
          double kz2 = kmag*(*(*(dkz[i] + j) + m))/(*(*(dkmag[i+1] + marker2[n2]) + mark2Copy[n2])+1.0e-10);


          double kiaw = sqrt(pow(kx2-kx1,2.0)+pow(kz2-kz1,2.0));
        // magnitude of the difference between the two vectors

          double ws = kiaw*cs;            // acoustic frequency, cs is a constant
          double omega1= omega;  // laser frequency difference. To start, just zero.
          double omega2 = omega;
   		    double eta = ((omega2-omega1)-(kx2-kx1)*u_flow[ix][iz])/(ws+1.0e-10);


         double efield1 = sqrt(8.*pi*1.0e7*(*(i_b1[ix]+iz))/c);             // initial electric field of rayi;lnlni46
        // double efield2 = sqrt(8.*pi*1.0e7*(*(i_b2[ix]+iz))/c);             // initial electric field of ray

         double P = (pow((iaw),2)*eta)/(pow((pow(eta,2)-1.0),2)+pow((iaw),2)*pow(eta,2));         // From Russ's paper
         //double gain1 = constant1*pow(efield2,2)*(ne/ncrit)*(1/iaw)*P;               //L^-1 from Russ's paper
         double gain2 = constant1*pow(efield1,2)*(ne/ncrit)*(1/iaw)*P;
                  //L^-1 from Russ's paper

                          // new energy of crossing (PROBE) ray (beam 2)
        if ( *(*(dkmag[i+1]+marker2[n2])+mark2Copy[n2]) >= 1.0*dx )
        {
          *(W2_new[ix]+iz) = *(W2[ix]+iz)*exp(-1*(*(W1[ix]+iz))*(*(*(dkmag[i+1]+marker2[n2])+mark2Copy[n2]))*gain2/sqrt(epsilon));
          *(*(W2_storage[ix]+iz)+n2) = *(W2_new[ix]+iz);
        }

        *(W1_new[ix]+iz) = *(W1[ix]+iz)*exp(-1*(*(W2[ix]+iz))*(*(*(dkmag[i]+j)+m))*gain2/sqrt(epsilon));

  // ENFORCE Energy conservation:
  //	       	                W1_new(ix,iz) = W1(ix,iz)-(W2_new(ix,iz)-W2(ix,iz));

                              *(*(W1_storage[ix]+iz)+raynum) = *(W1_new[ix]+iz);
  //				W1(ix,iz) = W1_new(ix,iz)


          }
  			}
      }
        if ( j % 20 == 0 )
        {
          cout << "     ..."<<(int)(100.*(1.0-(double)(j/double(1*nrays))))<<"%  remaining..."<<endl;
        }
  }


}
//  double clo = 0.99;
  //double chi = 1.0;
  cout << "Updating intensities due to CBET gains" << endl;
  double i_b1_new[nx][nz];
  double i_b2_new[nx][nz];
  for(int i = 0;i < nx; i++)
  {
    for(int j = 0; j < nz; j++)
    {
      i_b1_new[i][j] = i_b1[i][j];
      i_b2_new[i][j] = i_b2[i][j];
    }
  }
  for(int b = 0; b < beam - 1; b++)
  {
    for(int r = 0; r < nrays; r++)
    {
      for(int m = 0; m < ncrossings; m++)
      {
        if(*(*(*(boxes[b] + r) + m) + 1) == 0 || *(*(*(boxes[b] + r) + m) + 2) == 0)
        {
          break;
        }
          int ix = *(*(*(boxes[b] + r) + m) + 1);
          int iz = *(*(*(boxes[b] + r) + m) + 2);
          if(*(intersections[ix] + iz) != 0)
          {

            vector<int> nonzeros1;
            vector<int> nonzeros2;
            for(int l = 0; l < numstored;l++)
              {
                if(marked[ix][iz][l][0] != 0)
                {
                  nonzeros1.push_back(l);
                }
                if(marked[ix][iz][l][1] != 0)
                {
                  nonzeros2.push_back(l);
                }
              }
              int numrays1 = nonzeros1.size();
              int numrays2 = nonzeros2.size();
             int marker1[numrays1];
             int marker2[numrays2];
             int r2Copy[numrays2];
             int c2Copy[numrays2];
             for(int l = 0; l < fmax(numrays1, numrays2);l++)
             {
               if(l < numrays1)
               {
                 marker1[l] = marked[ix][iz][nonzeros1[l]][0];
               }
               if(l < numrays2)
               {
                 marker2[l] = marked[ix][iz][nonzeros1[l]][1];
                 r2Copy[l] = marker2[l];
                 c2Copy[l] = marker2[l];
               }
               for(int r = 0; r < numrays1;r++)
               {
                 if(marker1[r] == r)
                 {
                   ray1num = r;
                   break;
                 }
               }
              for(int n = 0; n < numrays2;n++)
              {
                for(int l = 0; l < ncrossings;l++)
                {
                  double ix2 = *(*(*(boxes[b+1]+r2Copy[n])+ l)+ 1);
                  double iz2 = *(*(*(boxes[b+1]+r2Copy[n])+ l)+ 2);
                  if(ix == ix2 && iz == iz2)
                  {
                    c2Copy[n] = l;
                    break;
                  }
                }
              }
              double fractional_change_1 = (-1.0*(1.0 - (*(W1_new[ix] + iz)/(*(W1_init[ix] + iz))) * (*(i_b1[ix] + iz))));
              double fractional_change_2 = (-1.0*(1.0 - (*(W2_new[ix] + iz)/(*(W2_init[ix] + iz))) * (*(i_b2[ix] + iz))));
              i_b1_new[ix][iz] += fractional_change_1;
              i_b2_new[ix][iz] += fractional_change_1;
              double x_prev_1 = *(x[ix] + iz);
              double z_prev_1 = *(z[ix] + iz);
              double x_prev_2 = *(x[ix] + iz);
              double z_prev_2 = *(z[ix] + iz);
              for(int l = m+1; l < ncrossings; l++)
              {
                int ix_next_1 = *(*(*(boxes[1]+r)+l)+1);
                int iz_next_1 = *(*(*(boxes[1]+r)+l)+2);
                double x_curr_1 = *(x[ix_next_1] + iz_next_1);
                double z_curr_1 = *(z[ix_next_1] + iz_next_1);
                if(ix_next_1 == 0 || iz_next_1 == 0)
                {
                  break;
                }else
                {
                  if(x_curr_1 != x_prev_1 || z_curr_1 != z_prev_1)
                  {
                    i_b1_new[ix_next_1][iz_next_1] +=fractional_change_1 * (*(*(present[ix] + iz) + 1)/(*(*(present[ix_next_1] + iz_next_1) + 1)));
                  }
                  x_prev_1 = x_curr_1;
                  z_prev_1 = z_curr_1;
                }
              }
              int n2 = fmin(ray1num, numrays2);

              for(int l = c2Copy[n2] + 1; l < ncrossings; l++)
              {
                int ix_next_2 = *(*(*(boxes[2]+r2Copy[n2])+l)+1);
                int iz_next_2 = *(*(*(boxes[2]+r2Copy[n2])+l)+2);
                int x_curr_2 = *(x[ix_next_2] + iz_next_2);
                int z_curr_2 = *(z[ix_next_2] + iz_next_2);
                if(ix_next_2 == 0 || iz_next_2 == 0)
                {
                  break;
                }else
                {
                  if(x_curr_2 != x_prev_2 || z_curr_2 != z_prev_2)
                  {
                    i_b2_new[ix_next_2][iz_next_2] += fractional_change_2*(*(*(present[ix]+iz)+1))/(*(*(present[ix_next_2]+iz_next_2)+2));
                  }
                  x_prev_2 = x_curr_2;
                  z_prev_2 = x_prev_2;
                }
              }
          }
        }
    }
    if(r % 20 == 0)
    {
      cout << "     ..." << int(100.*(1.0-(double(r)/double(1*nrays)))) << "%  remaining..."<<endl;
    }
  }
}
}

//use two threads here
void launch_ray_XZ(double x_init, double z_init, double kx_init, double kz_init)
{
//  printf("%f\n", z_init);
  //ofstream valTrack("ValTrack.txt");
  //Launch_Ray_XZ Array Declaration
  double myx[nt] = {0};
  double mytime[nt] = {0};
  span(mytime, dt, nt*dt, nt);
  double myz[nt] = {0};
  double mykx[nt]={0};
  double mykz[nt]= {0};
  double myvx[nt]={0};
  double myvz[nt]= {0};
  double amplitude_norm[nt]= {0};
  double markingx[nt] = {0};
  double markingz[nt] = {0};
  //double xbounds[nt];
  //double zbounds[nt];
  //double xbounds_double[nt];
  //double zbounds_double[nt];
  double nuei[nt] = {0};
  //Initializing Arrays
  for(int i = 0; i < nt; i++)
  {
    nuei[i] = 1.0;
  }
  myx[0] = x_init;
  myz[0] = z_init;
  for(int i = 0;i < nx;i++)
  {
    if(myx[0]- x[i][0] <= (0.5 +1.0e-10)*dx && myx[0]-x[i][0]>=-1*(0.5+1e-10)*dx)
    {
      thisx_0=i;
      thisx_00=i;
      break;
    }
  }
  for(int i = 0;i < nz;i++)
  {
    if(myz[0]-z[0][i] <= (0.5 +1.0e-10)*dx && myz[0]-z[0][i]>=-1*(0.5+1e-10)*dz)
    {
      thisz_0=i;
      thisz_00=i;
      break;
    }
  }
  double k = sqrt((pow(omega,2.0)-pow(*(wpe[thisx_0]+thisz_0),2.0)/pow(c,2.0)));
  double knorm = sqrt(pow(kx_init,2.0)+pow(kz_init,2.0));
  mykx[0]=(kx_init/knorm)*k;			// Normalized value for the ray's initial k_x
  mykz[0]=(kz_init/knorm)*k;			// Normalized value for the ray's initial k_z
  myvx[0] = pow(c,2.0)*mykx[0]/omega;                   // v_group, group velocity (dw/dk) from D(k,w).
  myvz[0] =  pow(c,2.0)*mykz[0]/omega;
/*
  cout <<  "omega: " << omega << "\n";
  cout << "wpe[thisx_0][thisz_0]: " << wpe[thisx_0][thisz_0] << "\n";
  cout << "thisz_0: " << thisz_0 << "\n";
  cout <<  "thisx_0: " << thisx_0 << "\n";
  cout << "x_init: " << x_init << "\n";
  cout << "z_init: " << z_init << "\n";
  cout << "kx_init: " << kx_init << "\n";
  cout << "kz_init: " << kz_init << "\n";
  cout <<  "MyKX[0]: " << mykx[0] << "\n";
  cout << "MyKZ[0]: " << mykz[0] << "\n";
  cout << "MyVX[0]: " << myvx[0] << "\n";
  cout << "MyVZ[0]: " << myvz[0] << "\n";
*/
  markingx[0] = thisx_0;
  markingz[0] = thisz_0;
//  xbounds[0] = myx[0];
//  zbounds[0] = myz[0];
//__________Time Stepping__________
  int numcrossing = 1;
  for(int i = 1; i < nt;i++)
  {
    myvz[i] = myvz[i-1] - pow(c,2.0)/(2.0*ncrit)*(*(dedendz[thisx_0])+thisz_0)*dt;
    myvx[i] = myvx[i-1] - pow(c,2.0)/(2.0*ncrit)*(*(dedendx[thisx_0])+thisz_0)*dt;
    myx[i] = myx[i-1] + myvx[i]*dt;
    myz[i] = myz[i-1] + myvx[i]*dt;
    int search_index_x = 1;
    int search_index_z = 1;
        int thisx_m = fmax(0, thisx_0-search_index_x);
        int thisx_p = fmin(nx, thisx_0+search_index_x);
        int thisz_m = fmax(0, thisz_0-search_index_z);
        int thisz_p = fmin(nz, thisz_0+search_index_z);
        /*
        valTrack <<  "1. thisx_m: " << thisx_m << "\n";
        valTrack <<  "1. thisx_p: " << thisx_p << "\n";
        valTrack <<  "1. thisz_m: " << thisz_m << "\n";
        valTrack <<  "1. thisz_p: " << thisz_p << "\n";
        */
        for(int j = thisx_m; j < thisx_p;j++)
        {
          if ( myx[i] - *(x[j]) <= (0.5+1.0e-10)*dx && myx[i] - *(x[j]) >= -(0.5+1.0e-10)*dx ){
                          thisx = j;
                          break;
              }
        }
      //  valTrack << "thisx: " << thisx << "\n";


        for(int j = thisz_m; j < thisz_p; j++)
        {
          if ( myz[i] - *(z[j]) <= (0.5+1.0e-10)*dz && myz[i] - *(z[j]) >= -(0.5+1.0e-10)*dz ){
                          thisz = j;
                          break;
              }
        }
        double linez[2]={myz[i-1], myz[i]};
        double linex[2]={myx[i-1], myx[i]};
        int lastx = 10000;
        int lastz = 10000;
        thisx_m = fmax(0, thisx_0-search_index_x);
        thisx_p = fmin(nx, thisx_0+search_index_x);
        thisz_m = fmax(0, thisz_0-search_index_z);
        thisz_p = fmin(nz, thisz_0+search_index_z);
        /*
        valTrack <<  "2. thisx_m: " << thisx_m << "\n";
        valTrack <<  "2. thisx_p: " << thisx_p << "\n";
        valTrack <<  "2. thisz_m: " << thisz_m << "\n";
        valTrack <<  "2. thisz_p: " << thisz_p << "\n";
        */
        interpFunc* lineFunc = new interpFunc(linez, linex, 2);

        for(int j = thisx_m; j < thisx_p;j++)
        {
          double currx = *(x[j])-dx/2;
          if((myx[i] > currx && myx[i-1] <= currx) || (myx[i] < currx && myx[i-1] >= currx))
          {

            double crossx = lineFunc->getInterp(currx);
            if(abs(crossx-lastz)>1.0e-20)
            {
              *(*(ints[beam]+raynum)+numcrossing)=uray[i];
              *(*(crossesx[beam]+raynum)+numcrossing) = currx;
              *(*(crossesz[beam]+raynum)+numcrossing) = crossx;


              if(myx[i] <= xmax+dx/2&&myx[i] >=xmin-dx/2)
              {
                **(*(boxes[beam]+raynum)+numcrossing) = thisx;
                *(*(*(boxes[beam]+raynum)+numcrossing) + 1)= thisz;
              }

              lastx = currx;
              numcrossing += 1;
              break;
            }
          }
        }
        for(int j = thisz_m; j < thisz_p;j++)
        {
          double currz = *(x[j])-dz/2;
          if((myz[i] > currz && myz[i-1] <= currz) || (myz[i] < currz && myz[i-1] > currz))
          {
             double crossz = lineFunc->getInterp(currz);
            if(abs(crossz-lastx)>1.0e-20)
            {
              *(*(ints[beam]+raynum)+numcrossing)=uray[i];
              *(*(crossesz[beam]+raynum)+numcrossing) = currz;
              *(*(crossesz[beam]+raynum)+numcrossing) = crossz;
              if(myz[i] <= zmax+dz/2&&myz[i] >=zmin-dz/2)
              {
                **(*(boxes[beam]+raynum)+numcrossing) = thisx;
                *(*(*(boxes[beam]+raynum)+numcrossing) + 1)= thisz;
              }
              lastz = currz;
              numcrossing += 1;
              break;
            }
          }
        }
        thisx_0 = thisx;
        thisz_0 = thisz;
        markingx[i] = thisx;
        markingz[i] = thisz;
        if(markingx[i] != markingx[i-1] && markingz[i] != markingz[i-1])
        {
          ztarg = *(z[thisx]+ thisz) - (dz/2.0);
          if(myvz[i] < 0.0)
          {
            ztarg = *(z[thisx]+ thisz) + (dz/2.0);
          }
          slope = (myz[i] - myz[i-1])/(myx[i] - myx[i-1]+1.0e-10);
      //    xbounds[i] = xtarg;
      //    zbounds[i] = ztarg;
          xtarg = *(x[thisx]+thisz)-(dx/2.0);
          if(myvx[i] >= 0.0)
          {
            xtarg = *(x[thisx]+thisz)+(dx/2.0);
          }
          slope = (myx[i] - myx[i-1])/(myz[i] - myz[i-1]+1.0e-10);
          ztarg = myz[i-1]+(xtarg-myx[i-1])/slope;
        //  xbounds_double[i] = xtarg;
        //  zbounds_double[i] = ztarg;
          for(int j = 0; j < numstored;j++)
          {
            if(marked[thisx_0][thisz_0][j][beam] == 0)
            {
              marked[thisx_0][thisz_0][j][beam] = raynum;
              *(*(present[thisx] + thisz) +beam) += 1.0;
              break;
            }
          }

        }else if(markingz[i] != markingz[i-1] )
        {
          ztarg = *(z[thisx] + thisz)-dz/2.0;
          if(myvz[i] < 0.0)
          {
            ztarg = *(z[thisx] + thisz)+dz/2.0;
          }
          slope = (myz[i] - myz[i-1])/(myx[i] - myx[i-1]+1.0e-10);
          xtarg = myx[i]+(ztarg-myz[i-1])/slope;
      //    xbounds[i]=xtarg;
        //  zbounds[i]=ztarg;
          for(int j = 0; j < numstored;j++)
          {
            if(marked[thisx_0][thisz_0][j][beam] == 0)
            {
              marked[thisx_0][thisz_0][j][beam] = raynum;
              *(*(present[thisx] + thisz) +beam) += 1.0;
              break;
            }
          }
        }else
        {
          xtarg = *(x[thisx] + thisz)-dx/2.0;
          if(myvx[i] < 0.0)
          {
            xtarg = *(x[thisx] + thisz)+dx/2.0;
          }
          slope = (myx[i] - myx[i-1])/(myz[i] - myz[i-1]+1.0e-10);
          ztarg = myz[i]+(xtarg-myx[i-1])/slope;
        //  xbounds[i]=xtarg;
        //  zbounds[i]=ztarg;
          for(int j = 0; j < numstored;j++)
          {
            if(marked[thisx_0][thisz_0][j][beam] == 0)
            {
            marked[thisx_0][thisz_0][j][beam] = raynum;
              *(*(present[thisx] + thisz) +beam) += 1.0;
              break;
            }
          }
        }
        uray[i] = uray[i-1];
  	    int increment = uray[i];
        double xp = (myx[i] - x[thisx][thisz]+dx/2.0)/dx;
  	    double zp = (myz[i] - z[thisx][thisz]+dz/2.0)/dz;
        cout <<"XP and ZP:  " << xp << " || " << zp << endl;
        /*
        valTrack <<  "xp: " << xp << "\n";
        valTrack <<  "zp: " << zp << "\n";
        valTrack <<  "myx[i]: " << myx[i] << "\n";
        valTrack <<  "myz[i]: " << myz[i] << "\n";
        valTrack <<  "x[thisx][thisz]: " << x[thisx][thisz] << "\n";
        valTrack <<  "z[thisx][thisz]: " << z[thisx][thisz] << "\n";
        */
        if ( xp >= 0 && zp >= 0 ){
        double dl = zp;
        double dm = xp;
        double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
        double a2 = (1.0-dl)*dm;		// green	: (x+1, z  )
        double a3 = dl*(1.0-dm);		// yellow	: (x  , z+1)
        double a4 = dl*dm;			// red 		: (x+1, z+1)
        *(*(edep[thisx+1]+thisz)+1) += a1*increment;	// blue
        *(*(edep[thisx+2]+thisz)+1) += a2*increment;	// green
        *(*(edep[thisx+1]+thisz)+2) += a3*increment;	// yellow
        *(*(edep[thisx+2]+thisz)+2) += a4*increment;	// red
      } else if ( xp < 0 && zp >= 0 ){
        double dl = zp;
        double dm = abs(xp);		// because xp < 0
        double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
        double a2 = (1.0-dl)*dm;		// green	: (x-1, z  )
        double a3 = dl*(1.0-dm);		// yellow	: (x  , z+1)
        double a4 = dl*dm;			// red 		: (x-1, z+1)
        *(*(edep[thisx+1]+thisz)+1) += a1*increment;	// blue
        *(*(edep[thisx]+thisz)+1) += a2*increment;	// green
        *(*(edep[thisx+1]+thisz)+2) += a3*increment;	// yellow
        *(*(edep[thisx]+thisz)+2) += a4*increment;	// red
      } else if ( xp >= 0 && zp < 0 ){
        //printf("%s\n", "Blach");
        double dl = abs(zp);		// because zp < 0
        double dm = xp;
        double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
        double a2 = (1.0-dl)*dm;		// green	: (x+1, z  )
        double a3 = dl*(1.0-dm);		// yellow	: (x  , z-1)
        double a4 = dl*dm;			// red 		: (x+1, z-1)
        (*(*(edep[thisx+1]+thisz)+1)) += a1*increment;	// blue
        (*(*(edep[thisx+2]+thisz)+1)) += a2*increment;	// green
        *(*(edep[thisx+1]+thisz)) += a3*increment;	// yellow
        *(*(edep[thisx+2]+thisz)) += a4*increment;	// red
      } else if ( xp < 0 && zp < 0 ){
        double dl = abs(zp);		// because zp < 0
        double dm = abs(xp);		// because xp < 0
        double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
        double a2 = (1.0-dl)*dm;		// green	: (x-1, z  )
        double a3 = dl*(1.0-dm);		// yellow	: (x  , z-1)
        double a4 = dl*dm;			// red 		: (x-1, z-1)
        (*(*(edep[thisx+1]+thisz)+1)) += a1*increment;	// blue
        (*(*(edep[thisx]+thisz)+1)) += a2*increment;	// green
        *(*(edep[thisx+1]+thisz)) += a3*increment;	// yellow
        *(*(edep[thisx]+thisz)) += a4*increment;	// red
      } else {
        double store = *(*(edep[thisx])+thisz);
        **(edep[thisx]+thisz) = store + (nuei[i] * (*(eden[thisx]+thisz))/ncrit * uray[i-1]*dt);
        cout << "***** ERROR in interpolation of laser deposition to grid!! *****" << endl;
        break;
      }
      amplitude_norm[i] = (pow(omega,2.0)-pow(*(wpe[thisx_00]+thisz_00),2.0))/(pow(omega,2.0)-pow(pow(*(wpe[thisx]+thisz),2.0),(1./4.)));
      mytime[i] = dt*i; // Sets the time, equal to the time step multiplied by the number of time steps

      // This will cause the code to stop following the ray once it escapes the extent of the plasma:
      if ( (myx[i] < (xmin-(dx/2.0))) | (myx[i] > (xmax+(dx/2.0))))
      {          // the "|" means "or" (symbol above the return key)
        int finalt = i-1;
        double* amp_norm = new double[finalt];
        double* rayx = new double[finalt];
        double* rayz = new double[finalt];
        for(int j = 0; j < finalt;j++)
        {
          amp_norm[j] = amplitude_norm[j];
          rayx[j] = myx[j];
          rayz[j] = myz[j];
        }
        break;                  // "breaks" out of the i loop once the if condition is satisfied
      } else if ( (myz[i] < (zmin-(dz/2.0))) | (myz[i] > (zmax+(dz/2.0)))){   // the "|" means "or" (symbol above the return key)
        int finalt = i-1;
        double* amp_norm = new double[finalt];
        double* rayx = new double[finalt];
        double* rayz = new double[finalt];
        for(int j = 0; j < finalt;j++)
        {
          amp_norm[j] = amplitude_norm[j];
          rayx[j] = myx[j];
          rayz[j] = myz[j];
        }

    }
}
}
void launchRays()
{
  cout << "Tracking Rays" << endl;
  double x0[nrays];
  double z0[nrays];
  double kx0[nrays];
  double kz0[nrays];
  double phase_x[nrays];
  double pow_x[nrays];
  span(phase_x,beam_min_z, beam_max_z, nrays);
  for(int i = 0; i < nrays; i++)
  {
    kx0[i] = 1.0;
    kz0[i] = -0.1;
    pow_x[i] = exp(-1*pow(pow(phase_x[i]/sigma,2.0),4.0/2.0));
  }
  double finalts[nrays][nbeams];

  interpFunc* powVPhase = new interpFunc(phase_x, pow_x, nrays);
  int beam = 1;
  cout << "BEAMNUM is" << beam << endl;
  span(z0, beam_min_z, beam_max_z, nrays);
  for(int i = 0; i < nrays; i++)
  {
    x0[i] = xmin-(dt/courant_mult*c*0.5);
    z0[i] += offset - (dz/2) - (dt/courant_mult*c*0.5);

  }
  for(int i = 0; i < nrays;i++)
  {
    raynum = i;
    uray[i] = uray_mult*powVPhase->getInterp(z0[i]);
    injected += uray[i];
    launch_ray_XZ(x0[i],z0[i],kx0[i],kz0[i]);
    *(finalts[i] + beam) = finalt;
    cout << kx0[i] << " || " << kz0[i] << '\n';

  }
  beam = 2;
  cout << "BEAMNUM is" << beam << endl;
  for(int i = 0; i < nrays;i++)
  {
    z0[i] = zmin-(dt/courant_mult*c*0.5);
    x0[i]-= (dx/2+dt/courant_mult*c*0.5);
    kx0[i] = 0.0;
    kz0[i] = 1.0;
    raynum = i;
    uray[i] = uray_mult*powVPhase->getInterp(z0[i]);
    injected += uray[i];
    launch_ray_XZ(x0[i],z0[i],kx0[i],kz0[i]);
    *(finalts[i] + beam) = finalt;
  }

printf("%s\n", "Got past launching the rays");
int if_intersection_diagnostics = 0;
for(int i = 0; i < nx; i++)
{
  for(int j = 0; j < nz; j++)
  {
    intersections[i][j] = 0;
  }
}
printf("%s\n", "RayCheck 1");
printf("%d\n", marked[0][0][0][0]);
for(int i = 1; i < nx;i++)
{
  for(int j = 1; j < nz; j++)
  {
    for(int m = 0; m < numstored; m++)
    {
      if(if_intersection_diagnostics == 1)
      {
        cout << "i, j, and, m are now"<<i<<j<<m <<endl;
      }
      if(marked[i][j][m][1] == 0)
      {
        break;
      }else
      {
        for(int l = 0; l < numstored; l++)
        {
          if(marked[i][j][m][2]==0)
          {
            break;
          }else
          {
            *(intersections[i]+j) += 1.0;
          }
        }
      }
    }
  }
}
printf("%s\n", "RayCheck 2");
}

int main(int argc, char const *argv[]) {
  printf("%s\n", "Check 1");
  initialize();
  //printf("%s\n", "Check 2");
  launchRays();
//  printf("%s\n", "Check 3");
  cbet();
  return 0;
}
