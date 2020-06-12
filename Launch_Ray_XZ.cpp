#include "implSim.h"
using namespace std;

void initializeArrXZ(double x_init, double z_init, double kx_init, double kz_init)
{
  //ofstream valTrack("ValTrack.txt");
  //Launch_Ray_XZ Array Declaration
  myx= new double[nt]{0.0};
  mytime = new double[nt]{0.0};
  span(mytime, dt, nt*dt, nt);
  myz = new double[nt]{0.0};
  mykx=new double[nt]{0.0};
  mykz= new double[nt]{0.0};
  myvx=new double[nt]{0.0};
  myvz= new double[nt]{0.0};
  amplitude_norm= new double[nt]{0.0};
  markingx = new double[nt]{0.0};
  markingz = new double[nt]{0.0};
  //double xbounds[nt];
  //double zbounds[nt];
  //double xbounds_double[nt];
  //double zbounds_double[nt];
  nuei = new double[nt]{0.0};
  //Initializing Arrays
  for(int i = 0; i < nt; i++)
  {
    nuei[i] = 1.0;
  }
  myx[0] = x_init;
  myz[0] = z_init;
  cout << scientific;
  for(int i = 0;i < nx;i++)
{
    if( myx[0] - x[i] <= (0.5+1.0e-10)*dx && myx[0] - x[i] >= -(0.5+1.0e-10)*dx )
  {
      thisx_0=i;
      //cout << "thisx_0: " << thisx_0<<endl;
      thisx_00=i;
      break;
    }
  }
  for(int i = 0;i < nz;i++)
  {
    if(myz[0] - z[i] <= (0.5+1.0e-10)*dz && myz[0] - z[i] >= -(0.5+1.0e-10)*dz )
    {
      thisz_0=i;
      thisz_00=i;
      break;
    }
  }
  k = sqrt((pow(omega,2.0)-pow(wpe[thisx_0][thisz_0],2.0))/pow(c,2.0));
  knorm = sqrt(pow(kx_init,2.0)+pow(kz_init,2.0));
  mykx[0]=(kx_init/knorm)*k;			// Normalized value for the ray's initial k_x
  mykz[0]=(kz_init/knorm)*k;			// Normalized value for the ray's initial k_z
  myvx[0] = pow(c,2.0)*mykx[0]/omega;                   // v_group, group velocity (dw/dk) from D(k,w).
  myvz[0] =  pow(c,2.0)*mykz[0]/omega;
  markingx[0] = thisx_0;
  markingz[0] = thisz_0;
//  xbounds[0] = myx[0];
//  zbounds[0] = myz[0];
}
void rayLaunch()
{
  //__________Time Stepping__________
    int numcrossing = 1;
    for(int i = 1; i < nt;i++)
    {
      myvz[i] = myvz[i-1] - pow(c,2.0)/(2.0*ncrit)*dedendz[thisx_0][thisz_0]*dt;
      myvx[i] = myvx[i-1] - pow(c,2.0)/(2.0*ncrit)*dedendx[thisx_0][thisz_0]*dt;
      myx[i] = myx[i-1] + myvx[i]*dt;
      myz[i] = myz[i-1] + myvz[i]*dt;
      /*
        cout << "k: " << k <<endl;
        cout << "knorm: " << knorm <<endl;
        cout << "pow(c,2.0)/(2.0*ncrit)*dedendx[thisx_0][thisz_0]*dt: " << pow(c,2.0)/(2.0*ncrit)*dedendx[thisx_0][thisz_0]*dt<<endl;
        cout << "dedendx[thisx_0][thisz_0]: " << dedendx[thisx_0][thisz_0] << endl;
        cout << "dedendz[thisx_0][thisz_0]: " << dedendz[thisx_0][thisz_0] << endl;
        cout << "mykx[i]: " << mykx[i] << endl;
        cout << "mykz[i]: " << mykz[i] << endl;
        cout << "myvx[i]: " << myvx[i] << endl;
        cout << "myx[i]: " << myx[i] << endl;
        cout << "myvz[i]: " << myvz[i] << endl;
        cout << "myz[i]: " << myz[i] << endl;
        cout << "myvx[0]: " << myvx[0] << endl;
        cout << "myx[0]: " << myx[0] << endl;
        cout << "myvz[0]: " << myvz[0] << endl;
        cout << "myz[0]: " << myz[0] << endl;
        cout << endl;
    */

      //cout <<"myx[i]: "<< myx[i] <<" ::myz[i]: "<< myz[i] <<" ::myvz[i]: "<< myvz[i] <<" ::myvx[i]: "<< myvx[i]<<endl;
      int search_index_x = nx;
      int search_index_z = nz;
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
        if ( myx[i] - x[j] <= (0.5+1.0e-10)*dx && myx[i] - x[j] >= -(0.5+1.0e-10)*dx )
        {
          thisx = j;
          break;
        }
      }
  //  valTrack << "thisx: " << thisx << "\n";


      for(int j = thisz_m; j < thisz_p; j++)
      {
        if ( myz[i] - z[j] <= (0.5+1.0e-10)*dz && myz[i] - z[j] >= -(0.5+1.0e-10)*dz )
        {
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
      for(int j = thisx_m; j < thisx_p;j++)
      {
        double currx = x[j]-dx/2;
        if((myx[i] > currx && myx[i-1] <= currx) || (myx[i] < currx && myx[i-1] >= currx))
        {
          double crossx =interp(linex, linez, currx, 2);
          if(abs(crossx-lastz)>1.0e-20)
          {
            *(*(ints[beam]+raynum)+numcrossing)=uray[i];
            *(*(crossesx[beam]+raynum)+numcrossing) = currx;
            *(*(crossesz[beam]+raynum)+numcrossing) = crossx;
            if(myx[i] <= xmax+dx/2 && myx[i] >= xmin-dx/2)
            {
              boxes[beam][raynum][numcrossing][0] = thisx;
              boxes[beam][raynum][numcrossing][1] = thisz;
              /*
              cout << "myx[i]: " << myx[i] << endl;
              cout << "xmax+dx/2: " << xmax+dx/2 << endl;
              cout << "xmin-dx/2: " << xmin-dx/2 << endl;
              cout << "i: " << i << endl;
              cout << "nt: " << nt << endl;
              cout << endl;
              */
            }
            lastx = currx;
            numcrossing += 1;
            break;
          }
        }
      }
          for(int j = thisz_m; j < thisz_p;j++)
          {
            double currz = x[j]-dz/2;
            if((myz[i] > currz && myz[i-1] <= currz) || (myz[i] < currz && myz[i-1] > currz))
            {
               double crossz = interp(linex, linez, currz, 2);
              if(abs(crossz-lastx)>1.0e-20)
              {
                ints[beam][raynum][numcrossing]=uray[i];
                crossesz[beam][raynum][numcrossing] = currz;
                crossesz[beam][raynum][numcrossing] = crossz;
                if(myz[i] <= zmax+dz/2&&myz[i] >=zmin-dz/2)
                {
                  boxes[beam][raynum][numcrossing][0] = thisx;
                  boxes[beam][raynum][numcrossing][1] = thisz;
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
            ztarg = z[thisz] - (dz/2.0);
            if(myvz[i] < 0.0)
            {
              ztarg = z[thisz] + (dz/2.0);
            }
            slope = (myz[i] - myz[i-1])/(myx[i] - myx[i-1]+1.0e-10);
        //    xbounds[i] = xtarg;
        //    zbounds[i] = ztarg;
            xtarg = x[thisx]-(dx/2.0);
            if(myvx[i] >= 0.0)
            {
              xtarg = x[thisx]+(dx/2.0);
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
                present[thisx][thisz][beam] += 1.0;
                break;
              }
            }

          }else if(markingz[i] != markingz[i-1] )
          {
            ztarg = z[thisz]-dz/2.0;
            if(myvz[i] < 0.0)
            {
              ztarg = z[thisz]+dz/2.0;
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
                present[thisx][thisz][beam] += 1.0;
                break;
              }
            }
          }else
          {
            xtarg = x[thisx]-dx/2.0;
            if(myvx[i] < 0.0)
            {
              xtarg = x[thisx]+dx/2.0;
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
                present[thisx][thisz][beam] += 1.0;
                break;
              }
            }
          }
          uray[i] = uray[i-1];
    	    double increment = uray[i];
          double xp = (myx[i] - (x[thisx]+dx/2.0))/dx;
    	    double zp = (myz[i] - (z[thisz]+dz/2.0))/dz;
          //cout << uray[0] << endl;
          //cout << i << endl;

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
          edep[thisx+1][thisz+1][beam] += a1*increment;	// blue
          edep[thisx][thisz+1][beam] += a2*increment;	// green
          edep[thisx+1][thisz+2][beam] += a3*increment;	// yellow
          edep[thisx][thisz+2][beam] += a4*increment;	// red
        } else if ( xp < 0 && zp >= 0 ){
          double dl = zp;
          double dm = abs(xp);		// because xp < 0
          double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
          double a2 = (1.0-dl)*dm;		// green	: (x-1, z  )
          double a3 = dl*(1.0-dm);		// yellow	: (x  , z+1)
          double a4 = dl*dm;			// red 		: (x-1, z+1)
          edep[thisx+1][thisz+1][beam] += a1*increment;	// blue
          edep[thisx][thisz+1][beam] += a2*increment;	// green
          edep[thisx+1][thisz+2][beam] += a3*increment;	// yellow
          edep[thisx][thisz+2][beam] += a4*increment;	// red
        } else if ( xp >= 0 && zp < 0 ){
          //printf("%s\n", "Blach");
          double dl = abs(zp);		// because zp < 0
          double dm = xp;
          double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
          double a2 = (1.0-dl)*dm;		// green	: (x+1, z  )
          double a3 = dl*(1.0-dm);		// yellow	: (x  , z-1)
          double a4 = dl*dm;			// red 		: (x+1, z-1)
          edep[thisx+1][thisz+1][beam] += a1*increment;	// blue
          edep[thisx+2][thisz+1][beam] += a2*increment;	// green
          edep[thisx+1][thisz][beam] += a3*increment;	// yellow
          edep[thisx+2][thisz][beam] += a4*increment;	// red
        } else if ( xp < 0 && zp < 0 ){
          double dl = abs(zp);		// because zp < 0
          double dm = abs(xp);		// because xp < 0
          double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
          double a2 = (1.0-dl)*dm;		// green	: (x-1, z  )
          double a3 = dl*(1.0-dm);		// yellow	: (x  , z-1)
          double a4 = dl*dm;			// red 		: (x-1, z-1)
          edep[thisx+1][thisz+1][beam] += a1*increment;	// blue
          edep[thisx][thisz+1][beam] += a2*increment;	// green
          edep[thisx+1][thisz][beam] += a3*increment;	// yellow
          edep[thisx][thisz][beam] += a4*increment;	// red
        } else {
          double store = edep[thisx][thisz][0];
          edep[thisx][thisz][0] = store + (nuei[i] * (*(eden[thisx]+thisz))/ncrit * uray[i-1]*dt);
          cout << "***** ERROR in interpolation of laser deposition to grid!! *****" << endl;
          break;
        }
        amplitude_norm[i] = (pow(omega,2.0)-pow(*(wpe[thisx_00]+thisz_00),2.0))/(pow(omega,2.0)-pow(pow(*(wpe[thisx]+thisz),2.0),(1./4.)));
        mytime[i] = dt*i; // Sets the time, equal to the time step multiplied by the number of time steps

        // This will cause the code to stop following the ray once it escapes the extent of the plasma:
        if ( (myx[i] < (xmin-(dx/2.0))) || (myx[i] > (xmax+(dx/2.0))))
        {          // the "|" means "or" (symbol above the return key)
          finalt = i-1;
          delete [] rayx;
          delete [] rayz;
          rayx = new double[finalt]{0.0};
          rayz = new double[finalt]{0.0};
          amp_norm = new double[finalt]{0.0};
          /*
          vector<double> plotx(myx, myx+finalt);
          vector<double> plotz(myz, myz+finalt);
          plt::plot(plotx, plotz, "ro");
          */
          for(int j = 0; j < finalt;j++)
          {
            amp_norm[j] = amplitude_norm[j];
            rayx[j] = myx[j];
            rayz[j] = myz[j];
          }
          delete [] amp_norm;
          break;                  // "breaks" out of the i loop once the if condition is satisfied
        } else if ( (myz[i] < (zmin-(dz/2.0))) || (myz[i] > (zmax+(dz/2.0)))){
             // the "|" means "or" (symbol above the return key)
          finalt = i-1;
          amp_norm = new double[finalt]{0.0};
          delete [] rayx;
          delete [] rayz;
          rayx = new double[finalt]{0.0};
          rayz = new double[finalt]{0.0};
          /*
          vector<double> plotx(myx, myx+finalt);
          vector<double> plotz(myz, myz+finalt);
          plt::plot(plotx, plotz, "ro");
          */
          for(int j = 0; j < finalt;j++)
          {
            amp_norm[j] = amplitude_norm[j];
            rayx[j] = myx[j];
            rayz[j] = myz[j];
          }
          delete [] amp_norm;
          break;
      }
  }
  delete [] myx;
  delete [] mytime;
  delete [] myz;
  delete [] mykx;
  delete [] mykz;
  delete [] myvx;
  delete [] nuei;
  delete [] myvz;
  delete [] amplitude_norm;
  delete [] markingx;
  delete [] markingz;
}






//use two threads here
void launch_ray_XZ(double x_init, double z_init, double kx_init, double kz_init)
{
  initializeArrXZ(x_init,z_init, kx_init, kz_init);
  rayLaunch();

}
