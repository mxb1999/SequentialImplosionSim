#include "implSim.h"
using namespace std;
//initializing necessary arrays for the calculation
void initializeArrXZ(double x_init, double z_init, double kx_init, double kz_init)
{
  //Launch_Ray_XZ Array Declaration
  thisx_0 = 0;
  thisz_0 = 0;
  mytime = new double[nt]{0.0};
  span(mytime, dt, nt*dt, nt);
  amplitude_norm= new double[nt]{0.0};
  markingx = new double[nt]{0.0};
  markingz = new double[nt]{0.0};
  nuei = new double[nt]{0.0};
  //Initializing Arrays
  for(int i = 0; i < nt; i++)
  {
    nuei[i] = 1.0;
  }
  myx[0] = x_init;

  myz[0] = z_init;
  cout << scientific;
  //determining the initial x values within the desired range to track the beam
  for(int i = 0;i < nx;i++)
  {
    if(myx[0] - x[i] <= ((0.5+1.0e-10)*dx + 1e-11) && myx[0] - x[i] >= -1*((0.5+1.0e-10)*dx  + 1e-11) )
  {
      thisx_0=i;
      thisx_00=i;
      break;
    }
  }
  //determining the initial z values within the desired range to track the beam

  for(int i = 0;i < nz;i++)
  {
    if(beam == 0)
    {
    //  cout <<  "myz[0]: " << myz[0] << " z[i]: " << z[i] <<"myz[0] - z[i]: " << myz[0] - z[i] << ": ((0.5+1.0e-10)*dz): " << ((0.5+1.0e-10)*dz) << endl;
    }
    if(myz[0] - z[i] <= ((0.5+1.0e-10)*dz + 1e-11) && myz[0] - z[i] >= -1*((0.5+1.0e-10)*dz + 1e-11) )
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

}
//Launching the rays
void rayLaunch()
{
  //__________Time Stepping__________
    int numcrossing = 0;
    //looping through time intervals
    for(int i = 1; i < nt;i++)
    {
      myvz[i] = myvz[i-1] - pow(c,2.0)/(2.0*ncrit)*dedendz[thisx_0][thisz_0]*dt;
      myvx[i] = myvx[i-1] - pow(c,2.0)/(2.0*ncrit)*dedendx[thisx_0][thisz_0]*dt;
      myx[i] = myx[i-1] + myvx[i]*dt;
      myz[i] = myz[i-1] + myvz[i]*dt;

      int search_index_x = 1;
      int search_index_z = 1;
      int thisx_m = fmax(0, thisx_0-search_index_x);
      int thisx_p = fmin(nx-1, thisx_0+search_index_x);
      int thisz_m = fmax(0, thisz_0-search_index_z);
      int thisz_p = fmin(nz-1, thisz_0+search_index_z);

      //assigning the current x value to be tracked

      for(int j = thisx_m; j <= thisx_p;j++)
      {

        if ( myx[i] - x[j] <= ((0.5+1.0e-10)*dx + 1e-11) && myx[i] - x[j] >= -1*((0.5+1.0e-10)*dx + 1e-11))
        {
          //cout << "thisx: " << thisx << " for beam = " << beam << endl;
          thisx = j;
          if(beam == 0)
          {
          //  cout << thisx << endl;
          }
          break;
        }
      }

      //assigning the current z value to be tracked
      for(int j = thisz_m; j <= thisz_p; j++)
      {
        if(i == 1 && beam == 1)
        {
          //cout<< " myz[i]: "<<  myz[i]<< " z[j]: " <<z[j]<< " (0.5+1.0e-10)*dz: "<< (0.5+1.0e-10)*dz << endl;
        }
        if (myz[i] - z[j] <= ((0.5+1.0e-10)*dz + 1e-11) && myz[i] - z[j] >= -1*((0.5+1.0e-10)*dz + 1e-11))
        {

          thisz = j;
          break;
         }


      }

      double linez[2]={myz[i-1], myz[i]};
      double linex[2]={myx[i-1], myx[i]};
      int lastx = 10000;
      int lastz = 10000;
      //iterating through the selected portions of the x spatial tracking arrays
      for(int j = thisx_m; j <= thisx_p;j++)
      {
        double currx = x[j]-dx/2;
        if((myx[i] > currx && myx[i-1] <= currx) || (myx[i] < currx && myx[i-1] >= currx))
        {
          double crossx =interp(linez, linex, currx, 2);

          if(abs(crossx-lastz)>1.0e-20)
          {

            ints[beam][raynum][numcrossing]=uray[i];
            crossesx[beam][raynum][numcrossing] = currx;
            crossesz[beam][raynum][numcrossing] = crossx;
            if(myx[i] <= xmax+dx/2 && myx[i] >= xmin-dx/2)
            {
            //  cout << myx[i] << endl;
            //  cout << "currx: " << currx << endl;
              boxes[beam][raynum][numcrossing][0] = thisx;
              boxes[beam][raynum][numcrossing][1] = thisz;
            }
            lastx = currx;
            numcrossing += 1;
            break;
          }
        }
      }
      if(beam == 1 )
      {

      //  cout << thisx << endl;
    //    cout << myz[i] << endl;
      //  cout << myz[i-1] << endl;
        //cout << currz << endl;
        //cout << endl;
      }
    //  if(beam == 0){
    //    cout << "thisx_m: " << thisx_m<< " thisx_p: " << thisx_p << endl;
      //  cout << "thisz_m: " << thisz_m<< " thisz_p: " << thisz_p << endl;}

      //iterating through the selected portions of the z spatial tracking arrays
        for(int j = thisz_m; j <= thisz_p;j++)
        {
          double currz = z[j]-dz/2;


          if((myz[i] > (currz) && myz[i-1] < (currz)) || (myz[i] < (currz+1e-11) && myz[i-1] > (currz)))
          {
             double crossz = interp(linex, linez, currz, 2);

            if(abs(crossz-lastx)>1.0e-20)
            {
              ints[beam][raynum][numcrossing]=uray[i];
              crossesz[beam][raynum][numcrossing] = currz;
              crossesz[beam][raynum][numcrossing] = crossz;

              if(myz[i] < (zmax+dz/2 +1e-11) && myz[i] > (zmin-dz/2+1e-11))
              {
                /*if(beam == 0)
                {
                  cout << "beam: " <<  beam << "  raynum: " << raynum << "  numcrossing: " << numcrossing << endl;
                  cout << "thisx: " << thisx << " thisz: " << thisz << endl;
                }
                */
                cout << "[" << thisx << "," << thisz << "]" << endl;
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
        //ensuring that the beams are not identical
        if(markingx[i] != markingx[i-1] && markingz[i] != markingz[i-1])
        {
          ztarg = z[thisz] - (dz/2.0);
          if(myvz[i] < 0.0)
          {
            ztarg = z[thisz] + (dz/2.0);
          }
          slope = (myz[i] - myz[i-1])/(myx[i] - myx[i-1]+1.0e-10);

          xtarg = x[thisx]-(dx/2.0);
          if(myvx[i] >= 0.0)
          {
            xtarg = x[thisx]+(dx/2.0);
          }
          slope = (myx[i] - myx[i-1])/(myz[i] - myz[i-1]+1.0e-10);
          ztarg = myz[i-1]+(xtarg-myx[i-1])/slope;

          for(int j = 0; j < numstored;j++)
          {
            if(marked[thisx][thisz][j][beam] == 0)
            {

              marked[thisx][thisz][j][beam] = raynum;
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
          for(int j = 0; j < numstored;j++)
          {
            if(marked[thisx][thisz][j][beam] == 0)
            {
              marked[thisx][thisz][j][beam] = raynum;
              present[thisx][thisz][beam] += 1.0;
              break;
            }
          }
        }else if (markingx[i] != markingx[i-1])
        {
          xtarg = x[thisx]-dx/2.0;
          if(myvx[i] < 0.0)
          {
            xtarg = x[thisx]+dx/2.0;
          }
          slope = (myx[i] - myx[i-1])/(myz[i] - myz[i-1]+1.0e-10);
          ztarg = myz[i]+(xtarg-myx[i-1])/slope;

          for(int j = 0; j < numstored;j++)
          {
            if(marked[thisx][thisz][j][beam] == 0)
            {
              marked[thisx][thisz][j][beam] = raynum;
              present[thisx][thisz][beam] += 1.0;
              break;
            }
          }
        }



        //altering energy deposition values due to the launched rays
        uray[i] = uray[i-1];
  	    double increment = uray[i];
        double xp = (myx[i] - (x[thisx]+dx/2.0))/dx;
        double zp = (myz[i] - (z[thisz]+dz/2.0))/dz;



          //  cout << "increment: " << increment << endl;
          //  cout << "xp: " << xp<< "  zp: " << zp << endl;
      if(i == 1)
         {
           //cout << "thisx + 1: " << thisx + 1<< "  thisz + 1: " << thisz + 1 << endl;
        //   cout << "thisx - 1 + 1: " << thisx - 1 + 1 << " thisz - 1 + 1: " << thisz-1 + 1 << endl;
        //   cout << "thisx + 1 + 1: " << thisx + 1 + 1 << " thisz + 1 + 1: " << thisz+1 + 1 << endl;
        //   cout << endl;
      }

      //  cout << endl;

        if ( xp >= 0 && zp >= 0 ){
          //cout << "1" << endl;
        double dl = zp;
        double dm = xp;
        double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
        double a2 = (1.0-dl)*dm;		// green	: (x+1, z  )
        double a3 = dl*(1.0-dm);		// yellow	: (x  , z+1)
        double a4 = dl*dm;			// red 		: (x+1, z+1)
        edep[thisx+1][thisz+1][beam] += a1*increment;	// blue
        edep[thisx+1+1][thisz+1][beam] += a2*increment;	// green
        edep[thisx+1][thisz+1+1][beam] += a3*increment;	// yellow
        edep[thisx+1+1][thisz+1+1][beam] += a4*increment;	// red
      } else if ( xp < 0 && zp >= 0 ){
        double dl = zp;
      // cout << "2" << endl;
        double dm = abs(xp);		// because xp < 0
        double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
        double a2 = (1.0-dl)*dm;		// green	: (x-1, z  )
        double a3 = dl*(1.0-dm);		// yellow	: (x  , z+1)
        double a4 = dl*dm;			// red 		: (x-1, z+1)
        /*if(thisx == 0)
        {
          edep[nx-1][thisz][beam] += a2*increment;	// green
          edep[nx-1][thisz+1][beam] += a4*increment;
        }else
        {
          edep[thisx-1][thisz][beam] += a2*increment;	// green
          edep[thisx-1][thisz+1][beam] += a4*increment;	// red
        }*/
        edep[thisx+1][thisz+1][beam] += a1*increment;	// blue
        edep[thisx-1+1][thisz+1][beam] += a2*increment;	// green
        edep[thisx-1+1][thisz+1+1][beam] += a4*increment;	// red
        edep[thisx+1][thisz+1+1][beam] += a3*increment;	// yellow

      } else if ( xp >= 0 && zp < 0 ){
        //printf("%s\n", "Blach");
    //   cout << "3" << endl;
        double dl = abs(zp);		// because zp < 0
        double dm = xp;
        double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
        double a2 = (1.0-dl)*dm;		// green	: (x+1, z  )
        double a3 = dl*(1.0-dm);		// yellow	: (x  , z-1)
        double a4 = dl*dm;			// red 		: (x+1, z-1)
        /*if(thisz == 0)
        {
          edep[thisx][nz-1][beam] += a3*increment;	// yellow
          edep[thisx+1][nz-1][beam] += a4*increment;	// red
        }else{
          edep[thisx][thisz-1][beam] += a3*increment;	// yellow
          edep[thisx+1][thisz-1][beam] += a4*increment;	// red
        }*/
        edep[thisx+1][thisz-1+1][beam] += a3*increment;	// yellow
        edep[thisx+1+1][thisz-1+1][beam] += a4*increment;	// red
        edep[thisx+1][thisz+1][beam] += a1*increment;	// blue
        edep[thisx+1+1][thisz+1][beam] += a2*increment;	// green
      } else if ( xp < 0 && zp < 0 ){
//        cout << "4" << endl;
        double dl = abs(zp);		// because zp < 0
        double dm = abs(xp);		// because xp < 0
        double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
        double a2 = (1.0-dl)*dm;		// green	: (x-1, z  )
        double a3 = dl*(1.0-dm);		// yellow	: (x  , z-1)
        double a4 = dl*dm;			// red 		: (x-1, z-1)
      /*  if(thisx == 0 && thisz == 0)
        {

          edep[thisx+1][nz+1-1][beam] += a3*increment;	// yellow
          edep[nx+1-1][thisz+1][beam] += a2*increment;	// green
          edep[nx-1+1][nz-1+1][beam] += a4*increment;	// red
        }else if(thisx == 0)
        {
          edep[nx-1][thisz][beam] += a2*increment;	// green
          edep[nx-1][thisz-1][beam] += a4*increment;	// red
        }else if(thisz == 0)
        {
          edep[thisx][nz-1][beam] += a3*increment;	// yellow
          edep[thisx-1][nz-1][beam] += a4*increment;	// red
        }else
        {
        edep[thisx][thisz-1][beam] += a3*increment;	// yellow
        edep[thisx-1][thisz][beam] += a2*increment;	// green
        edep[thisx-1][thisz-1][beam] += a4*increment;	// red
        }*/
        edep[thisx+1][thisz+1][beam] += a1*increment;	// blue
        edep[thisx+1][thisz+1-1][beam] += a3*increment;	// yellow
        edep[thisx+1-1][thisz+1][beam] += a2*increment;	// green
        edep[thisx+1-1][thisz+1-1][beam] += a4*increment;	// red
      } else {
        double store = edep[thisx][thisz][0];
        edep[thisx][thisz][0] = store + (nuei[i] * (*(eden[thisx]+thisz))/ncrit * uray[i-1]*dt);
        cout << "***** ERROR in interpolation of laser deposition to grid!! *****" << endl;
        break;
      }
      amplitude_norm[i] = (pow(omega,2.0)-pow(*(wpe[thisx_00]+thisz_00),2.0))/(pow(omega,2.0)-pow(pow(*(wpe[thisx]+thisz),2.0),(1./4.)));
      mytime[i] = dt*i;
      //checking if the ray is "out of bounds"
      /*
      if(i == 1133)
      {
        cout << endl;
        for(int j = 0; j < nt; j++)
        {
          cout << myx[j] << " :: " << j << endl;
        }
        cout << endl;
      }
*/
      if ( (myx[i] < (xmin-(dx/2.0))) || (myx[i] > (xmax+(dx/2.0))))
      {
        finalt = i-1;
        delete [] rayx;
        delete [] rayz;
        rayx = new double[finalt]{0.0};
        rayz = new double[finalt]{0.0};
        amp_norm = new double[finalt]{0.0};
        for(int j = 0; j < finalt;j++)
        {
          amp_norm[j] = amplitude_norm[j];
          rayx[j] = myx[j];
          rayz[j] = myz[j];
        }
        delete [] amp_norm;
      /*  if(beam == 0)
        {cout << "X Break at: " << i << endl;
        cout << "myx[i]: " << myx[i] << endl;
        cout << "xmin: " << xmin << endl;
        cout << "xmax: " << xmax << endl;
        cout << endl;}*/
        break;                  // "breaks" out of the i loop once the if condition is satisfied
      } else if ( (myz[i] < (zmin-(dz/2.0))) || (myz[i] > (zmax+(dz/2.0)))){
           // the "|" means "or" (symbol above the return key)
        finalt = i-1;
        amp_norm = new double[finalt]{0.0};
        delete [] rayx;
        delete [] rayz;
        rayx = new double[finalt]{0.0};
        rayz = new double[finalt]{0.0};
        for(int j = 0; j < finalt;j++)
        {
          amp_norm[j] = amplitude_norm[j];
          rayx[j] = myx[j];
          rayz[j] = myz[j];
        }
        delete [] amp_norm;
      /*  if(beam == 0){
          cout << "Z Break at: " << i << endl;
          cout << "myz[i]: " << myz[i] << endl;
          cout << "zmin: " << zmin << endl;
          cout << "zmax: " << zmax << endl;
          cout << endl;}*/

        break;
    }
  }
  delete [] mytime;
  delete [] nuei;
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
