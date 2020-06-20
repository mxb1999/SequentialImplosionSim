#include "implSim.h"

using namespace std;
//initializing arrays and constants for the CBET subroutines
void initializeArr()
{
  cs = 1e2*sqrt(ec*(Z*Te_eV+3.0*Ti_eV)/mi_kg);	// acoustic wave speed, approx. 4e7 cm/s in this example
    for(int i = 0; i < nx;i++)
    {
      for(int j = 0; j < nz;j++)
      {

        i_b1[i][j] = edep[i][j][0];
        i_b2[i][j] = edep[i][j][1];
        u_flow[i][j] = machnum[i][j]*cs;
        W1[i][j] = sqrt(1.0-eden[i][j]/ncrit)/double(rays_per_zone);
        W2[i][j] = sqrt(1.0-eden[i][j]/ncrit)/double(rays_per_zone);

        W1_new[i][j] = W1[i][j];
        W2_new[i][j] = W2[i][j];
        W1_init[i][j] = W1_new[i][j];
        W2_init[i][j] = W2_new[i][j];
      }
    }
    for(int i = 0; i < nbeams;i++)
    {
      for(int j = 0; j < nrays; j++)
      {
        for(int m = 0; m < ncrossings-1; m++)
        {
          dkx[i][j][m] = crossesx[i][j][m+1]-crossesx[i][j][m];
          dkz[i][j][m] = crossesz[i][j][m+1]-crossesz[i][j][m];
          dkmag[i][j][m] = sqrt(pow(dkx[i][j][m],2.0)+pow(dkz[i][j][m],2.0));
        }
      }
    }
}

void gain_CBET()
{
  double constant1 = (pow(estat,2.0))/(4*(1.0e3*me)*c*omega*kb*Te*(1+3*Ti/(Z*Te)));
  cout << "Calculating CBET gains" << endl;
  initializeArr();
  //performs the CBET calculation for each crossing of each ray calculated for the beams
  for(int i = 0; i < nbeams-1;i++)
  {
    for(int j = 0; j < nrays; j++)
    {
      for(int m = 0; m < ncrossings;m++)
      {
        //marked array has dimensions nx nz numstored nbeams
        int ix = boxes[i][j][m][0];
        int iz = boxes[i][j][m][1];
          if(intersections[ix][iz] != 0)
          {
            vector<int> nonzeros1;
            vector<int> nonzeros2;
            int numrays1 = 0;
            int numrays2 = 0;
            //finding the nonzero locations within marked
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
            //transferring the locations of the nonzero values to reference arrays
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
              int ix2 = boxes[i+1][mark2Copy[n]][q][0];
              int iz2 = boxes[i+1][mark2Copy[n]][q][1];
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
          double ne = eden[ix][iz];
          double epsilon = 1.0-ne/ncrit;
          double kmag = (omega/c)*sqrt(epsilon);         // magnitude of wavevector

          double kx1 = kmag * dkx[i][j][m] / (dkmag[i][j][m] + 1.0e-10);
          double kx2 = kmag * dkx[i+1][marker2[n2]][mark2Copy[n2]] / (dkmag[i+1][marker2[n2]][mark2Copy[n2]] + 1.0e-10);
          double kz1 = kmag * dkz[i][j][m]/(dkmag[i][j][m]+1.0e-10);
          double kz2 = kmag * dkz[i+1][marker2[n2]][mark2Copy[n2]] / (dkmag[i+1][marker2[n2]][mark2Copy[n2]] + 1.0e-10);


          double kiaw = sqrt(pow(kx2-kx1,2.0)+pow(kz2-kz1,2.0));
        // magnitude of the difference between the two vectors

          double ws = kiaw*cs;            // acoustic frequency, cs is a constant
          double omega1= omega;  // laser frequency difference. To start, just zero.
          double omega2 = omega;
          double eta = ((omega2-omega1)-(kx2-kx1)*u_flow[ix][iz])/(ws+1.0e-10);


         double efield1 = sqrt(8.*pi*1.0e7*i_b1[ix][iz]/c);             // initial electric field of rayi;lnlni46
        // double efield2 = sqrt(8.*pi*1.0e7*(*(i_b2[ix]+iz))/c);             // initial electric field of ray

         double P = (pow((iaw),2)*eta)/(pow((pow(eta,2)-1.0),2)+pow((iaw),2)*pow(eta,2));         // From Russ's paper
         //double gain1 = constant1*pow(efield2,2)*(ne/ncrit)*(1/iaw)*P;               //L^-1 from Russ's paper
         double gain2 = constant1*pow(efield1,2)*(ne/ncrit)*(1/iaw)*P;
                  //L^-1 from Russ's paper

                          // new energy of crossing (PROBE) ray (beam 2)
          //updating arrays based upon the calculated magnitude changes
          if (dkmag[i+1][marker2[n2]][mark2Copy[n2]] >= 1.0*dx)
          {
            W2_new[ix][iz] = W2[ix][iz]*exp(-1*W1[ix][iz]*dkmag[i+1][marker2[n2]][mark2Copy[n2]]*gain2/sqrt(epsilon));
            W2_storage[ix][iz][n2] = W2_new[ix][iz];
            W1_new[ix][iz] = W1[ix][iz]*exp(1*W2[ix][iz]*dkmag[i][j][m]*gain2/sqrt(epsilon));
            W1_storage[ix][iz][raynum] = W1_new[ix][iz];
          }

          }
        }
      }
      if ( j % 20 == 0 )
      {
        cout << "     ..."<<(int)(100.*(1.0-(double)(j/double(1*nrays))))<<"%  remaining..."<<endl;
      }
    }


  }
}
//updating CBET Intensities
void update_CBET()
{
  cout << "Updating CBET Intensities" << endl;
  i_b1_new = new double*[nx];
  i_b2_new = new double*[nx];
  for(int i = 0; i < nx; i++)
  {
    i_b1_new[i] = new double[nz];
    i_b2_new[i] = new double[nz];
    for(int j = 0; j < nz; j++)
    {
      i_b1_new[i][j] = i_b1[i][j];
      i_b2_new[i][j] = i_b2[i][j];
    }
  }
  //perform calculation for each crossing of each ray calculated for the beams in launch rays
  for(int i = 0; i < nbeams-1;i++)
  {
    for(int j = 0; j < nrays; j++)
    {
      for(int m = 0; m < ncrossings; m++)
      {
        int ix = boxes[i][j][m][0];
        int iz = boxes[i][j][m][1];
        //if two beams intersect
        if(intersections[ix][iz] != 0)
        {
          //find nonzero locations of rays
          vector<int> nonzero1;
          vector<int> nonzero2;
          int numrays1 = 0;
          int numrays2 = 0;
          for(int q = 0; q < numstored; q++)
          {
            if(marked[ix][iz][q][0] != 0)
            {
              nonzero1.push_back(q);
              numrays1++;
            }
            if(marked[ix][iz][q][1] != 0)
            {
              nonzero2.push_back(q);
              numrays2++;
            }
          }
          int marker1[numrays1]{0};
          int marker2[numrays2]{0};
          int ccopy[numrays2]{0};
          int rcopy[numrays2]{0};
          for(int q = 0; q < fmax(numrays1, numrays2);q++)
          {
            if(q < numrays1)
            {
              marker1[q] = marked[ix][iz][nonzero1[q]][0];
            }
            if(q < numrays2)
            {
              marker2[q] = marked[ix][iz][nonzero2[q]][1];
              //make copies of the marker array to be modified
              ccopy[q] = marker2[q];
              rcopy[q] = marker2[q];
            }
          }
          int ray1num = 0;
          for(int r = 0; r < numrays1; r++)
          {
            if(marker1[r] == j)
            {
              ray1num = r;
              break;
            }
          }

          for(int n = 0; n < numrays2; n++)
          {
            for(int q = 0; q < ncrossings;q++)
            {
              int ix2 = boxes[i+1][rcopy[n]][q][0];
              int iz2 = boxes[i+1][rcopy[n]][q][1];
              if(ix == ix2 && iz == iz2)
              {
                ccopy[n] = q;
                break;
              }
            }
          }
          //calculate the fractional change in CBET field to be applied
          double fractional_change_1 = (-1.0 * (1.0 - (W1_new[ix][iz]/W1_init[ix][iz])) * i_b1[ix][iz]);
          double fractional_change_2 = (-1.0 * (1.0 - (W2_new[ix][iz]/W2_init[ix][iz])) * i_b2[ix][iz]);
          /*cout <<fractional_change_1<< endl;
          cout << fractional_change_2 << endl;
          cout <<i_b1[ix][iz]<< endl;
          cout << i_b2[ix][iz] << endl;
          cout <<1.0 - W1_new[ix][iz]/W1_init[ix][iz]<< endl;
          cout << 1.0 - W2_new[ix][iz]/W2_init[ix][iz] << endl;
          cout << endl;
          */
          //apply the change to the field at the current xz position
          i_b1_new[ix][iz] += fractional_change_1;
          i_b2_new[ix][iz] += fractional_change_2;
          double x_prev_1 = x[ix];
          double z_prev_1 = z[iz];
          double x_prev_2 = x[ix];
          double z_prev_2 = z[iz];
          for(int l = m+1;l < ncrossings;l++)
          {
            int ix_next_1 = boxes[0][j][l][0];
            int iz_next_1 = boxes[0][j][l][1];
            double x_curr_1 = x[ix_next_1];
            double z_curr_1 = z[iz_next_1];

            if(ix_next_1 == 0 || iz_next_1 == 0)
            {
              break;
            }else
            {

              if(x_curr_1 != x_prev_1 || z_curr_1 != z_prev_1)
              {
                i_b1_new[ix_next_1][iz_next_1] += fractional_change_1 * (present[ix][iz][0]/present[ix_next_1][iz_next_1][0]);
              }
              x_prev_1 = x_curr_1;
              z_prev_1 = z_curr_1;
            }
          }
          int n2 = fmin(ray1num, numrays2 - 1);

          for(int q = ccopy[n2]+1;q < ncrossings; q++)
          {
            int ix_next_2 = boxes[1][rcopy[n2]][q][0];
            int iz_next_2 = boxes[1][rcopy[n2]][q][1];
            double x_curr_2 = x[ix_next_2];
            double z_curr_2 = z[iz_next_2];
            if(ix_next_2 == 0 || iz_next_2 == 0)
            {
              break;
            }else
            {

              if(x_curr_2 != x_prev_2 || z_curr_2 != z_prev_2)
              {
                i_b2_new[ix_next_2][iz_next_2] += fractional_change_2 * (present[ix][iz][0]/present[ix_next_2][iz_next_2][1]);
              }
              x_prev_2 = x_curr_2;
              z_prev_2 = z_curr_2;
            }

          }

        }
      }
      if ( j % 20 == 0 )
      {
        cout << "     ..."<<(int)(100.*(1.0-(double)(j/double(1*nrays))))<<"%  remaining..."<<endl;
      }
    }
  }

  i_bplot = new double*[nx];
  i_b_newplot = new double*[nx];
  edepplot = new double*[nx];
  edenplot = new double*[nx];
  for(int i = 0; i < nx;i++)
  {
    i_bplot[i] = new double[nz];
    i_b_newplot[i] = new double[nz];
    edepplot[i] = new double[nz]{0.0};
    edenplot[i] = new double[nz]{0.0};
    for(int j = 0; j < nz;j++)
    {
      edenplot[i][j] = eden[i][j]/ncrit;
      i_bplot[i][j] = 8.53e-10*sqrt(i_b1[i][j]+i_b2[i][j]+1.0e-10)*(1.053/3.0);
      i_b_newplot[i][j] = 8.53e-10*sqrt(fmax(1.0e-10,i_b1_new[i][j])+fmax(1.0e-10,i_b2_new[i][j]))*(1.053/3.0);


      for(int m = 0;m<nbeams;m++)
      {
        edepplot[i][j]+=edep[i][j][m];
      }
    }
  }

}
void cbet()
{
  gain_CBET();
  update_CBET();
  for(int i = 0; i < nx; i++)
  {
    for(int j = 0; j < nz; j++)
    {

      if(i_b1_new[i][j] != i_b1[i][j]){
        //  cout << i_b2_new[i][j] << endl;

       cout << "i: " << i << " j: " << j<< endl;
      //cout << i_b2[i][j] << endl;
  //  cout << endl;
      }
    }
  }
  for(int i = 0; i < nbeams;i++)
  {
    for(int j = 0; j < nrays; j++)
    {
      for(int m = 0; m < ncrossings; m++)
      {
        if(boxes[i][j][m][0] != 0 || boxes[i][j][m][1] != 0)
        {
        //  cout << boxes[i][j][m][0] << "  " <<boxes[i][j][m][1] << " || ";
        }
      }
    }
  }
  {

  }

}
