#include "implSim.h"

using namespace std;
//initializing arrays and constants for the CBET subroutines
void initializeArr()
{
  cs = 1e2*sqrt(ec*(Z*Te_eV+3.0*Ti_eV)/mi_kg);	// acoustic wave speed, approx. 4e7 cm/s in this example
  //#pragma omp parallel for num_threads(threads)
    for(int i = 0; i < nx;i++)
    {
      for(int j = 0; j < nz;j++)
      {
        for(int m = 0; m < nbeams; m++)
        {
          i_b[m][i][j] =  edep[i][j][m];
          //cout << edep[i][j][m] << endl;
        }
        u_flow[i][j] = machnum[i][j]*cs;
        W1[i][j] = sqrt(1.0-eden[i][j]/ncrit)/double(rays_per_zone);
        W2[i][j] = sqrt(1.0-eden[i][j]/ncrit)/double(rays_per_zone);
        W1_new[i][j] = W1[i][j];
        W2_new[i][j] = W2[i][j];
        W1_init[i][j] = W1_new[i][j];
        W2_init[i][j] = W2_new[i][j];
      }
    }
  cout << scientific;
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
      }//nbeams nrays ncrossings
    }
}
void gain_CBETSeq()
{
  double constant1 = (pow(estat,2.0))/(4*(1.0e3*me)*c*omega*kb*Te*(1+3*Ti/(Z*Te)));
  cout << "Calculating CBET gains" << endl;
  initializeArr();

  //performs the CBET calculation for each crossing of each ray calculated for the beams
  int storeDup[nx][nz]{0};
  for(int i = 0; i < nbeams-1;i++)
  {
    //#pragma omp parallel for num_threads(threads)
    for(int j = 0; j < nrays; j++)
    {
      for(int m = 0; m < ncrossings;m++)
      {
        //marked array has dimensions nx nz numstored nbeams
        int ix = boxes[i][j][m][0] - 1;
        int iz = boxes[i][j][m][1] - 1;
          if((ix >= 0 && iz >= 0))
          {
            if(intersections[ix][iz] != 0)
            {

            storeDup[ix][iz]++;
            vector<int> nonzeros1;
            vector<int> nonzeros2;
            int numrays1 = 0;
            int numrays2 = 0;
            //finding the nonzero locations within marked
            for(int q = 0; q < nrays; q++)
            {
              if(truemark[0*nrays+j][ix*nz+iz][q] != 0)
              {
                nonzeros1.push_back(q);
                numrays1++;
              }
              if(truemark[1*nrays+j][ix*nz+iz][q] != 0)
              {
                nonzeros2.push_back(q);
                numrays2++;
              }
            }

            int mark2Copy[numrays2];
            //transferring the locations of the nonzero values to reference arrays
            for(int l = 0; l < numrays2; l++)
            {
                mark2Copy[l] = nonzeros2[l];
            }
          for(int r = 0; r < numrays1;r++)
          {
            if(nonzeros1[r] == j)
            {
              ray1num = r;
              break;
            }
          }
          for(int n = 0; n < numrays2; n++)
          {
            for(int q = 0; q <ncrossings; q++)
            {
              int ix2 = boxes[i+1][nonzeros2[n]][q][0]- 1;
              int iz2 = boxes[i+1][nonzeros2[n]][q][1]- 1;
              if ( ix == ix2 && iz == iz2 )
              {
                mark2Copy[n] = q;
                break;
              }
            }
          }

        int n2limit = fmin(present[ix][iz][0],numrays2);

        for ( int n2 = 0; n2 < n2limit; n2++)
        {
          double ne = eden[ix][iz];
          double epsilon = 1.0-ne/ncrit;
          double kmag = (omega/c)*sqrt(epsilon);         // magnitude of wavevector

          double kx1 = kmag * dkx[i][j][m] / (dkmag[i][j][m] + 1.0e-10);
          double kx2 = kmag * dkx[i+1][nonzeros2[n2]][mark2Copy[n2]] / (dkmag[i+1][nonzeros2[n2]][mark2Copy[n2]] + 1.0e-10);
          double kz1 = kmag * dkz[i][j][m]/(dkmag[i][j][m]+1.0e-10);
          double kz2 = kmag * dkz[i+1][nonzeros2[n2]][mark2Copy[n2]] / (dkmag[i+1][nonzeros2[n2]][mark2Copy[n2]] + 1.0e-10);

          double kiaw = sqrt(pow(kx2-kx1,2.0)+pow(kz2-kz1,2.0));
        // magnitude of the difference between the two vectors
        cout << nonzeros2[n2]<< endl;

          double ws = kiaw*cs;            // acoustic frequency, cs is a constant
          double omega1= omega;  // laser frequency difference. To start, just zero.
          double omega2 = omega;
          double eta = ((omega2-omega1)-(kx2-kx1)*u_flow[ix][iz])/(ws+1.0e-10);


         double efield1 = sqrt(8.*pi*1.0e7*i_b[0][ix][iz]/c);             // initial electric field of rayi;lnlni46
        // double efield2 = sqrt(8.*pi*1.0e7*(*(i_b[1][ix]+iz))/c);             // initial electric field of ray

         double P = (pow((iaw),2)*eta)/(pow((pow(eta,2)-1.0),2)+pow((iaw),2)*pow(eta,2));         // From Russ's paper
         //double gain1 = constant1*pow(efield2,2)*(ne/ncrit)*(1/iaw)*P;               //L^-1 from Russ's paper
         double gain2 = constant1*pow(efield1,2)*(ne/ncrit)*(1/iaw)*P;
                  //L^-1 from Russ's paper
                          // new energy of crossing (PROBE) ray (beam 2)
          //updating arrays based upon the calculated magnitude changes
          //Chronological Issue for Parallelization: W1_new and W2_new can be written to at ix and iz multiple times'

          if (dkmag[i+1][nonzeros2[n2]][mark2Copy[n2]] >= 1.0*dx)
          {
              #pragma omp atomic write
              W2_new[ix][iz] = W2[ix][iz]*exp(-1*W1[ix][iz]*dkmag[i+1][nonzeros2[n2]][mark2Copy[n2]]*gain2/sqrt(epsilon));
              #pragma omp atomic write
              W1_new[ix][iz] = W1[ix][iz]*exp(1*W2[ix][iz]*dkmag[i][j][m]*gain2/sqrt(epsilon));

              W2_storage[ix][iz][j] =W2_new[ix][iz];
          }

          }
        }}
      }
      if ( j % 20 == 0 )
      {
        //cout << "     ..."<<(int)(100.*(1.0-(double)(j/double(1*nrays))))<<"%  remaining..."<<endl;
      }
    }
  }
  /*
  for(int i = 0; i < nx; i++)
  {
    for(int j = 0; j < nz; j++)
    {
      for(int m = 0; m < nrays; m++)
      {
        if(W1_new[i][j] != W1[i][j])
        {
          W1_new[i][j] = W1_storage[i][j][m];
        }
        if(W2_new[i][j] != W2[i][j])
        {
          W2_new[i][j] = W2_storage[i][j][m];
        }
      }
    }
  }
  */
}
//updating CBET Intensities
double*** update_CBETSeq(double*** i_b_in)
{
  //cout << "Updating CBET Intensities" << endl;
  double*** i_b_temp = new double**[nbeams];
  for(int q = 0; q < nbeams; q++)
  {
    i_b_temp[q] = new double*[nx];
    for(int i = 0; i < nx; i++)
    {
      i_b_temp[q][i] = new double[nz];
    }
  }
  for(int q = 0; q < nbeams; q++)
  {
    for(int i = 0; i < nx; i++)
    {
      for(int j = 0; j < nz; j++)
      {
        i_b_temp[q][i][j] = i_b_in[q][i][j];
      }
    }
  }

  //perform calculation for each crossing of each ray calculated for the beams in launch rays
  for(int i = 0; i < nbeams-1;i++)
  {
    //#pragma omp parallel for num_threads(threads)
    for(int j = 0; j < nrays; j++)
    {
      for(int m = 0; m < ncrossings; m++)
      {
        if(boxes[i][j][m][0] == 0 || boxes[i][j][m][1] == 0)
        {
          break;
        }

        int ix = boxes[i][j][m][0]- 1;
        int iz = boxes[i][j][m][1]- 1;
        //if two beams intersect
        if(intersections[ix][iz] != 0)
        {
          //find nonzero locations of rays
          vector<int> nonzero1;
          vector<int> nonzero2;
          int numrays1 = 0;
          int numrays2 = 0;
          //#pragma omp parallel for num_threads(6)
          for(int q = 0; q < nrays; q++)
          {
            if(truemark[0*nrays+j][ix*nz+iz][q] != 0)
            {
              nonzero1.push_back(q);
              numrays1++;
            }
            if(truemark[1*nrays+j][ix*nz+iz][q] != 0)
            {
              nonzero2.push_back(q);
              numrays2++;
            }
          }
          int ccopy[numrays2]{0};
          int rcopy[numrays2]{0};
          int max = fmax(numrays1, numrays2);
        //  #pragma omp parallel for num_threads(6)
          for(int q = 0; q < max;q++)
          {
            if(q < numrays1)
            {
              nonzero1[q] = q;
            }
            if(q < numrays2)
            {
              nonzero2[q] = q;
              //make copies of the marker array to be modified
              ccopy[q] = nonzero2[q];
              rcopy[q] = nonzero2[q];
            }
          }
          int ray1num = 0;
          for(int r = 0; r < numrays1; r++)
          {
            if(nonzero1[r] == j)
            {
              ray1num = r;
              break;
            }
          }
          for(int n = 0; n < numrays2; n++)
          {
            for(int q = 0; q < ncrossings;q++)
            {
              int ix2 = boxes[i+1][rcopy[n]][q][0]- 1;
              int iz2 = boxes[i+1][rcopy[n]][q][1]- 1;
              if(ix == ix2 && iz == iz2)
              {
                ccopy[n] = q;
                break;
              }
            }
          }
          //calculate the fractional change in CBET field to be applied
          double fractional_change[2] = {(-1.0 * (1.0 - (W1_new[ix][iz]/W1_init[ix][iz])) * i_b_in[0][ix][iz]),
              (-1.0 * (1.0 - (W2_new[ix][iz]/W2_init[ix][iz])) * i_b_in[1][ix][iz])};

          for(int j = 0; j < nbeams; j++)
          {
            #pragma omp atomic update
            i_b_temp[j][ix][iz] += fractional_change[j];
          }
          double x_prev_1 = x[ix];
          double z_prev_1 = z[iz];
          double x_prev_2 = x[ix];
          double z_prev_2 = z[iz];
          for(int l = m+1;l < ncrossings;l++)
          {
            /*
            if(i_b_in[0][ix][iz] != 0 && fractional_change[0] / i_b_in[0][ix][iz] < convergence)
            {
              cout <<"check 3 1" << endl;
              break;
            }else if(i_b_in[0][ix][iz] != 0)
            {
              cout <<"check 3 2" << endl;
                maxDev = fmax(maxDev, fractional_change[0] / i_b_in[0][ix][iz]);
            }
            */
            int ix_next_1 = boxes[0][j][l][0]- 1;
            int iz_next_1 = boxes[0][j][l][1]- 1;


              if(ix_next_1 == -1 || iz_next_1 == -1)
              {
                break;
              }else
              {
                double x_curr_1 = x[ix_next_1];
                double z_curr_1 = z[iz_next_1];
                if(x_curr_1 != x_prev_1 || z_curr_1 != z_prev_1)
                {
                    #pragma omp atomic update
                    i_b_temp[0][ix_next_1][iz_next_1] += fractional_change[0] * ((double)present[ix][iz][0])/present[ix_next_1][iz_next_1][0];
                }
                x_prev_1 = x_curr_1;
                z_prev_1 = z_curr_1;
              }

          }
          int n2 = fmin(ray1num, numrays2 - 1);
          for(int q = ccopy[n2]+1;q < ncrossings; q++)
          {
            /*
            if(i_b_in[1][ix][iz] != 0 && fractional_change[1] / i_b_in[1][ix][iz] < convergence)
            {
              break;
            }else if(i_b_in[1][ix][iz] != 0){
                maxDev = fmax(maxDev, fractional_change[1] / i_b_in[1][ix][iz]);
            }
            */
            int ix_next_2 = boxes[1][rcopy[n2]][q][0]- 1;
            int iz_next_2 = boxes[1][rcopy[n2]][q][1]- 1;
            if(ix_next_2 == -1 || iz_next_2 == -1)
            {
              break;
            }else
            {
              double x_curr_2 = x[ix_next_2];
              double z_curr_2 = z[iz_next_2];
              if(x_curr_2 != x_prev_2 || z_curr_2 != z_prev_2)
              {
                  #pragma omp atomic update
                  i_b_temp[1][ix_next_2][iz_next_2] += fractional_change[1] * ((double)present[ix][iz][0])/(present[ix_next_2][iz_next_2][1]);
              }
              x_prev_2 = x_curr_2;
              z_prev_2 = z_curr_2;
            }
            //getchar();
          }

        }
      }

    }
  }

  return i_b_temp;
}

void cbet()
{
/*
  int runs = 1000;
  double sdev1[runs];
  double sdev2[runs];
  int cnt1[runs];
  int cnt2[runs];
  for(int m = 0; m < runs; m++)
  {
    int count1 = 0;
    int count2 = 0;
    gain_CBETSeq();
    update_CBETSeq();
    vector<double> store1;
    vector<double> store2;
    double ib1seq[nx][nz]{0};
    double ib2seq[nx][nz]{0};
    for(int i = 0; i < nx;i++)
    {
      for(int j = 0; j < nx;j++)
      {
        ib1seq[i][j] = i_b1_new[i][j];
        ib2seq[i][j] = i_b2_new[i][j];
      }
    }
    gain_CBETPar();
    update_CBETPar();
    for(int i = 0; i < nx;i++)
    {
      for(int j = 0; j < nz;j++)
      {
        if(i_b1_new[i][j] != ib1seq[i][j])
        {
          store1.push_back(i_b1_new[i][j]-ib1seq[i][j]);
          count1++;
        }
        if(i_b2_new[i][j] != ib2seq[i][j])
        {
          store2.push_back(i_b2_new[i][j]-ib2seq[i][j]);
          count2++;
        }
      }
    }
    sdev1[m] = stdDev(store1, count1);
    sdev2[m] = stdDev(store2, count2);
    cnt1[m] = count1;
    cnt2[m] = count2;
    cout << "______" << m << "______" << endl;
  }
  double devsum1 = 0;
  double devsum2 = 0;
  int cntsum1 = 0;
  int cntsum2 = 0;
  for(int i = 0; i < runs; i++)
  {
    devsum1 += sdev1[i];
    devsum2 += sdev2[i];
    cntsum1 += cnt1[i];
    cntsum2 += cnt2[i];
  }
  double avgDev1 = devsum1/runs;
  double avgDev2 = devsum2/runs;
  int avgcnt1 = cntsum1/runs;
  int avgcnt2 = cntsum2/runs;
  double avgpct1 = double(avgcnt1)/(nx*nz);
  double avgpct2 = double(avgcnt2)/(nx*nz);
  cout << "______________________________________________________________________" << endl;
  cout << "Average standard deviation for i_b1_newPar - i_b1_newSeq: " << avgDev1 << " for " << runs << " runs" << endl;
  cout << "Average standard deviation for i_b2_newPar - i_b2_newSeq: " << avgDev2 << " for " << runs << " runs" << endl;
  cout << "Percent of entries wrong for i_b1_newPar: " << avgpct1 << " for " << runs << " runs" << endl;
  cout << "Percent of entries wrong for i_b2_newPar: " << avgpct2 << " for " << runs << " runs" << endl;
*/


    auto start = chrono::high_resolution_clock::now();
    gain_CBETSeq();
    auto inter = chrono::high_resolution_clock::now();
    double prev = 1;
    maxDev = 0;
    double*** i_pass;
    i_pass = update_CBETSeq(i_b);

   while(abs(maxDev - prev) > convergence)
    {
      prev = maxDev;
      i_pass = update_CBETSeq(i_pass);
   }

    auto stop = chrono::high_resolution_clock::now();
    cout << "Gain CBET: " << chrono::duration_cast<chrono::milliseconds>(inter-start).count() << " ms" << endl;
    cout << "Update CBET: " << chrono::duration_cast<chrono::milliseconds>(stop-inter).count() << " ms" << endl;
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
        #pragma omp parallel for num_threads(threads)
        for(int j = 0; j < nz;j++)
        {
          edenplot[i][j] = eden[i][j]/ncrit;
          i_bplot[i][j] = 8.53e-10*sqrt(i_b[0][i][j]+i_b[1][i][j]+1.0e-10)*(1.053/3.0);
          i_b_newplot[i][j] = 8.53e-10*sqrt(fmax(1.0e-10,i_pass[0][i][j])+fmax(1.0e-10,i_pass[1][i][j]))*(1.053/3.0);
          for(int m = 0;m<nbeams;m++)
          {
            edepplot[i][j]+=edep[i][j][m];
          }
        }
      }
      /*
      for(int i = 0; i < nx; i++)
  {
    for(int j = 0; j < nz; j++)
    {
	    if(i_pass[1][i][j] != 0)
	    {
	      cout << i_pass[1][i][j] << endl;
	    }
    }
  }
  */
}
