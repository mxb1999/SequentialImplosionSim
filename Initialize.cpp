#include "implSim.h"
#include "declarations.h"
using namespace std;
//dynamically allocate and initialize the arrays


void initialize()
{

  auto start = chrono::high_resolution_clock::now();

  //dynamic allocation (using pointers to access arrays so the stack is not filled)
  intersections = new double*[nx]; //nx nz
  marked = new int*[nbeams*numstored]; //numstored nbeams nx nz
  dedendx = new double*[nx]; //nx nz
  truemark = new int**[nbeams*nrays];
  dedendz = new double*[nx]; //nx nz
  x = new double[nx]; //nx nz
  z = new double[nz]; //nx nz
  eden = new double*[nx]; //nx nz
  edep = new double**[nx+2]; //nx+2 nz+2 nbeams
  present = new int**[nx]; //nx nz nbeams
  machnum = new double*[nx]; //nx nz
  boxes = new int***[nbeams]; //nbeams nrays ncrossings 2
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
  i_b = new double**[nbeams]; //nbeams nx nz
  wpe = new double*[nx]; //nx nz
  crossesz = new double**[nbeams]; //nbeams nrays ncrossings
  crossesx = new double**[nbeams]; //nbeams nrays ncrossings
  ints = new int**[nbeams]; //nbeams nrays ncrossings

  for(int i = 0; i < nbeams*nrays;i++)
  {
    truemark[i] = new int*[nx*nz];
    for(int j = 0; j < nx*nz;j++)
    {
      truemark[i][j] = new int[nrays]{0};
    }
  }
  auto check1 = chrono::high_resolution_clock::now();
  for(int i = 0; i < nx+2; i++)
  {
    edep[i] = new double*[nz+2];
      for(int j = 0; j < nz+2; j++)
      {
        edep[i][j] = new double[nbeams];
      }
  }
  auto interim = chrono::high_resolution_clock::now();
  for(int i = 0; i < numstored*nbeams; i++)
  {
    marked[i] = new int[nx*nz];
  }
  for(int i = 0; i < nx; i++)
  {
    intersections[i] = new double[nz];
    eden[i] = new double[nz];
    machnum[i] = new double[nz];
    dedendx[i] = new double[nz];
    dedendz[i] = new double[nz];
    present[i] = new int*[nz];
    //marked[i] = new int**[nz];
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
    for(int j = 0; j < nz; j++)
    {
        present[i][j] = new int[nbeams]{0};
        W1_storage[i][j] = new double[nrays];
        W2_storage[i][j] = new double[nrays];
      }
  }
  auto check2 = chrono::high_resolution_clock::now();
  auto tester = chrono::high_resolution_clock::now();
    for(int i = 0; i < nbeams; i++)
    {
      boxes[i] = new int**[nrays];
      dkx[i] = new double*[nrays];
      dkz[i] = new double*[nrays];
      dkmag[i] = new double*[nrays];
      crossesz[i] = new double*[nrays];
      crossesx[i] = new double*[nrays];
      ints[i] = new int*[nrays];
      i_b[i] = new double*[nx];
        for(int j = 0; j < nrays; j++)
        {
          i_b[i][j] = new double[nz]{0};
          dkx[i][j] = new double[ncrossings-1]{0};
          dkz[i][j] = new double[ncrossings-1]{0};
          dkmag[i][j] = new double[ncrossings-1]{0};
          crossesz[i][j] = new double[ncrossings]{0};
          crossesx[i][j] = new double[ncrossings]{0};
          boxes[i][j] = new int*[ncrossings]{0};
          ints[i][j] = new int[ncrossings]{0};
          #pragma omp parallel for num_threads(2)
            for(int m = 0; m < ncrossings;m++)
            {
              boxes[i][j][m] = new int[2]{0};
            }
        }
    }
  auto check3 = chrono::high_resolution_clock::now();

  cout << "Setting initial conditions for ray tracker..." <<endl;
  cout << "nrays per beam is"<< nrays <<endl;

  //Calculating the initial energy density, wpe, and machnum values
  span(x, xmin, xmax, nx);
  span(z, zmin, zmax, nz);
  #pragma omp parallel for num_threads(12)
    for(int i = 0; i < nx;i++)
    {
      for(int j = 0; j < nz;j++)
      {
        eden[i][j] = fmax(0.0,((0.3*ncrit-0.1*ncrit)/(xmax-xmin))*(x[i]-xmin)+(0.1*ncrit));
        wpe[i][j] = sqrt(eden[i][j]*1e6*pow(ec,2.0)/(me*e0));
        machnum[i][j] = fmax(0.0,(((-0.4)-(-2.4))/(xmax-xmin))*(x[i]-xmin))+(-2.4);
      }
    }
  #pragma omp parallel for num_threads(12)
    for(int i = 0; i < nx-1; i++)
    {
        for(int j = 0; j < nz-1; j++)
        {
          dedendx[i][j] = (eden[i+1][j]-eden[i][j])/(x[i+1]-x[i]);
          dedendz[i][j] = (eden[i][j+1]-eden[i][j])/(z[j+1]-z[j]);
        }
    }
  auto check4 = chrono::high_resolution_clock::now();
  int max = fmax(nx, nz);
  #pragma omp parallel for num_threads(12)
    for(int i = 0; i < max;i++)
    {
      if(i < nx)
      {
        dedendz[i][nz-1] = dedendz[i][nz-2];

      }
      if(i < nz)
      {
        dedendx[nx-1][i] = dedendx[nx-2][i];
      }
    }
  auto check5 = chrono::high_resolution_clock::now();
  cout << "Initialize CPU Time 1: " << chrono::duration_cast<chrono::milliseconds>(check1-start).count() << " ms" << endl;
  cout << "Test Time: " << chrono::duration_cast<chrono::milliseconds>(tester-check2).count() << " ms" << endl;
  cout << "Initialize CPU Time 2: " << chrono::duration_cast<chrono::milliseconds>(check2-interim).count() << " ms" << endl;
  cout << "Initialize CPU Time 3: " << chrono::duration_cast<chrono::milliseconds>(check3-check2).count() << " ms" << endl;
  cout << "Initialize CPU Time 4: " << chrono::duration_cast<chrono::milliseconds>(check4-check3).count() << " ms" << endl;
  cout << "Initialize CPU Time 5: " << chrono::duration_cast<chrono::milliseconds>(check4-check3).count() << " ms" << endl;

}
