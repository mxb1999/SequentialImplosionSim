#include "implSim.h"


using namespace std;
//dynamically allocate and initialize the arrays
void initialize()
{
  uray = new double[nt];
  for(int i = 0; i < nt; i++)
  {
    uray[i] = 1.0;
  }
  //dynamic allocation (using pointers to access arrays so the stack is not filled)
  intersections = new double*[nx]; //nx nz
  marked = new int***[nx]; //nx nz numstored nbeams
  dedendx = new double*[nx]; //nx nz
  dedendz = new double*[nx]; //nx nz
  x = new double[nx]{0.0}; //nx nz
  z = new double[nz]{0.0}; //nx nz
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
  i_b1 = new double*[nx];//nx+2 nz+2
  i_b2 = new double*[nx];//nx+2 nz+2
  wpe = new double*[nx]; //nx nz
  crossesz = new double**[nbeams]; //nbeams nrays ncrossings
  crossesx = new double**[nbeams]; //nbeams nrays ncrossings
  ints = new int**[nbeams]; //nbeams nrays ncrossings
  myx= new double[nt]{0.0};
  myz = new double[nt]{0.0};
  mykx=new double[nt]{0.0};
  mykz= new double[nt]{0.0};
  myvx=new double[nt]{0.0};
  myvz= new double[nt]{0.0};
  for(int i = 0; i < nx+2; i++)
  {
    if(i < nx)
    {
      intersections[i] = new double[nz]{0.0};
      eden[i] = new double[nz]{0.0};
      machnum[i] = new double[nz]{0.0};
      dedendx[i] = new double[nz]{0.0};
      dedendz[i] = new double[nz]{0.0};
      present[i] = new int*[nz];
      marked[i] = new int**[nz];
      W1_storage[i] = new double*[nz];
      W2_storage[i] = new double*[nz];
      u_flow[i] = new double[nz]{0.0};
      W1[i] = new double[nz]{0.0};
      W2[i] = new double[nz]{0.0};
      W1_init[i] = new double[nz]{0.0};
      W2_init[i] = new double[nz]{0.0};
      W1_new[i] = new double[nz]{0.0};
      W2_new[i] = new double[nz]{0.0};
      wpe[i] = new double[nz]{0.0};
      marked[i] = new int**[nz];
      i_b1[i] = new double[nz]{0.0};
      i_b2[i] = new double[nz]{0.0};
    }

    edep[i] = new double*[nz+2];
    for(int j = 0; j < nz+2; j++)
    {
      if(j < nz && i < nx)
      {
        marked[i][j] = new int*[numstored];
        present[i][j] = new int[nbeams]{0};
        W1_storage[i][j] = new double[numstored]{0.0};
        W2_storage[i][j] = new double[numstored]{0.0};
        for(int m = 0; m < numstored;m++)
        {
          marked[i][j][m] = new int[nbeams]{0};
        }
      }
      edep[i][j] = new double[nbeams]{0.0};

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
      dkx[i][j] = new double[ncrossings-1]{0.0};
      dkz[i][j] = new double[ncrossings-1]{0.0};
      dkmag[i][j] = new double[ncrossings-1]{0.0};
      crossesz[i][j] = new double[ncrossings]{0.0};
      crossesx[i][j] = new double[ncrossings]{0.0};
      boxes[i][j] = new int*[nx*3];
      ints[i][j] = new int[ncrossings]{0};
      for(int m = 0; m < nx*3;m++)
      {
        boxes[i][j][m] = new int[2]{0};
      }
    }
  }
  cout << "Setting initial conditions for ray tracker..." <<endl;
  cout << "nrays per beam is"<< nrays <<endl;

  //Calculating the initial energy density, wpe, and machnum values
  span(x, xmin, xmax, nx);
  span(z, zmin, zmax, nz);
  for(int i = 0; i < nx;i++)
  {
    for(int j = 0; j < nz;j++)
    {
      eden[i][j] = fmax(0.0,((0.3*ncrit-0.1*ncrit)/(xmax-xmin))*(x[i]-xmin)+(0.1*ncrit));
      wpe[i][j] = sqrt(eden[i][j]*1e6*pow(ec,2.0)/(me*e0));
      machnum[i][j] = fmax(0.0,(((-0.4)-(-2.4))/(xmax-xmin))*(x[i]-xmin))+(-2.4);
    }
  }

  for(int i = 0; i < nx-1; i++)
  {
    for(int j = 0; j < nz-1; j++)
    {
      dedendx[i][j] = (eden[i+1][j]-eden[i][j])/(x[i+1]-x[i]);
      dedendz[i][j] = (eden[i][j+1]-eden[i][j])/(z[j+1]-z[j]);
    }
  }
  for(int i = 0; i < fmax(nx,nz);i++)
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
  //plt::plot(xplot, zplot);
//  plt::show();
}
