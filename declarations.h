#include "vector.h"


//Launch Ray Values
double cs;
int ray1num;
int beam;
double maxDev;
//Pointers for necessary arrays
double* rayx;//nt
int*** truemark;//nx nz nbeams nrays needs to store all of the ray intersections. Boxes was previously used in order to store the location, and marked a certain number of rays. To generalize the algorithm, another must be used.

double* rayz;//nt
double* amp_norm;//nt
double** intersections; //nx nz
int** marked; //nx nz numstored nbeams
double** dedendx; //nx nz
double** dedendz; //nx nz
double* x; //nx
double* z; //nz
double** eden; //nx nz

double*** edep; //nx+2 nz+2 nbeams
int*** present; //nx nz nbeams
double** machnum; //nx nz
int**** boxes; //nbeams nrays nx*3 2
double*** W1_storage; //nx nz numstored
double*** W2_storage; //nx nz numstored
double** u_flow; //nx nz
double*** dkx; //nbeams nrays 2
double*** dkz; //nbeams nrays 2
double*** dkmag; //nbeams nrays 2
double** W1;//nx nz
double** W2;//nx nz
double** W1_init;//nx nz
//Launch_Ray_XZ specific arrays (all have a length of nt)


//CBET specific arrays
double** W2_init;//nx nz
double** W1_new;//nx nz
double** W2_new;//nx n
double*** i_b; //nbeams nx nz
double*** i_b_new;//nbeams nx nz
double** wpe; //nx nz
double*** crossesz; //nbeams nrays ncrossings
double*** crossesx; //nbeams nrays ncrossings
int*** ints; //nbeams nrays ncrossings
//arrays used only for plotting
double** i_bplot;//nx nz
double** i_b_newplot;//nx nz
double** edenplot; //the array is eden/ncrit,  nx nz
double** edepplot; //nx nz
